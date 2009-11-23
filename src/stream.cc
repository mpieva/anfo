//    Copyright 2009 Udo Stenzel
//    This file is part of ANFO
//
//    ANFO is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    Anfo is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with Anfo.  If not, see <http://www.gnu.org/licenses/>.

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include "compress_stream.h"
#include "stream.h"
#include "util.h"

extern "C" {
#include "fastlz.h"
}

#include <google/protobuf/repeated_field.h>

#include <set>
#include <iostream>

#if HAVE_FCNTL_H
#include <fcntl.h>
#endif

namespace std {
	bool operator < ( const output::Read& p, const output::Read& q )
	{ return p.SerializeAsString() < q.SerializeAsString() ; }

	bool operator < ( const config::Policy& p, const config::Policy& q )
	{ return p.SerializeAsString() < q.SerializeAsString() ; }
}

namespace streams {

using namespace google::protobuf::io ;
using namespace output ;
using namespace std ;

void transfer( Stream& in, Stream& out ) 
{
	out.put_header( in.fetch_header() ) ;
	while( in.get_state() == Stream::have_output && out.get_state() == Stream::need_input )
		out.put_result( in.fetch_result() ) ;
	out.put_footer( in.fetch_footer() ) ;
}

int anfo_reader__num_files_ = 0 ;

namespace {
	string basename( const string& s )
	{
		string::size_type p = s.rfind( '/' ) ;
		return p != string::npos ? s.substr(p+1) : s ;
	}
} ;

FastqReader::FastqReader( std::auto_ptr< google::protobuf::io::ZeroCopyInputStream > is, bool solexa_scores, char origin ) 
	: is_( is ), sol_scores_(solexa_scores), origin_(origin) { read_next_message() ; }

AnfoReader::AnfoReader( std::auto_ptr< google::protobuf::io::ZeroCopyInputStream > is, const std::string& name ) : is_( is ), name_( name )
{
	std::string magic ;
	CodedInputStream cis( is_.get() ) ;

	if( !cis.ReadString( &magic, 4 ) || magic != "ANFO" ) {
		throw ParseError( name_ + " is not an ANFO file" ) ;
	} else {
		uint32_t tag ;
		if( cis.ReadVarint32( &tag ) && tag == mk_msg_tag( 1 ) && cis.ReadVarint32( &tag ) ) {
			int lim = cis.PushLimit( tag ) ;
			if( hdr_.ParseFromCodedStream( &cis ) ) {
				sanitize( hdr_ ) ;
				cis.PopLimit( lim ) ;
				++anfo_reader__num_files_ ;
				read_next_message( cis, name_ ) ;
				return ;
			}
		}
		throw ParseError( "deserialization error in header of " + name_ ) ;
	}
	foot_.set_exit_code(1) ;
}

Result AnfoReader::fetch_result()
{
	Result r ;
	swap( r, res_ ) ;
	CodedInputStream cis( is_.get() ) ;
	read_next_message( cis, name_ ) ;
	return r ;
}

namespace {
	void upgrade_cigar( google::protobuf::RepeatedField<unsigned int>& n, const string& o )
	{
		for( unsigned i = 0 ; i != o.size() ; ++i )
			if(      (uint8_t)o[i] < 128 ) n.Add( mk_cigar( Hit::Match,  (unsigned)(uint8_t)o[i]       ) ) ;
			else if( (uint8_t)o[i] < 192 ) n.Add( mk_cigar( Hit::Insert, (unsigned)(uint8_t)o[i] - 128 ) ) ;
			else                           n.Add( mk_cigar( Hit::Delete, (unsigned)(uint8_t)o[i] - 192 ) ) ;
	}

	Hit upgrade( const OldHit& o ) 
	{
		Hit h ;
		h.set_genome_name( o.has_genome_name() ? o.genome_name() : string() ) ;
		h.set_sequence( o.sequence() ) ;
		h.set_start_pos( o.start_pos() ) ;
		h.set_aln_length( o.aln_length() ) ;
		h.set_score( (int)( 0.5 + o.score() / log(10.0) ) ) ;
		if( o.has_evalue() ) h.set_evalue( o.evalue() ) ;
		if( o.has_taxid() ) h.set_taxid( o.taxid() ) ;
		if( o.has_taxid_species() ) h.set_taxid_species( o.taxid_species() ) ;
		if( o.has_taxid_order() ) h.set_taxid_order( o.taxid_order() ) ;
		// if( o.has_genome_file() ) h.set_genome_file( o.genome_file() ) ;
		upgrade_cigar( *h.mutable_cigar(), o.cigar() ) ;
		return h ;
	}

	Result upgrade( const OldResult& o )
	{
		Result rs ;
		{
			Read &r = *rs.mutable_read() ;
			r.set_seqid( o.has_seqid() ? o.seqid() : string() ) ;
			if( o.has_description() ) r.set_description( o.description() ) ;
			r.set_sequence( o.sequence() ) ;
			if( o.has_quality() ) r.set_quality( o.quality() ) ;
			if( o.has_trim_left() ) r.set_trim_left( o.trim_left() ) ;
			if( o.has_trim_right() ) r.set_trim_right( o.trim_right() ) ;
		}

		rs.mutable_member()->MergeFrom( o.member() ) ;

		if( o.has_best_hit() ) *rs.add_hit() = upgrade( o.best_hit() ) ;
		if( o.has_best_to_genome() ) {
			Hit &h = *rs.add_hit() ;
			h = upgrade( o.best_to_genome() ) ;
			if( o.has_diff_to_next() ) h.set_diff_to_next( o.diff_to_next() ) ;
			if( o.has_diff_to_next_chromosome() ) h.set_diff_to_next_chromosome( o.diff_to_next_chromosome() ) ;
			if( o.has_diff_to_next_chromosome_class() ) h.set_diff_to_next_chromosome_class( o.diff_to_next_chromosome_class() ) ;
		}

		{
			AlnStats &a = *rs.mutable_aln_stats() ;
			a.set_reason( o.reason() ) ;

			if( o.has_num_raw_seeds() ) a.set_num_raw_seeds( o.num_raw_seeds() ) ;
			if( o.has_num_grown_seeds() ) a.set_num_grown_seeds( o.num_grown_seeds() ) ;
			if( o.has_num_clumps() ) a.set_num_clumps( o.num_clumps() ) ;
			if( o.has_num_useless() ) a.set_num_useless( o.num_useless() ) ;

			if( o.has_open_nodes_after_alignment() ) a.set_open_nodes_after_alignment( o.open_nodes_after_alignment() ) ;
			if( o.has_closed_nodes_after_alignment() ) a.set_closed_nodes_after_alignment( o.closed_nodes_after_alignment() ) ;
			if( o.has_tracked_closed_nodes_after_alignment() ) a.set_tracked_closed_nodes_after_alignment( o.tracked_closed_nodes_after_alignment() ) ;
		}

		if( o.has_diff_to_next_species() ) rs.set_diff_to_next_species( o.diff_to_next_species() ) ;
		if( o.has_diff_to_next_order() ) rs.set_diff_to_next_order( o.diff_to_next_order() ) ;
		return rs ;
	}
} ;

void Stream::read_next_message( google::protobuf::io::CodedInputStream& cis, const std::string& name )
{
	state_ = invalid ;
	uint32_t tag = 0 ;
	if( cis.ExpectAtEnd() ) {
		throw ParseError( name + " ended unexpectedly" ) ;
	}
	else if( (tag = cis.ReadTag()) )
	{
		uint32_t size = 0 ;
		if( cis.ReadVarint32( &size ) ) 
		{
			int lim = cis.PushLimit( size ) ;
			OldResult ores ;

			if( tag == mk_msg_tag( 4 ) && res_.ParseFromCodedStream( &cis ) )
			{
				cis.PopLimit( lim ) ;
				state_ = have_output ;
				sanitize( *res_.mutable_read() ) ;
				return ;
			}
			if( tag == mk_msg_tag( 2 ) && ores.ParseFromCodedStream( &cis ) )
			{
				cis.PopLimit( lim ) ;
				state_ = have_output ;
				res_ = upgrade( ores ) ;
				sanitize( *res_.mutable_read() ) ;
				return ;
			}
			if( tag == mk_msg_tag( 3 ) && foot_.ParseFromCodedStream( &cis ) )
			{
				cis.PopLimit( lim ) ;
				state_ = end_of_stream ;
				return ; 
			}

			throw ParseError( "deserialization error in " + name ) ;
		}
	}
}

AnfoWriter::AnfoWriter( google::protobuf::io::ZeroCopyOutputStream *zos, const char* fname ) : o_( zos ), name_( fname ), wrote_(0)
{
	o_.WriteRaw( "ANFO", 4 ) ;
}

AnfoWriter::AnfoWriter( int fd, const char* fname, bool expensive )
	: zos_( compress_any( expensive, new FileOutputStream( fd ) ) )
	, o_( zos_.get() ), name_( fname ), wrote_(0)
{
	o_.WriteRaw( "ANFO", 4 ) ;
}

AnfoWriter::AnfoWriter( const char* fname, bool expensive )
	: zos_( compress_any( expensive, new FileOutputStream(
					throw_errno_if_minus1( creat( fname, 0666 ), "opening", fname ) ) ) )
	, o_( zos_.get() ), name_( fname ), wrote_(0)
{
	o_.WriteRaw( "ANFO", 4 ) ;
}

void AnfoWriter::put_result( const Result& r )
{
	write_delimited_message( o_, 4, r ) ; 
	++wrote_ ;
	if( wrote_ % 1024 == 0 )
	{
		stringstream s ;
		s << name_ << ": " << wrote_ << " msgs" ;
		chan_( Console::info, s.str() ) ;
	}
}

void ChunkedWriter::init() 
{
	buf_.resize( 1024*1024 ) ;
	aos_.reset( new ArrayOutputStream( &buf_[0], buf_.size() ) ) ;
	CodedOutputStream o( zos_.get() ) ;
	o.WriteRaw( "ANF1", 4 ) ;
}

ChunkedWriter::ChunkedWriter( std::auto_ptr< ZeroCopyOutputStream > zos, int l, const char* fname ) :
	zos_( zos ), name_( fname ), wrote_(0), method_( method_of(l) ), level_( level_of(l) ) { init() ; }
ChunkedWriter::ChunkedWriter( int fd, int l, const char* fname ) :
	zos_( new FileOutputStream( fd ) ), name_( fname ), wrote_(0), method_( method_of(l) ), level_( level_of(l) ) { init() ; }
ChunkedWriter::ChunkedWriter( const char* fname, int l ) :
	zos_( new FileOutputStream( throw_errno_if_minus1( creat( fname, 0666 ), "opening", fname ) ) ),
	name_( fname ), wrote_(0), method_( method_of(l) ), level_( level_of(l) ) { init() ; }

void ChunkedWriter::flush_buffer( unsigned needed ) 
{
	uint32_t uncomp_size = aos_->ByteCount() ;
	if( uncomp_size == 0 ) return ;
	if( needed && buf_.size() - uncomp_size >= needed+8 ) return ;
	aos_.reset( 0 ) ;

	CodedOutputStream o( zos_.get() ) ;
	o.WriteLittleEndian32( uncomp_size ) ;	// uncompressed size

	if( uncomp_size < 16 || method_ == none ) {
		// very small chunk: confuses compressors, so write uncompressed
		o.WriteLittleEndian32( uncomp_size ) ;	// compressed size & compression method
		o.WriteRaw( &buf_[0], uncomp_size ) ;
	}
	else
	{
		vector< char > tmp ;
		uint32_t comp_size ;
		uLong comp_size_l ;
		switch( method_ )
		{
			case fastlz:
				tmp.resize( uncomp_size * 21 / 20 + 67 ) ;
				comp_size = fastlz_compress_level( level_, &buf_[0], uncomp_size, &tmp[0] ) ;
				break ;

			case gzip:
				comp_size_l = compressBound( uncomp_size ) ;
				tmp.resize( comp_size_l ) ;
				if( Z_OK != compress2( (Bytef*)&tmp[0], &comp_size_l, (const Bytef*)&buf_[0], uncomp_size, level_ ) )
					throw "cannot happen!  overflow in compress2" ;
				comp_size = comp_size_l ;
				break ;

			case bzip:
				comp_size = uncomp_size * 101 / 100 + 601 ;
				tmp.resize( comp_size ) ;
				if( BZ_OK != BZ2_bzBuffToBuffCompress( &tmp[0], &comp_size, &buf_[0], uncomp_size, level_, 0, 0 ) )
					throw "cannot happen!  overflow in BZ2_bzBuffToBuffCompress" ;
				break ;

			default:
				throw "cannot happen!  unknown compression method" ;
		}

		o.WriteLittleEndian32( (uint32_t)(method_) << 28 | comp_size ) ;	// compressed size & compression method
		o.WriteRaw( &tmp[0], comp_size ) ;
	}
	aos_.reset( new ArrayOutputStream( &buf_[0], buf_.size() ) ) ;
}

void ChunkedWriter::put_header( const Header& h )
{
	write_delimited_message( aos_.get(), 1, h ) ;
	Stream::put_header( h ) ;
	flush_buffer() ;
}

void ChunkedWriter::put_footer( const Footer& f ) 
{
	flush_buffer() ;
	int64_t footer_start = zos_->ByteCount() ;
	write_delimited_message( aos_.get(), 3, f ) ;
	flush_buffer() ;
	CodedOutputStream cos( zos_.get() ) ;
	std::cerr << "footer chunk starts at " << footer_start << std::endl ;
	cos.WriteLittleEndian64( footer_start ) ;
	Stream::put_footer( f ) ;
}

void ChunkedWriter::put_result( const Result& r )
{
	flush_buffer( r.ByteSize() ) ;
	write_delimited_message( aos_.get(), 4, r ) ; 
	++wrote_ ;
	if( wrote_ % 1024 == 0 )
	{
		stringstream s ;
		s << name_ << ": " << wrote_ << " msgs" ;
		chan_( Console::info, s.str() ) ;
	}
}

ChunkedWriter::~ChunkedWriter()
{
	flush_buffer() ;
}

bool ChunkedReader::get_next_chunk() 
{
	if( ais_.get() && ais_->ByteCount() < buf_.size() ) return true ;
	ais_.reset( 0 ) ;

	CodedInputStream cis( is_.get() ) ;
	if( cis.ExpectAtEnd() ) return false ;

	uint32_t uncomp_size, comp_size ;
	if( !cis.ReadLittleEndian32( &uncomp_size ) || !cis.ReadLittleEndian32( &comp_size ) ) 
		throw ParseError( "couldn't read chunk header from " + name_ ) ;

	int m = comp_size >> 28 ;
	comp_size &= ~(~0 << 28) ;

	std::cerr << "found chunk of " << comp_size << " bytes, originally " << uncomp_size << ", method " 
		<< (m == ChunkedWriter::none ? "none" : m == ChunkedWriter::fastlz ? "fastlz" : m == ChunkedWriter::gzip ? "gzip" : m ==ChunkedWriter::bzip ? "bzip" : "unknown" ) << std::endl ;

	vector< char > tmp ;
	buf_.resize( uncomp_size ) ;
	tmp.resize( comp_size ) ;
	if( !cis.ReadRaw( &tmp[0], comp_size ) ) 
		throw ParseError( "couldn't read  whole chunk from " + name_ ) ;

	uLongf dlen = uncomp_size ;
	switch( m )
	{
		case ChunkedWriter::none:
			if( comp_size != uncomp_size ) throw "size of uncompressed chunk is wrong" ;
			buf_.swap( tmp ) ;
			break ;

		case ChunkedWriter::fastlz:
			if( (int)uncomp_size != fastlz_decompress( &tmp[0], comp_size, &buf_[0], uncomp_size ) )
				throw "FastLZ decompression failed" ;
			break ;

		case ChunkedWriter::gzip: 
			if( Z_OK != uncompress( (Bytef*)&buf_[0], &dlen, (const Bytef*)&tmp[0], comp_size ) 
					|| uncomp_size != dlen )
				throw "GZip decompression failed" ;
			break ;

		case ChunkedWriter::bzip:
			if( BZ_OK != BZ2_bzBuffToBuffDecompress( &buf_[0], &uncomp_size, &tmp[0], comp_size, 0, 0 ) )
				throw "BZip2 decompression failed" ;
			break ;

		default:
			throw "unknown compression method" ;
	}

	ais_.reset( new ArrayInputStream( &buf_[0], uncomp_size ) ) ;
	return true ;
}

ChunkedReader::ChunkedReader( std::auto_ptr< google::protobuf::io::ZeroCopyInputStream > is, const std::string& name ) : is_( is ), name_( name )
{
	{
		std::string magic ;
		CodedInputStream cis( is_.get() ) ;
		if( !cis.ReadString( &magic, 4 ) || magic != "ANF1" ) 
			throw ParseError( name_ + " is not an ANFO file" ) ;
	}
	if( !get_next_chunk() ) throw ParseError( "EOF before header in " + name_ ) ;

	{
		CodedInputStream cis( ais_.get() ) ;
		uint32_t tag ;
		if( !cis.ReadVarint32( &tag ) || tag != mk_msg_tag( 1 ) || !cis.ReadVarint32( &tag ) )
			throw ParseError( "couldn't read message tag in " +name_ ) ;

		int lim = cis.PushLimit( tag ) ;
		if( !hdr_.ParseFromCodedStream( &cis ) ) 
			throw ParseError( "deserialization error in header of " + name_ ) ;
		cis.PopLimit( lim ) ;
	}

	sanitize( hdr_ ) ;
	++anfo_reader__num_files_ ;

	if( !get_next_chunk() ) throw ParseError( "EOF before payload in " + name_ ) ;
	
	CodedInputStream cis2( ais_.get() ) ;
	read_next_message( cis2, name_ ) ;
}

Result ChunkedReader::fetch_result()
{
	Result r ;
	swap( r, res_ ) ;
	if( get_next_chunk() )
	{
		CodedInputStream cis( ais_.get() ) ;
		read_next_message( cis, name_ ) ;
	}
	else state_ = end_of_stream ;
	return r ;
}

Footer ChunkedReader::fetch_footer()
{
	// skip backwards link
	is_->Skip( 8 ) ;
	return Stream::fetch_footer() ;
}

template <typename E> void nub( google::protobuf::RepeatedPtrField<E>& r )
{
	set<E> s ;
	size_t b = 0, e = r.size(), o = 0 ;
	for( ; b != e ; ++b ) 
	{
		if( s.find( r.Get(b) ) == s.end() )
		{
			s.insert( r.Get(b) ) ;
			if( b != 0 ) swap( *r.Mutable(o), *r.Mutable(b) ) ;
			++o ;
		}
	}
	for( ; o != b ; ++o ) r.RemoveLast() ;
}
template <typename E> void nub( google::protobuf::RepeatedField<E>& r )
{
	set<E> s ;
	size_t b = 0, e = r.size(), o = 0 ;
	for( ; b != e ; ++b ) 
	{
		if( s.find( r.Get(b) ) == s.end() )
		{
			s.insert( r.Get(b) ) ;
			if( b != 0 ) swap( *r.Mutable(o), *r.Mutable(b) ) ;
			++o ;
		}
	}
	for( ; o != b ; ++o ) r.RemoveLast() ;
}

void merge_sensibly( Header& lhs, const Header& rhs )
{
	if( !lhs.IsInitialized() ) lhs = rhs ;
	else {
		bool no_task_id = !lhs.has_sge_task_id() || (rhs.has_sge_task_id() && lhs.sge_task_id() != rhs.sge_task_id()) ;
		bool no_job_id = !lhs.has_sge_job_id() || (rhs.has_sge_job_id() && lhs.sge_job_id() != rhs.sge_job_id()) ;

		bool keep_sort = lhs.is_sorted_by_name() == rhs.is_sorted_by_name() && lhs.has_is_sorted_by_coordinate() == rhs.has_is_sorted_by_coordinate() && (!lhs.has_is_sorted_by_coordinate() || lhs.is_sorted_by_coordinate() == rhs.is_sorted_by_coordinate() ) ;

		lhs.MergeFrom( rhs ) ;
		if( no_task_id ) lhs.clear_sge_task_id() ;
		if( no_job_id ) lhs.clear_sge_job_id() ;
		if( !keep_sort ) { lhs.clear_is_sorted_by_name() ; lhs.clear_is_sorted_by_coordinate() ; }
	}
	sanitize( lhs ) ;
}

void sanitize( Header& hdr )
{
	nub( *hdr.mutable_sge_slicing_index() ) ;
	nub( *hdr.mutable_command_line() ) ;
	nub( *hdr.mutable_config()->mutable_genome_path() ) ;
	nub( *hdr.mutable_config()->mutable_policy() ) ;
	if( hdr.has_sge_slicing_stride() ) 
	{
		set< int > indices ;
		for( int i = 0 ; i != hdr.sge_slicing_index_size() ; ++i )
			indices.insert( hdr.sge_slicing_index(i) ) ;
		
		if( indices.size() == hdr.sge_slicing_stride() )
		{
			hdr.clear_sge_slicing_index() ;
			hdr.clear_sge_slicing_stride() ;
		}
	}
	hdr.clear_was_sorted_by_coordinate() ;
}


//! \brief sanity check, in case we are dealt broken files.
//! This only sanitizes dangerous breakage.  Right now we fix out of
//! range trim points and quality strings of the wrong length.  The
//! fix is always to discard an optional field.
//!
//! \todo Check likelihood arrays and discard them in case of length
//!       mismatch.
void sanitize( Read& rd ) 
{
	unsigned l = rd.sequence().length() ;
	if( rd.has_quality() && rd.quality().length() != l ) rd.clear_quality() ;
	if( rd.has_trim_right() && rd.trim_right() > l ) rd.clear_trim_right() ;
	if( rd.trim_left() > l ) rd.clear_trim_left() ;
}

//! \brief merges aln_stats by adding them up
void merge_sensibly( AlnStats& lhs, const AlnStats& rhs )
{
	// reason: do the sensible thing
	if( lhs.reason() != aligned && rhs.reason() != no_policy && rhs.reason() != no_seeds )
	{
		if( rhs.reason() == aligned ) lhs.clear_reason() ;
		else if( lhs.reason() == no_seeds || lhs.reason() == no_policy ) lhs.set_reason( rhs.reason() ) ;
		else if( lhs.reason() == too_many_seeds && rhs.reason() == bad_alignment ) lhs.set_reason( rhs.reason() ) ;
	}

	// everything else: add it up
	lhs.set_num_raw_seeds( lhs.num_raw_seeds() + rhs.num_raw_seeds() ) ;
	lhs.set_num_grown_seeds( lhs.num_grown_seeds() + rhs.num_grown_seeds() ) ;
	lhs.set_num_clumps( lhs.num_clumps() + rhs.num_clumps() ) ;
	lhs.set_num_useless( lhs.num_useless() + rhs.num_useless() ) ;
	lhs.set_open_nodes_after_alignment( lhs.open_nodes_after_alignment() + rhs.open_nodes_after_alignment() ) ;
	lhs.set_closed_nodes_after_alignment( lhs.closed_nodes_after_alignment() + rhs.closed_nodes_after_alignment() ) ;
	lhs.set_tracked_closed_nodes_after_alignment( lhs.tracked_closed_nodes_after_alignment() + rhs.tracked_closed_nodes_after_alignment() ) ;
}

//! \brief merges two hits by keeping the better one
void merge_sensibly( Hit& lhs, const Hit& rhs )
{
	// take better hit, recalculate diff_to_next{,_chromosome{,_class}}
	if( lhs.score() <= rhs.score() )
	{
		// left is better
		if( !lhs.has_diff_to_next() || lhs.score() + lhs.diff_to_next() > rhs.score() )
			lhs.set_diff_to_next( rhs.score() - lhs.score() ) ;

		//! \todo diff to chromosome, chromosome class? dunno...
	}
	else
	{
		// right is better
		if( !rhs.has_diff_to_next() || rhs.score() + rhs.diff_to_next() > lhs.score() )
		{
			int d = lhs.score() - rhs.score() ;
			lhs = rhs ;
			lhs.set_diff_to_next( d ) ;
		}
		else lhs = rhs ;

		//! \todo diff to chromosome, chromosome class? dunno...
	}
}

//! \brief merges two results, keeping the best hit
void merge_sensibly( Result& lhs, const Result& rhs )
{
	// How to merge what...
	// - read: all equal, no merging needed

	// - member: concatenate and remove doubles
	lhs.mutable_member()->MergeFrom( rhs.member() ) ;
	nub( *lhs.mutable_member() ) ;

	// - hits: merge those for the same genome, concatenate the rest
	for( int j = 0 ; j != rhs.hit_size() ; ) {
		for( int i = 0 ; i != lhs.hit_size() ; ++i ) {
			if( rhs.hit(j).genome_name() == lhs.hit(i).genome_name() ) {
				merge_sensibly( *lhs.mutable_hit(i), rhs.hit(j) ) ;
				goto next ;
			}
		}
		*lhs.add_hit() = rhs.hit(j) ;
next:
		++j ;
	}

	// - aln_stats: just add them up, this is for debugging only anyway
	merge_sensibly( *lhs.mutable_aln_stats(), rhs.aln_stats() ) ; 

	// - diff_to_next_species, diff_to_next_order: TODO
}

//! \brief merges two footers
//! Anything unexpected is simply merged, the resulting exit code is the
//! logical or of the two inputs.
void merge_sensibly( output::Footer& lhs, const output::Footer& rhs )
{
	int exit_code = lhs.exit_code() | rhs.exit_code() ;
	lhs.MergeFrom( rhs ) ;
	lhs.set_exit_code( exit_code ) ;
}

bool has_hit_to( const output::Result& r, const char* g )
{
	if( !g ) return r.hit_size() > 0 ;

	for( int i = 0 ; i != r.hit_size() ; ++i )
		if( r.hit(i).genome_name() == g )
			return true ;

	return false ;
}

const output::Hit& hit_to( const output::Result& r, const char* g )
{
	if( !g ) {
		if( r.hit_size() ) {
			const output::Hit *h = &r.hit(0) ;
			for( int i = 1 ; i != r.hit_size() ; ++i )
				if( r.hit(i).score() < h->score() )
					h = &r.hit(i) ; 
			return *h ;
		}
	}
	else 
		for( int i = 0 ; i != r.hit_size() ; ++i )
			if( r.hit(i).genome_name() == g )
				return r.hit(i) ;
	throw "hit_to: no suitable hit" ;
}

output::Hit* mutable_hit_to( output::Result* r, const char* g )
{
	if( !g ) {
		if( r->hit_size() ) {
			output::Hit *h = r->mutable_hit(0) ;
			for( int i = 1 ; i != r->hit_size() ; ++i )
				if( r->hit(i).score() < h->score() )
					h = r->mutable_hit(i) ; 
			return h ;
		}
	}
	else
		for( int i = 0 ; i != r->hit_size() ; ++i )
			if( r->hit(i).genome_name() == g )
				return r->mutable_hit(i) ;

	Hit *h = r->add_hit() ;
	if( g ) h->set_genome_name( g ) ;
	return h ;
}


bool ScoreFilter::xform( Result& r )
{
	int ix_in = 0, ix_out = 0 ;
	while( ix_in != r.hit_size() )
	{
		// keep hits if we're actually looking for a specific genome and
		// they hit the wrong one or if their score is good (small) enough
		if( ( genome_ && *genome_ && r.hit(ix_in).genome_name() != genome_ ) ||
				( slope_ * ( len_from_bin_cigar( r.hit(ix_in).cigar() )
							 - intercept_ ) >= r.hit(ix_in).score() ) )
		{
			if( ix_in != ix_out ) *r.mutable_hit(ix_out) = r.hit(ix_in) ;
			++ix_out ;
		}
		++ix_in ;
	}

	if( !ix_out )
	{
		r.mutable_aln_stats()->set_reason( bad_alignment ) ;
		r.clear_hit() ;
	}
	else while( ix_out != r.hit_size() ) r.mutable_hit()->RemoveLast() ;
	return true ;
}

bool MapqFilter::xform( Result& r ) {
	return has_hit_to( r, g_ ) &&
		( !hit_to( r, g_).has_diff_to_next() ||
		  hit_to( r, g_ ).diff_to_next() >= minmapq_ ) ;
}

bool LengthFilter::xform( Result& r ) {
	int len = ( r.read().has_trim_right() ? r.read().trim_right() : r.read().sequence().size() ) - r.read().trim_left() ;
	if( r.hit_size() && len < minlength_ )
	{
		r.clear_hit() ;
		if( r.has_aln_stats() ) r.mutable_aln_stats()->set_reason( no_policy ) ;
	}
	return true ;
}

bool HitFilter::xform( Result& r ) 
{
	if( has_hit_to( r, g_ ) )
		return !s_ || !*s_ || hit_to( r, g_ ).sequence() == s_ ;

	return false ;
}

bool Subsample::xform( Result& ) 
{
	return f_ >= drand48() ;
}

bool RmdupStream::is_duplicate( const Result& lhs, const Result& rhs ) 
{
	if( !has_hit_to( lhs, g_ ) || !has_hit_to( rhs, g_ )
			|| lhs.read().sequence().size() != rhs.read().sequence().size() )
		return false ;

	const output::Hit &l = hit_to( lhs, g_ ), &r = hit_to( rhs, g_ ) ;

	return l.genome_name() == r.genome_name() && l.sequence() == r.sequence()
		&& l.start_pos() == r.start_pos() && l.aln_length() == r.aln_length() ;
}

//! \todo How do we deal with ambiguity codes?  What's the meaning of
//!       their quality scores anyway?
void RmdupStream::add_read( const Result& rhs ) 
{
	// if the new member is itself a cluster, we add its member reads,
	// not the single synthetic one
	if( rhs.member_size() ) 
		for( int i = 0 ; i != rhs.member_size() ; ++i )
			*cur_.add_member() = rhs.member(i) ;
	else
		*cur_.add_member() = rhs.read() ;

	for( size_t i = 0 ; i != rhs.read().sequence().size() ; ++i )
	{
		int base = -1 ;
		switch( rhs.read().sequence()[i] ) {
			case 'a': case 'A': base = 0 ; break ;
			case 'c': case 'C': base = 1 ; break ;
			case 't': case 'T':
			case 'u': case 'U': base = 2 ; break ;
			case 'g': case 'G': base = 3 ; break ;
		}

		Logdom qual = Logdom::from_phred( rhs.read().has_quality() ? rhs.read().quality()[i] : 30 ) ;
		for( int j = 0 ; j != 4 ; ++j )
			// XXX distribute errors sensibly
			quals_[j].at(i) *= j != base ? qual / 3 : 1 - qual ;
	}
}

//! \brief returns the next result
//! This also has to change state.  If we haven't seen a footer, we
//! request more input.  Else we're at end of stream.
Result RmdupStream::fetch_result() 
{
	Result r ;
	swap( r, res_ ) ;
	state_ = foot_.IsInitialized() ? end_of_stream : state_ = need_input ;
	return r ;
}

namespace {
    void limit_quality( Read& r, uint8_t maxq ) 
    {
		if( !r.has_quality() ) return ;
        for( size_t i = 0 ; i != r.quality().size() ; ++i )
            if( r.quality()[i] > maxq )
					(*r.mutable_quality())[i] = maxq ;
    }
} ;

//! \brief signals end of input stream
//! After the footer no more input is possible.  We were waiting for
//! input, so no output was available.  cur_ might contain valid data,
//! and if so, we call a consensus and offer it as output.  Else we
//! signal end of stream.
void RmdupStream::put_footer( const Footer& f ) { 
	Stream::put_footer( f ) ;
	if( cur_.IsInitialized() ) 
	{
		state_ = have_output ;
		call_consensus() ;
		swap( cur_, res_ ) ;
	}
}

//! \brief receives a result record and merges it if appropriate
//!
//! There are the following possibilities what to do here:
//! - A result with a bad alignment is passed through (means it is
//!   stored in res_ and output becomes available).
//! - If cur_ is invalid, the result is stored there and cur_ becomes
//!   valid.
//! - A result with correct coordinates (according to is_duplicate) is
//!   directly merged into cur_, no output becomes available.
//! - Anything else cannot be merged, so a consensus is called, cur_
//!   moves to res_, next moves to cur_, and output becomes available.
//!
//! \todo In principle, search for duplicates can be repeated (e.g.
//!       after sorting on a different genome coordinate), and this is
//!       supported; it is not really tested, though.

void RmdupStream::put_result( const Result& next ) 
{
	if( !has_hit_to( next, g_ ) || hit_to( next, g_ ).score() >
			slope_ * ( len_from_bin_cigar( hit_to( next, g_ ).cigar() ) - intercept_ ) )
	{
		// bad alignment -- this one passes through without merging
		// we clamp qualities, though
		res_ = next ;
		Read &r = *res_.mutable_read() ;
        limit_quality( r, maxq_ ) ;
		state_ = have_output ;
	}
	else if( !cur_.IsInitialized() ) {
		cur_ = next ;
		for( size_t i = 0 ; i != 4 ; ++i )
		{
			quals_[i].clear() ;
			quals_[i].resize( cur_.read().sequence().size() ) ;
		}
	}
	else if( is_duplicate( cur_, next ) )
	{
		// Merge them.  If cur is a plain result, turn it into a
		// degenerate merged one first...
		if( cur_.member_size() == 0 )
		{
			add_read( cur_ ) ;
			cur_.mutable_read()->set_seqid( "C_" + cur_.read().seqid() ) ;
			cur_.mutable_read()->clear_description() ;
		}
		// Merge the new one.  No state change necessary, we continue to
		// request input.
		add_read( next ) ;
	}
	else
	{
		// Nothing to match.  Call a consensus for cur_, then move it to
		// res_.  State that output is available, store new result in
		// cur_.
		call_consensus() ;
		swap( res_, cur_ ) ;
		cur_ = next ;
		for( size_t i = 0 ; i != 4 ; ++i )
		{
			quals_[i].clear() ;
			quals_[i].resize( cur_.read().sequence().size() ) ;
		}
		state_ = have_output ;
	}
}

void RmdupStream::call_consensus()
{
	if( !cur_.member_size() ) {
        limit_quality( *cur_.mutable_read(), maxq_ ) ;
        return ;
    }

	cur_.mutable_read()->clear_sequence() ;
	cur_.mutable_read()->clear_quality() ;
	for( size_t i = 0 ; i != quals_[0].size() ; ++i )
	{
		// select base with highest quality
		size_t m = 0 ;
		for( size_t j = 1 ; j != 4 ; ++j )
			if( quals_[j].at( i ) > quals_[m].at( i ) ) m = j ;

		// but calculate the error probability from _all_others_ to
		// retain precision

		Logdom denom = Logdom::null(), num = Logdom::null() ;
		for( size_t j = 0 ; j != 4 ; ++j )
		{
			denom += quals_[j].at( i ) ;
			if( j != m ) num += quals_[j].at( i ) ;
		}

		int qscore = (num/denom).to_phred() ;
		if( qscore > maxq_ ) qscore = maxq_ ;

		cur_.mutable_read()->mutable_sequence()->push_back( m["ACTG"] ) ;
		cur_.mutable_read()->mutable_quality()->push_back( qscore ) ;
	}
}

void ConcatStream::add_stream( StreamHolder s )
{
	merge_sensibly( hdr_, s->fetch_header() ) ;
	if( s->get_state() == have_output )
	{
		streams_.push_back( s ) ;
		state_ = have_output ;
	}
	else
	{
		if( state_ == invalid ) state_ = end_of_stream ;
		merge_sensibly( foot_, s->fetch_footer() ) ;
	}
}

Result ConcatStream::fetch_result()
{
	Result r = streams_[0]->fetch_result() ;
	if( streams_[0]->get_state() != have_output )
	{
		merge_sensibly( foot_, streams_[0]->fetch_footer() ) ;
		streams_.pop_front() ;
	}
	if( streams_.empty() ) state_ = end_of_stream ;
	return r ;
}

bool QualFilter::xform( Result& r )
{
	if( r.read().has_quality() ) 
		for( size_t i = 0 ; i != r.read().sequence().size() && i != r.read().quality().size() ; ++i )
			if( r.read().quality()[i] < q_ ) (*r.mutable_read()->mutable_sequence())[i] = 'N' ;
	return true ;
}

namespace {
	class StreamWithProgress : public google::protobuf::io::ZeroCopyInputStream 
	{
		private:
			std::auto_ptr< google::protobuf::io::ZeroCopyInputStream > is_ ;
			std::string name_ ;
			int64_t total_, read_ ;
			Chan out_ ;

			bool check( bool p ) {
				if( !p ) out_.close() ;
				else if( is_->ByteCount() >> 16 != read_ >> 16 )
				{
					read_ = is_->ByteCount() ;
					stringstream s ;
					s << name_ << ": " << read_ ;
					if( total_ ) s << '/' << total_ ;
					s << " Bytes" ;
					if( total_ ) s << " (" << read_*100/total_ << "%)" ;
					out_( Console::info, s.str() ) ;
				}
				return p ;
			}

		public:
			StreamWithProgress( std::auto_ptr< google::protobuf::io::ZeroCopyInputStream > is, const char* name, int64_t total )
				: is_( is ), name_( basename( name ) ), total_( total ), read_( 0 ) {}

			virtual bool Next( const void **data, int *size ) { return check( is_->Next( data, size ) ) ; }
			virtual void BackUp( int count ) { is_->BackUp( count ) ; }
			virtual bool Skip( int count ) { return check( is_->Skip( count ) ) ; }
			virtual int64_t ByteCount() const { return is_->ByteCount() ; }


	} ;

	StreamHolder make_input_stream_( std::auto_ptr< google::protobuf::io::ZeroCopyInputStream > is, const char* name, bool solexa_scores, char origin )
	{
		// peek into stream, but put it back.  then check magic numbers and
		// create the right stream
		const void* p ; int l ;
		if( !is->Next( &p, &l ) ) return new Stream ;
		is->BackUp( l ) ;

		const uint8_t* q = (const uint8_t*)p ;
		if( l >= 4 && q[0] == 'A' && q[1] == 'N' && q[2] == 'F' && q[3] == 'O' )
			return new AnfoReader( is, name ) ;

		else if( l >= 4 && q[0] == 'A' && q[1] == 'N' && q[2] == 'F' && q[3] == '1' )
			return new ChunkedReader( is, name ) ;

		else if( l >= 3 && q[0] == 'B' && q[1] == 'Z' && q[2] == 'h' )
		{
			std::auto_ptr< google::protobuf::io::ZeroCopyInputStream > bs( new BunzipStream( is ) ) ;
			return make_input_stream_( bs, name, solexa_scores, origin ) ;
		}
		else if( l >= 2 && q[0] == 31 && q[1] == 139 )
		{
			std::auto_ptr< google::protobuf::io::ZeroCopyInputStream > zs( new InflateStream( is ) ) ;
			return make_input_stream_( zs, name, solexa_scores, origin ) ;
		}
		else return new FastqReader( is, solexa_scores, origin ) ;
	}
} ;


StreamHolder make_input_stream( const char *name, bool solexa_scores, char origin )
{
	return name && *name && strcmp( name, "-" ) 
		? make_input_stream( throw_errno_if_minus1(
				open( name, O_RDONLY ), "opening ", name ), name, solexa_scores, origin ) 
		: make_input_stream( dup( 0 ), "<stdin>", solexa_scores, origin ) ;
}

StreamHolder make_input_stream( int fd, const char *name, bool solexa_scores, char origin )
{
	struct stat st ;
	std::auto_ptr< google::protobuf::io::FileInputStream > s( new google::protobuf::io::FileInputStream( fd ) ) ;
	s->SetCloseOnDelete( true ) ;
	return make_input_stream( std::auto_ptr< google::protobuf::io::ZeroCopyInputStream >(s), name, fstat( fd, &st ) ? -1 : st.st_size, solexa_scores, origin ) ;
}

StreamHolder make_input_stream( std::auto_ptr< google::protobuf::io::ZeroCopyInputStream > is, const char *name,
		int64_t total, bool solexa_scores, char origin )
{
	std::auto_ptr< google::protobuf::io::ZeroCopyInputStream > s( new StreamWithProgress( is, name, total ) ) ;
	return make_input_stream_( s, name, solexa_scores, origin ) ;
}

} ; // namespace
