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
#include "misc_streams.h"
#include "stream.h"
#include "util.h"

extern "C" {
#include "fastlz.h"
}

#include <google/protobuf/repeated_field.h>

#include <algorithm>
#include <cctype>
#include <numeric>
#include <set>

#if HAVE_FCNTL_H
#include <fcntl.h>
#endif

#include <iostream>

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

int transfer( Stream& in, Stream& out ) 
{
	out.put_header( in.fetch_header() ) ;
	while( in.get_state() == Stream::have_output && out.get_state() == Stream::need_input )
		out.put_result( in.fetch_result() ) ;
	out.put_footer( in.fetch_footer() ) ;
	return out.fetch_footer().exit_code() ;
}

int anfo_reader__num_files_ = 0 ;

namespace {
	string basename( const string& s )
	{
		string::size_type p = s.rfind( '/' ) ;
		return p != string::npos ? s.substr(p+1) : s ;
	}
} ;

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
                if( state_ == end_of_stream ) is_.reset( 0 ) ;
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

		rs.mutable_members()->MergeFrom( o.member() ) ;

		if( o.has_best_hit() ) *rs.add_hit() = upgrade( o.best_hit() ) ;
		if( o.has_best_to_genome() ) {
			Hit &h = *rs.add_hit() ;
			h = upgrade( o.best_to_genome() ) ;
			if( o.has_diff_to_next() ) h.set_diff_to_next( o.diff_to_next() ) ;
			if( o.has_diff_to_next_chromosome() ) h.set_diff_to_next_chromosome( o.diff_to_next_chromosome() ) ;
			if( o.has_diff_to_next_chromosome_class() ) h.set_diff_to_next_chromosome_class( o.diff_to_next_chromosome_class() ) ;
		}

		{
			AlnStats &a = *rs.add_aln_stats() ;
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
				sanitize( res_ ) ;
				return ;
			}
			if( tag == mk_msg_tag( 2 ) && ores.ParseFromCodedStream( &cis ) )
			{
				cis.PopLimit( lim ) ;
				state_ = have_output ;
				res_ = upgrade( ores ) ;
				sanitize( res_ ) ;
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

/*
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
*/

void ChunkedWriter::init() 
{
    // sensible buffer size: big enough to make compression worthwhile,
    // small enough that BZip won't split it again
	buf_.resize( 850000 ) ;
	aos_.reset( new ArrayOutputStream( &buf_[0], buf_.size() ) ) ;
	CodedOutputStream o( zos_.get() ) ;
	o.WriteRaw( "ANF1", 4 ) ;
}

ChunkedWriter::ChunkedWriter( const pair< ZeroCopyOutputStream*, string >& p, int l ) :
	zos_( p.first ), name_( p.second ), wrote_(0), method_( method_of(l) ), level_( level_of(l) ) { init() ; }
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

#if HAVE_LIBZ && HAVE_ZLIB_H
			case gzip:
				comp_size_l = compressBound( uncomp_size ) ;
				tmp.resize( comp_size_l ) ;
				if( Z_OK != compress2( (Bytef*)&tmp[0], &comp_size_l, (const Bytef*)&buf_[0], uncomp_size, level_ ) )
					throw "cannot happen!  overflow in compress2" ;
				comp_size = comp_size_l ;
				break ;
#endif

#if HAVE_LIBBZ2 && HAVE_BZLIB_H
			case bzip:
                // docu says 5% is overkill.  didn't work with 1% reserve, though...
				comp_size = uncomp_size * 21 / 20 + 601 ;
				tmp.resize( comp_size ) ;
				if( BZ_OK != BZ2_bzBuffToBuffCompress( &tmp[0], &comp_size, &buf_[0], uncomp_size, level_, 0, 0 ) )
					throw "cannot happen!  overflow in BZ2_bzBuffToBuffCompress" ;
				break ;
#endif

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
	cos.WriteLittleEndian64( footer_start ) ;
	Stream::put_footer( f ) ;

	stringstream ss ;
	ss << name_ << ": footer chunk starts at " << footer_start ;
	console.output( Console::notice, ss.str() ) ;
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
	if( ais_.get() && (unsigned)ais_->ByteCount() < buf_.size() ) return true ;
	ais_.reset( 0 ) ;

	CodedInputStream cis( is_.get() ) ;
	if( cis.ExpectAtEnd() ) return false ;

	uint32_t uncomp_size, comp_size ;
	if( !cis.ReadLittleEndian32( &uncomp_size ) || !cis.ReadLittleEndian32( &comp_size ) ) 
		throw ParseError( "couldn't read chunk header from " + name_ ) ;

	int m = comp_size >> 28 ;
	comp_size &= ~(~0 << 28) ;

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
#if HAVE_LIBZ && HAVE_ZLIB_H
			if( Z_OK != uncompress( (Bytef*)&buf_[0], &dlen, (const Bytef*)&tmp[0], comp_size ) 
					|| uncomp_size != dlen )
				throw "GZip decompression failed" ;
#else
			throw "GZip'ed chunk found, but no zlib support present." ;
#endif
			break ;

		case ChunkedWriter::bzip:
#if HAVE_LIBBZ2 && HAVE_BZLIB_H
			if( BZ_OK != BZ2_bzBuffToBuffDecompress( &buf_[0], &uncomp_size, &tmp[0], comp_size, 0, 0 ) )
				throw "BZip2 decompression failed" ;
#else
			throw "Bzip'ed chunk found, but no libbz2 support present." ;
#endif
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
        if( state_ == end_of_stream ) 
        {
            // skip backwards link (just because...)
            is_->Skip( 8 ) ;
            // delete and close(!) input(-file)
            is_.reset( 0 ) ;
            --anfo_reader__num_files_ ;
        }
	}
	else state_ = end_of_stream ;
	return r ;
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
			if( b != o ) swap( *r.Mutable(o), *r.Mutable(b) ) ;
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
			if( b != o ) swap( *r.Mutable(o), *r.Mutable(b) ) ;
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

		bool keep_sort = lhs.is_sorted_by_name() == rhs.is_sorted_by_name() &&
			lhs.has_is_sorted_by_all_genomes() == rhs.has_is_sorted_by_all_genomes() &&
			lhs.is_sorted_by_coordinate_size() == rhs.is_sorted_by_coordinate_size() &&
			equal( lhs.is_sorted_by_coordinate().begin(), lhs.is_sorted_by_coordinate().end(), rhs.is_sorted_by_coordinate().begin() ) ;

		lhs.MergeFrom( rhs ) ;
		if( no_task_id ) lhs.clear_sge_task_id() ;
		if( no_job_id ) lhs.clear_sge_job_id() ;
		if( !keep_sort ) {
			lhs.clear_is_sorted_by_name() ;
			lhs.clear_is_sorted_by_coordinate() ;
			lhs.clear_is_sorted_by_all_genomes() ;
		}
	}
	sanitize( lhs ) ;
}

void sanitize( Header& hdr )
{
	nub( *hdr.mutable_command_line() ) ;
	nub( *hdr.mutable_config()->mutable_policy() ) ;
	hdr.clear_was_sorted_by_coordinate() ;
	if( hdr.is_sorted_by_coordinate_size() == 1 && hdr.is_sorted_by_coordinate(0) == "" )
	{
		hdr.clear_is_sorted_by_coordinate() ;
		hdr.set_is_sorted_by_all_genomes( true ) ;
	}
}


//! \brief sanity check, in case we are dealt broken files.
//! This only sanitizes dangerous breakage.  Right now we fix out of
//! range trim points and quality strings of the wrong length.  The
//! fix is always to discard an optional field.
//!
//! \todo Check likelihood arrays and discard them in case of length
//!       mismatch.
void sanitize_read( Read& rd ) 
{
	unsigned l = rd.sequence().length() ;
	if( rd.has_quality() && rd.quality().length() != l ) rd.clear_quality() ;
	if( rd.has_trim_right() && rd.trim_right() > l ) rd.clear_trim_right() ;
	if( rd.trim_left() > l ) rd.clear_trim_left() ;
    if( rd.has_description() && rd.description().empty() ) rd.clear_description() ;
}

void sanitize( Result& r ) 
{
    for( int i = 0 ; i != r.hit_size() ; ++i )
    {
        string& s = *r.mutable_hit(i)->mutable_genome_name() ;
        for( size_t j = 0 ; j != s.size() ; ++j )
            s[j] = tolower( s[j] ) ;
    }
    sanitize_read( *r.mutable_read() ) ;
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
	lhs.mutable_members()->MergeFrom( rhs.members() ) ;
	lhs.set_nmembers( lhs.nmembers() + rhs.nmembers() ) ;
	nub( *lhs.mutable_members() ) ;
	if( lhs.members_size() >= 8 ) {
		lhs.set_nmembers( lhs.nmembers() + lhs.members_size() ) ;
		lhs.clear_members() ;
	}

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

	// - aln_stats: just concat them, this is for debugging only anyway
	for( int i = 0 ; i != rhs.aln_stats_size() ; ++i ) *lhs.add_aln_stats() = rhs.aln_stats(i) ;

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


const output::Hit* hit_to( const output::Result& r )
{
	if( r.hit_size() ) {
		const output::Hit *h = &r.hit(0) ;
		for( int i = 1 ; i != r.hit_size() ; ++i )
			if( r.hit(i).score() < h->score() )
				h = &r.hit(i) ; 
		return h ;
	}
	return 0 ;
}

const output::Hit* hit_to( const output::Result& r, const string& g )
{
	for( int i = 0 ; i != r.hit_size() ; ++i )
		if( r.hit(i).genome_name() == g )
			return &r.hit(i) ;
	return 0 ;
}

output::Hit* mutable_hit_to( output::Result* r )
{
	if( r->hit_size() ) {
		output::Hit *h = r->mutable_hit(0) ;
		for( int i = 1 ; i != r->hit_size() ; ++i )
			if( r->hit(i).score() < h->score() )
				h = r->mutable_hit(i) ; 
		return h ;
	}

	return r->add_hit() ;
}
	
output::Hit* mutable_hit_to( output::Result* r, const string& g )
{
	for( int i = 0 ; i != r->hit_size() ; ++i )
		if( r->hit(i).genome_name() == g )
			return r->mutable_hit(i) ;

	Hit *h = r->add_hit() ;
	h->set_genome_name( g ) ;
	return h ;
}

void HitFilter::put_result( const Result& res )
{
	res_ = res ;

	int ix_in = 0, ix_out = 0 ;
	while( ix_in != res_.hit_size() )
	{
		// keep hits if we're actually looking for a specific genome and
		// they hit the wrong one or if they actually fit our predicate
		if( ( !gs_.empty() && !contains( gs_, res_.hit(ix_in).genome_name() ) ) || keep( res_.hit(ix_in) ) )
		{
			if( ix_in != ix_out ) swap( *res_.mutable_hit(ix_out), *res_.mutable_hit(ix_in) ) ;
			++ix_out ;
		}
		++ix_in ;
	}

	if( !ix_out ) res_.clear_hit() ;
	else while( ix_out != res_.hit_size() ) res_.mutable_hit()->RemoveLast() ;
	
	state_= have_output ;
}

bool OnlyGenome::xform( Result& res )
{
	int ix_in = 0, ix_out = 0 ;
	while( ix_in != res_.hit_size() )
	{
		// keep hits if we're actually looking for a specific genome and
		// they hit the wrong one or if they actually fit our predicate
		if( contains( gs_, res.hit(ix_in).genome_name() ) )
		{
			if( ix_in != ix_out ) swap( *res.mutable_hit(ix_out), *res.mutable_hit(ix_in) ) ;
			++ix_out ;
		}
		++ix_in ;
	}

	if( !ix_out ) res.clear_hit() ;
	else while( ix_out != res_.hit_size() ) res.mutable_hit()->RemoveLast() ;
	return true ;
}


bool ScoreFilter::keep( const Hit& h )
{ return slope_ * ( len_from_bin_cigar( h.cigar() ) - intercept_ ) >= h.score() ; }

static int effective_length( const Read& rd )
{
    return rd.has_trim_right()
        ? rd.trim_right() - rd.trim_left()
        : rd.sequence().length() - rd.trim_left() ;
}

bool TotalScoreFilter::xform( Result& r )
{
    int score = 0 ;
    for( int i = 0 ; i != r.hit_size() ; ++i )
        if( contains( gs_, r.hit(i).genome_name() ) )
            score += r.hit(i).score() ;
    return slope_ * ( effective_length( r.read() ) - intercept_ ) >= score ;
}

bool MapqFilter::keep( const Hit& h )
{ return !h.has_diff_to_next() || h.diff_to_next() >= minmapq_ ; }

bool QualFilter::xform( Result& h )
{
	const Read& r = h.read() ;
	return r.has_quality() && accumulate(
			r.quality().begin(),
			r.quality().end(),
			static_cast<int>( 0 ),
			plus<int>() ) 
		>= r.quality().size() * minqual_ ;
}

bool LengthFilter::xform( Result& r )
{
	int len = ( r.read().has_trim_right() ? r.read().trim_right() : r.read().sequence().size() ) - r.read().trim_left() ;
	if( r.hit_size() && len < minlength_ ) r.clear_hit() ;
	return true ;
}

namespace {
	bool good_hit( const Hit& h, const vector<string>& gs, const vector<string>& ss )
	{
		return ( gs.empty() || contains( gs, h.genome_name() ) ) 
			&& ( ss.empty() || contains( ss, h.sequence() ) ) ;
	}
} ;

bool RequireHit::xform( Result& r ) 
{
	for( int i = 0 ; i != r.hit_size() ; ++i )
		if( good_hit( r.hit(i), gs_, ss_ ) ) return true ;
	return false ;
}

bool RequireBestHit::xform( Result& r ) { const Hit *h = hit_to( r ) ; return h && good_hit( *h, gs_, ss_ ) ; }


bool Subsample::xform( Result& ) 
{
	return f_ >= drand48() ;
}

bool RmdupStream::is_duplicate( const Result& lhs, const Result& rhs ) const
{
	const output::Hit *l = hit_to( lhs, gs_.begin(), gs_.end() ), *r = hit_to( rhs, gs_.begin(), gs_.end() ) ;
	by_genome_coordinate comp ;
	return l && r && 
		!comp.compare( l, r, lhs.read(), rhs.read() ) &&
		!comp.compare( r, l, rhs.read(), lhs.read() ) ;
}

//! \todo How do we deal with ambiguity codes?  What's the meaning of
//!       their quality scores anyway?
void RmdupStream::add_read( const Result& rhs ) 
{
	// if the new member is itself a cluster, we add its member reads,
	// not the single synthetic one
	if( rhs.nmembers() )
	{
		cur_.set_nmembers( rhs.nmembers() + cur_.nmembers() + cur_.members_size() ) ;
		cur_.clear_members() ;
	}
	else if( rhs.members_size() ) 
		for( int i = 0 ; i != rhs.members_size() ; ++i )
			*cur_.add_members() = rhs.members(i) ;
	else
		*cur_.add_members() = rhs.read() ;

	if( cur_.members_size() >= 8 ) {
		cur_.set_nmembers( cur_.nmembers() + cur_.members_size() ) ;
		cur_.clear_members() ;
	}

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

inline int RmdupStream::max_score( const Hit* h ) const
{ return int( 0.5 + slope_ * ( len_from_bin_cigar( h->cigar() ) - intercept_ ) ) ; }

//! \brief receives a result record and merges it if appropriate
//!
//! There are the following possibilities what to do here:
//! - A result with a bad alignment, but correct coordinates is passed
//!   through (means it is stored in res_ and output becomes available).
//! - If cur_ is invalid, the result is stored there and cur_ becomes
//!   valid.
//! - A result with correct coordinates (according to is_duplicate) is
//!   directly merged into cur_, no output becomes available.
//! - Anything else cannot be merged, so a consensus is called, cur_
//!   moves to res_, next moves to cur_, and output becomes available.

void RmdupStream::put_result( const Result& next ) 
{
	const Hit *h = hit_to( next, gs_.begin(), gs_.end() ) ;
	const Hit *h0 = hit_to( cur_, gs_.begin(), gs_.end() ) ;

	// first check: is anything buffered?
	if( !cur_.IsInitialized() ) 
	{
		// empty buffer.  if we now get a good alignment, we store it.
		// anything else passed through
		if( h && h->score() <= max_score( h ) )
		{
			cur_ = next ;
			for( size_t i = 0 ; i != 4 ; ++i )
			{
				quals_[i].clear() ;
				quals_[i].resize( cur_.read().sequence().size() ) ;
			}
		}
		else
		{
			// bad alignment or none at all -- this one passes through
			// without merging.  we clamp qualities, though
			res_ = next ;
			Read &r = *res_.mutable_read() ;
			limit_quality( r, maxq_ ) ;
			state_ = have_output ;
		}
	}
	// something is stored.  If it has a bad alignment, we need to get
	// rid of it, which makes room to store the next alignment.  (Very
	// annoying, but it's a good thing we never need more than one slot
	// for buffering...)
	else if( !h0 || h0->score() > max_score( h0 ) )
	{
		swap( res_, cur_ ) ;
		Read &r = *res_.mutable_read() ;
		limit_quality( r, maxq_ ) ;
		state_ = have_output ;

		cur_ = next ;
		for( size_t i = 0 ; i != 4 ; ++i )
		{
			quals_[i].clear() ;
			quals_[i].resize( cur_.read().sequence().size() ) ;
		}
	}
	// we got something stored, and it is known to be eligible for
	// merging.  if we get a matched alignment, we need to check if the
	// new one is any good...
	else if( h && is_duplicate( cur_, next ) ) 
	{
		if( h->score() <= max_score( h ) ) 
		{
			// a good one, we need to actually merge it.  If cur_ is a
			// plain result, turn it into a degenerate merged one
			// first...
			if( cur_.members_size() + cur_.nmembers() == 0 )
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
			// the new one is crap.  Pass it on, keep whatever is in
			// the accumulator.
			res_ = next ;
			Read &r = *res_.mutable_read() ;
			limit_quality( r, maxq_ ) ;
			state_ = have_output ;
		}
	}
	else
	{
		// we got something that has wrong coordinates or no alignment
		// at all.  Call a consensus for cur_, then move it to res_.
		// State that output is available, and store new result in cur_.
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
	if( !cur_.members_size() && !cur_.nmembers() ) {
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

//! \todo Tends to return the header from the first input stream, but
//! will accumulate further headers (and throw them away).  This is
//! necessary, because we cannot look ahead to further headers.  Clearly
//! a better solution would be great, but doesn't readily present
//! itself.  Could adding arbitrary header information to the footer
//! work?
Header ConcatStream::fetch_header()
{
	for( state_ = invalid ; state_ == invalid ; )
	{
		if( streams_.empty() ) state_ = end_of_stream ;
        else {
            merge_sensibly( hdr_, streams_[0]->fetch_header() ) ;
            if( streams_[0]->get_state() == have_output ) state_ = have_output ;
            else 
            {
                merge_sensibly( foot_, streams_[0]->fetch_footer() ) ;
                streams_.pop_front() ;
            }
        }
	}
	return hdr_ ;
}

Result ConcatStream::fetch_result()
{
    Result r = streams_[0]->fetch_result() ;
	for( state_ = invalid ; state_ == invalid ; )
	{
		if( streams_.empty() ) state_ = end_of_stream ;
        else if( streams_[0]->get_state() == have_output ) state_ = have_output ;
        else 
        {
            merge_sensibly( foot_, streams_[0]->fetch_footer() ) ;
            streams_.pop_front() ;
            if( !streams_.empty() ) 
                merge_sensibly( hdr_, streams_[0]->fetch_header() ) ;
        }
	}
    return r ;
}

bool QualMasker::xform( Result& r )
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
			StreamWithProgress( std::auto_ptr< google::protobuf::io::ZeroCopyInputStream > is, const string& name, int64_t total )
				: is_( is ), name_( basename( name ) ), total_( total ), read_( 0 ) {}

			virtual bool Next( const void **data, int *size ) { return check( is_->Next( data, size ) ) ; }
			virtual void BackUp( int count ) { is_->BackUp( count ) ; }
			virtual bool Skip( int count ) { return check( is_->Skip( count ) ) ; }
			virtual int64_t ByteCount() const { return is_->ByteCount() ; }


	} ;

	bool magic( const void *p, int l, const char* sig )
	{
		for( const char* q = (const char*)p ; *sig ; ++q, --l, ++sig )
			if( !l || *q != *sig ) return false ;
		return true ;
	}
	bool is_crap( const void *p, int l ) 
	{
		return l >= 4 && *((const uint32_t*)p) == 0 ;
	}
	bool is_fastq( const void *p, int l ) 
	{
		const uint8_t* q = (const uint8_t*)p ;
		return l >= 3 && (*q == '>' || *q == '@') && isprint( q[1] ) && !magic( p, l, "@HD" ) ;
	}
	bool is_sam( const void *p, int l ) 
	{
		const uint8_t* q = (const uint8_t*)p ;
		if( l < 8 ) return false ;
		for( int i = 0 ; i != 8 ; ++i )
			if( !isprint( q[i] ) ) return false ;
		return true ;
	}

	StreamHolder make_input_stream_( std::auto_ptr< google::protobuf::io::ZeroCopyInputStream > is, const string& name, bool solexa_scores, char origin, const string& genome )
	{
		// peek into stream, but put it back.  then check magic numbers and
		// create the right stream
		const void* p ; int l ;
		if( !is->Next( &p, &l ) ) return new Stream() ; // gives an empty stream
		is->BackUp( l ) ;

		if( magic( p, l, "ANFO" ) ) {
			console.output( Console::info, name + ": linear ANFO file" ) ;
			return new AnfoReader( is, name ) ;
		}
		if( magic( p, l, "ANF1" ) ) {
			console.output( Console::info, name + ": chunked ANFO file" ) ;
			return new ChunkedReader( is, name ) ;
		}
		if( magic( p, l, ".sff" ) ) {
			console.output( Console::info, name + ": SFF file" ) ;
			return new SffReader( is, name ) ;
		}
		if( magic( p, l, "BAM\x01" ) ) {
			console.output( Console::info, name + ": BAM file" ) ;
			if( genome.empty() )
				console.output( Console::warning, name + ": no genome set for BAM parsing" ) ;
			return new BamReader( is, name, genome ) ;
		}
		if( magic( p, l, "BZh" ) )
		{
#if HAVE_LIBBZ2 && HAVE_BZLIB_H
			console.output( Console::info, name + ": BZip compressed" ) ;
			std::auto_ptr< google::protobuf::io::ZeroCopyInputStream > bs( new BunzipStream( is ) ) ;
			return make_input_stream_( bs, name, solexa_scores, origin, genome ) ;
#else
			throw "found BZip'ed file, but have no libbz2 support" ;
#endif
		}
		if( magic( p, l, "\x1f\x8b" ) )
		{
#if HAVE_LIBZ && HAVE_ZLIB_H
			console.output( Console::info, name + ": GZip compressed" ) ;
			std::auto_ptr< google::protobuf::io::ZeroCopyInputStream > zs( new InflateStream( is ) ) ;
			return make_input_stream_( zs, name, solexa_scores, origin, genome ) ;
#else
			throw "found GZip'ed file, but have no zlib support" ;
#endif
		}
		if( is_crap( p, l ) )
		{
			throw name + ": corrupt file" ;
		}
		if( is_fastq( p, l ) ) 
		{
			console.output( Console::info, name + ": probably FastA or FastQ" ) ;
			return new FastqReader( is, solexa_scores, origin ) ;
		}
		if( is_sam( p, l ) ) 
		{
			console.output( Console::notice, name + ": unknown, probably SAM" ) ;
			if( genome.empty() )
				console.output( Console::warning, name + ": no genome set for SAM parsing" ) ;
			return new SamReader( is, name, genome ) ;
		}
		throw name + ": weird unrecognized file" ;
	}
} ;

UniversalReader::UniversalReader(
				const std::string& name,
				google::protobuf::io::ZeroCopyInputStream* is,
				bool solexa_scores,
				int origin,
				const string& genome
				)
    : is_( is ), name_( name ), str_(), solexa_scores_( solexa_scores ), origin_( origin ), genome_( genome )
{
    // we cannot open the file just yet, but we can check if it is there
    if( !is_.get() ) throw_errno_if_minus1( access( name_.c_str(), R_OK ), "accessing ", name_.c_str() ) ;
}

Header UniversalReader::fetch_header() 
{
	if( !str_ ) {
		if( !is_.get() ) {
			int fd = throw_errno_if_minus1( open( name_.c_str(), O_RDONLY ), "opening (UniversalReader) ", name_.c_str() ) ;
			std::auto_ptr< google::protobuf::io::FileInputStream > s( new google::protobuf::io::FileInputStream( fd ) ) ;
			s->SetCloseOnDelete( true ) ;
			struct stat st ;
			if( fstat( fd, &st ) ) is_ = s ;
			else is_.reset( new StreamWithProgress( 
						std::auto_ptr<google::protobuf::io::ZeroCopyInputStream>(s), name_, st.st_size ) ) ;
		}
		str_ = make_input_stream_( is_, name_, solexa_scores_, origin_, genome_ ) ;
	}
	return str_->fetch_header() ;
}

uint8_t SffReader::read_uint8()
{
	while( !buf_size_ )
		if( !is_->Next( &buf_, &buf_size_ ) )
			throw "BamReader: premature end of file in " + name_ ;

	uint8_t x = *(const char*)buf_ ;
	buf_ = (const char*)buf_ + 1 ;
	--buf_size_ ;
	return x ;
}

uint16_t SffReader::read_uint16() { return ((uint16_t)read_uint8() << 8) | read_uint8() ; }
uint32_t SffReader::read_uint32() { return ((uint32_t)read_uint16() << 16) | read_uint16() ; }
void SffReader::read_string( unsigned l, string* s ) {
	s->clear() ;
	while( l ) {
		if( !buf_size_ && !is_->Next( &buf_, &buf_size_ ) )
			throw "BamReader: premature end of file in " + name_ ;

		while( buf_size_ && l ) 
		{
			s->push_back( *(char*)buf_ ) ;
			--l ;
			--buf_size_ ;
			buf_ = (const char*)buf_ + 1 ;
		}
	}
}
void SffReader::skip( int l ) 
{
	while( l ) {
		if( !buf_size_ && !is_->Next( &buf_, &buf_size_ ) )
			throw "BamReader: premature end of file in " + name_ ;
		
		unsigned k = min( l, buf_size_ ) ;
		l -= k ;
		buf_size_ -= k ;
		buf_ = (char*)buf_ + k ;
	}
}

Header SffReader::fetch_header()
{
	// common header as per section 13.3.8.1 of 454 manual
	if( read_uint32() != 0x2E736666 || read_uint32() != 1 ) 
		throw "SffReader: " + name_ + " has wrong magic number or version" ;
	skip( 12 ) ;										// index offset and length
	remaining_ = read_uint32() ;						// number of reads
	unsigned hdr_length = read_uint16() ;				// header length
	read_uint16() ;										// length of key (ignored, is trimmed)
	number_of_flows_ = read_uint16() ;
	if( read_uint8() != 1 )
		throw "SffReader: " + name_ + " has unknown flow format code" ;
	// rest is ignored, will lump it into the padding
	skip( hdr_length - 31 ) ;

	state_ = remaining_ ? have_output : end_of_stream ;
	hdr_.Clear() ;
	return hdr_ ;
}

Result SffReader::fetch_result()
{
	res_.Clear() ;

	// read header as per section 13.3.8.2 of 454 manual
	unsigned hdr_length = read_uint16() ;
	unsigned name_length = read_uint16() ;
	unsigned num_bases = read_uint32() ;
	unsigned clip_qual_left = read_uint16() ;
	unsigned clip_qual_right = read_uint16() ;
	unsigned clip_adapter_left = read_uint16() ;
	unsigned clip_adapter_right = read_uint16() ;
	read_string( name_length, res_.mutable_read()->mutable_seqid() ) ;		// read name
	skip( hdr_length - 16 - name_length ) ;									// padding

	// read data as per section 13.3.8.3 of 454 manual
	skip( number_of_flows_ * 2 ) ;											// flow values
	skip( num_bases ) ;														// flow indexes
	read_string( num_bases, res_.mutable_read()->mutable_sequence() ) ;		// bases
	read_string( num_bases, res_.mutable_read()->mutable_quality() ) ;		// quality
	skip (7 - ((number_of_flows_ * 2 + num_bases * 3 -1) % 8)) ;			// eight_byte_padding

	if( clip_qual_left || clip_adapter_left ) 
		res_.mutable_read()->set_trim_left( max( clip_qual_left, clip_adapter_left ) -1 ) ;

	if( clip_qual_right && clip_adapter_right ) 
		res_.mutable_read()->set_trim_right( res_.read().sequence().size() - min( clip_qual_right, clip_adapter_right ) ) ;
	else if( clip_qual_right )
		res_.mutable_read()->set_trim_right( res_.read().sequence().size() - clip_qual_right ) ;
	else if( clip_adapter_right )
		res_.mutable_read()->set_trim_right( res_.read().sequence().size() - clip_adapter_right ) ;

	if( !--remaining_ )
    {
        state_ = end_of_stream ;
        is_.reset(0) ;
    }
	return res_ ;
}


uint8_t BamReader::read_uint8()
{
	while( !buf_size_ )
		if( !is_->Next( &buf_, &buf_size_ ) )
			throw "BamReader: premature end of file in " + name_ ;

	uint8_t x = *(const char*)buf_ ;
	buf_ = (const char*)buf_ + 1 ;
	--buf_size_ ;
	return x ;
}

uint16_t BamReader::read_uint16() { return read_uint8() | ((uint16_t)read_uint8() << 8) ; }
uint32_t BamReader::read_uint32() { return read_uint16() | ((uint32_t)read_uint16() << 16) ; }
void BamReader::read_string( unsigned l, string* s ) {
	s->clear() ;
	while( l ) {
		if( !buf_size_ && !is_->Next( &buf_, &buf_size_ ) )
			throw "BamReader: premature end of file in " + name_ ;

		while( buf_size_ && l ) 
		{
			s->push_back( *(char*)buf_ ) ;
			--l ;
			--buf_size_ ;
			buf_ = (const char*)buf_ + 1 ;
		}
	}
}

BamReader::BamReader( auto_ptr< google::protobuf::io::ZeroCopyInputStream > is, const string& name, const string& genome ) : 
	is_( is ), name_( name ), genome_( genome ), buf_(0), buf_size_(0)
{
	if( read_uint32() != 0x014d4142 ) throw "missing BAM magic number" ;
	string sam_header ;
	read_string( read_uint32(), &sam_header ) ;
	unsigned nrefseq = read_uint32() ;
	refseqs_.resize( nrefseq ) ;
	for( unsigned i = 0 ; i != nrefseq ; ++i )
	{
		read_string( read_uint32()-1, &refseqs_[i] ) ;
		read_uint8() ; // NUL terminator
		read_uint32() ; // chromosome length
	}
	while( buf_ && !buf_size_ ) if( !is_->Next( &buf_, &buf_size_ ) ) buf_ = 0 ;
	state_ = buf_ && buf_size_ ? have_output : end_of_stream ;
}

Result BamReader::fetch_result()
{
	res_.Clear() ;
	Hit *h = res_.add_hit() ;

	unsigned bsize = read_uint32() ;
	unsigned refid = read_uint32() ;
	if( refid < refseqs_.size() ) h->set_sequence( refseqs_[ refid ] ) ;
	if( !genome_.empty() ) h->set_genome_name( genome_ ) ;
	h->set_start_pos( read_uint32() ) ;
	int namelen = read_uint8() ;
	int mapq = read_uint8() ;
	if( mapq < 255 ) h->set_diff_to_next( mapq ) ;
	read_uint16() ; // BIN
	int cigarlen = read_uint16() ;
	int flags = read_uint16() ;
	int seqlen = read_uint32() ;
	read_uint32() ; // MATE RID
	read_uint32() ; // MATE POS
	read_uint32() ; // INS SIZE
	read_string( namelen-1, res_.mutable_read()->mutable_seqid() ) ;
	read_uint8() ;
	unsigned aln_len = 0 ;
	for( int i = 0 ; i != cigarlen ; ++i ) {
		int cig = read_uint32() ;
		h->add_cigar( cig ) ;
		switch( cigar_op( cig ) )
		{
			case Hit::Match:
			case Hit::Mismatch:
			case Hit::Delete:
			case Hit::Skip:
				aln_len += cigar_len( cig ) ;
				break ;

			case Hit::Insert:
			case Hit::SoftClip:
			case Hit::HardClip:
			case Hit::Pad:
				break ;
		}
	}
	h->set_aln_length( flags & 0x10 ? -aln_len : aln_len ) ;
	if( flags & 0x10 ) reverse( h->mutable_cigar()->begin(), h->mutable_cigar()->end() ) ;

	if( flags & 0x10 ) {
		string seq ;
		for( int i = 0 ; i != seqlen ; ++i )
		{
			static char bases_rev[] = "=TG.C...A......N" ;
			uint8_t code = read_uint8() ;
			seq.push_back( bases_rev[ (code >> 4) & 0xf ] ) ;
			if( ++i == seqlen ) break ;
			seq.push_back( bases_rev[ code & 0xf ] ) ;
		}
		res_.mutable_read()->mutable_sequence()->assign( seq.rbegin(), seq.rend() ) ;
	}
	else {
		for( int i = 0 ; i != seqlen ; ++i )
		{
			static char bases_fwd[] = "=AC.G...T......N" ;
			uint8_t code = read_uint8() ;
			res_.mutable_read()->mutable_sequence()->push_back( bases_fwd[ (code >> 4) & 0xf ] ) ;
			if( ++i == seqlen ) break ;
			res_.mutable_read()->mutable_sequence()->push_back( bases_fwd[ code & 0xf ] ) ;
		}
	}

	string qual ;
	read_string( seqlen, &qual ) ;
	bool has_qual = false ;
	for( unsigned i = 0 ; i != qual.size() ; ++i )
		if( (uint8_t)qual[i] != 0xff ) has_qual = true ;
	
	if( has_qual ) {
		if( flags & 0x10 ) res_.mutable_read()->mutable_quality()->assign(
				qual.rbegin(), qual.rend() ) ;
		else res_.mutable_read()->set_quality( qual ) ;
	}

	bsize -= 32 + namelen + 4*cigarlen + ((seqlen+1)/2) + seqlen ;
	while( bsize ) 
	{
		char tag0 = read_uint8() ;
		char tag1 = read_uint8() ;
		char type = read_uint8() ;
		bsize -= 3 ;

		int ival = 0 ;
		string sval ;

		switch( type )
		{
			case 'A':
			case 'c': case 'C': ival = read_uint8() ; bsize -= 1 ; break ;
			case 's': case 'S': ival = read_uint16() ; bsize -= 2 ; break ;
			case 'i': case 'I': ival = read_uint32() ; bsize -= 4 ; break ;
			case 'f': throw "can't decode float" ;
			case 'Z': --bsize ; while( char c = read_uint8() ) { sval.push_back( c ) ; --bsize ; } break ;
			case 'H': throw "can't decode HEX string" ;
			default: throw string( "can't decode shit " ) + tag0 + tag1 + ':' + type ;
		}
		if( tag0 == 'U' && tag1 == 'Q' ) h->set_score( ival ) ;
		if( tag0 == 'A' && tag1 == 'S' && !h->has_score() ) h->set_score( ival ) ;
	}
	if( flags & 4 ) res_.mutable_hit()->RemoveLast() ;

	while( buf_ && !buf_size_ ) if( !is_->Next( &buf_, &buf_size_ ) ) buf_ = 0 ;
	state_ = buf_ && buf_size_ ? have_output : end_of_stream ;
	return res_ ;
}

} ; // namespace

std::pair< PipeInputStream*, std::string > make_PipeInputStream( const std::string& p )
{
	console.output( Console::notice, "piping from " + p ) ;

	int fds[2] ;
	throw_errno_if_minus1( pipe( fds ), "creating pipe" ) ;

	pid_t chld = throw_errno_if_minus1( fork(), "forking pipe process" ) ;
	if( chld == 0 ) {
		throw_errno_if_minus1( dup2( fds[1], 1 ), "duplicating file descriptor" ) ;
		if( fds[1] != 1 ) throw_errno_if_minus1( close( fds[1] ), "closing fd" ) ;
		throw_errno_if_minus1( close( fds[0] ), "closing fd" ) ;
		const char *c = p.c_str() ;
		while( *c && isspace( *c ) ) ++c ;
		execl( "/bin/sh", "sh", "-c", c, (char*)0 ) ;
	}

	throw_errno_if_minus1( close( fds[1] ), "closing fd" ) ;
	return std::make_pair( new PipeInputStream( fds[0], chld ), "<pipe>" ) ;
}

std::pair< PipeOutputStream*, std::string > make_PipeOutputStream( const std::string& p )
{
	console.output( Console::notice, "piping to " + p ) ;

	int fds[2] ;
	throw_errno_if_minus1( pipe( fds ), "creating pipe" ) ;

	pid_t chld = throw_errno_if_minus1( fork(), "forking pipe process" ) ;
	if( chld == 0 ) {
		throw_errno_if_minus1( dup2( fds[0], 0 ), "duplicating file descriptor" ) ;
		if( fds[0] != 0 ) throw_errno_if_minus1( close( fds[0] ), "closing fd" ) ;
		throw_errno_if_minus1( close( fds[1] ), "closing fd" ) ;
		const char *c = p.c_str() ;
		while( *c && isspace( *c ) ) ++c ;
		execl( "/bin/sh", "sh", "-c", c, (char*)0 ) ;
	}

	throw_errno_if_minus1( close( fds[0] ), "closing fd" ) ;
	return std::make_pair( new PipeOutputStream( fds[1], chld ), "<pipe>" ) ;
}

zero_copy_output_buf::~zero_copy_output_buf() { sync() ; }

int zero_copy_output_buf::sync()
{
	if( epptr() != pptr() ) {
		os_->BackUp( epptr() - pptr() ) ;
		setp( 0, 0 ) ;
	}
	return 0 ;
}

zero_copy_output_buf::int_type zero_copy_output_buf::overflow( zero_copy_output_buf::int_type c )
{
	sync() ;
	void *buf ;
	int len ;
	if( !os_->Next( &buf, &len ) ) return traits_type::eof() ;
	setp( static_cast<char_type*>( buf ), static_cast<char_type*>( buf ) + len ) ;
	return traits_type::eq_int_type( c, traits_type::eof() )
		? traits_type::not_eof( c ) : sputc( c ) ;
}

