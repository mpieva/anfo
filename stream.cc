#include "compress_stream.h"
#include "stream.h"
#include "util.h"

#include <google/protobuf/repeated_field.h>

#include <iostream>
#include <set>

#if HAVE_FCNTL_H
#include <fcntl.h>
#endif

namespace std {
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

int AnfoReader::num_files_ = 0 ;

AnfoReader::AnfoReader( const std::string& name, bool quiet )
	: iis_( throw_errno_if_minus1( open( name.c_str(), O_RDONLY ), "opening ", name.c_str() ) )
	, zis_( decompress( &iis_ ) ), name_( name ), quiet_( quiet )
{
	initialize() ;
}

AnfoReader::AnfoReader( int fd, const std::string& name, bool quiet )
	: iis_( fd ), zis_( decompress( &iis_ ) ), name_( name ), quiet_( quiet )
{
	initialize() ;
}

void AnfoReader::initialize()
{
	iis_.SetCloseOnDelete( true ) ;
	std::string magic ;
	CodedInputStream cis( zis_.get() ) ;

	if( !cis.ReadString( &magic, 4 ) || magic != "ANFO" ) {
		if( !quiet_ ) clog << "\033[K" << name_ << ": not an ANFO file" << endl ;
	} else {
		uint32_t tag ;
		if( cis.ReadVarint32( &tag ) && tag == 10 && cis.ReadVarint32( &tag ) ) {
			int lim = cis.PushLimit( tag ) ;
			if( hdr_.ParseFromCodedStream( &cis ) ) {
				sanitize( hdr_ ) ;
				cis.PopLimit( lim ) ;
				++num_files_ ;
				read_next_message( cis ) ;
				return ;
			}
		}
		if( !quiet_ ) clog << "\033[K" << name_ << ": deserialization error in header" << endl ;
	}
	foot_.set_exit_code(1) ;
}

Result AnfoReader::fetch_result()
{
	Result r ;
	swap( r, res_ ) ;
	CodedInputStream cis( zis_.get() ) ;
	read_next_message( cis ) ;
	return r ;
}

void AnfoReader::read_next_message( CodedInputStream& cis )
{
	state_ = invalid ;
	uint32_t tag = 0 ;
	if( cis.ExpectAtEnd() ) {
		if( !quiet_ ) clog << "\033[K" << name_ << ": unexpected end of stream" << endl ;
	}
	else if( (tag = cis.ReadTag()) )
	{
		uint32_t size = 0 ;
		if( cis.ReadVarint32( &size ) ) 
		{
			int lim = cis.PushLimit( size ) ;
			if( tag == 18 && res_.ParseFromCodedStream( &cis ) )
			{
				cis.PopLimit( lim ) ;
				state_ = have_output ;
				return ;
			}
			if( tag == 26 && foot_.ParseFromCodedStream( &cis ) )
			{
				cis.PopLimit( lim ) ;
				state_ = end_of_stream ;
				return ; 
			}

			if( !quiet_ ) clog << "\033[K" << name_ << ": deserialization error" << endl ;
		}
	}
}

AnfoWriter::AnfoWriter( ZeroCopyOutputStream *zos ) : o_( zos )
{
	o_.WriteRaw( "ANFO", 4 ) ;
}

AnfoWriter::AnfoWriter( int fd, bool expensive )
	: fos_( new FileOutputStream( fd ) )
	, zos_( expensive ? compress_small( fos_.get() ) : compress_fast(  fos_.get() ) )
	, o_( zos_.get() )
{
	o_.WriteRaw( "ANFO", 4 ) ;
}

AnfoWriter::AnfoWriter( const char* fname, bool expensive )
	: fos_( new FileOutputStream( throw_errno_if_minus1( creat( fname, 0777 ), "opening", fname ) ) )
	, zos_( expensive ? compress_small( fos_.get() ) : compress_fast(  fos_.get() ) )
	, o_( zos_.get() )
{
	fos_->SetCloseOnDelete( true ) ;
	o_.WriteRaw( "ANFO", 4 ) ;
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
	bool no_task_id = !lhs.has_sge_task_id() || (rhs.has_sge_task_id() && lhs.sge_task_id() != rhs.sge_task_id()) ;
	bool no_job_id = !lhs.has_sge_job_id() || (rhs.has_sge_job_id() && lhs.sge_job_id() != rhs.sge_job_id()) ;

	lhs.MergeFrom( rhs ) ;
	nub( *lhs.mutable_sge_slicing_index() ) ;
	nub( *lhs.mutable_command_line() ) ;
	nub( *lhs.mutable_config()->mutable_genome_path() ) ;
	nub( *lhs.mutable_config()->mutable_policy() ) ;

	if( no_task_id ) lhs.clear_sge_task_id() ;
	if( no_job_id ) lhs.clear_sge_job_id() ;
	sanitize( lhs ) ;
}

void sanitize( Header& hdr )
{
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
}

//! \brief merges two results, keeping the best hit
void merge_sensibly( Result& lhs, const Result& rhs )
{
	// How to merge what...
	// - seqid, description, sequence, trimpoints: all equal, no merging needed

	// - reason: do the sensible thing...
	if( lhs.reason() != aligned && rhs.reason() != no_policy && rhs.reason() != no_seeds )
	{
		if( rhs.reason() == aligned ) lhs.clear_reason() ;
		else if( lhs.reason() == no_seeds || lhs.reason() == no_policy ) lhs.set_reason( rhs.reason() ) ;
		else if( lhs.reason() == too_many_seeds && rhs.reason() == bad_alignment ) lhs.set_reason( rhs.reason() ) ;
	}

	// - num_xxx: just add them
	lhs.set_num_raw_seeds( lhs.num_raw_seeds() + rhs.num_raw_seeds() ) ;
	lhs.set_num_grown_seeds( lhs.num_grown_seeds() + rhs.num_grown_seeds() ) ;
	lhs.set_num_clumps( lhs.num_clumps() + rhs.num_clumps() ) ;

	// - best_hit, diff_to_next_species, diff_to_next_order: TODO
	
	// - best_to_genome: take better hit, recalculate diff_to_next{,_chromosome{,_class}}
	if( rhs.has_best_to_genome() )
	{
		if( lhs.has_best_to_genome() )
		{
			// two hits, this is work...
			if( lhs.best_to_genome().score() <= rhs.best_to_genome().score() )
			{
				// left is better
				if( !lhs.has_diff_to_next() ||
						lhs.best_to_genome().score() + lhs.diff_to_next() > rhs.best_to_genome().score() )
					lhs.set_diff_to_next( rhs.best_to_genome().score() - lhs.best_to_genome().score() ) ;

				//! \todo diff to chromosome, chromosome class? dunno...
			}
			else
			{
				// right is better
				if( !rhs.has_diff_to_next() ||
						rhs.best_to_genome().score() + rhs.diff_to_next() > lhs.best_to_genome().score() )
					lhs.set_diff_to_next( lhs.best_to_genome().score() - rhs.best_to_genome().score() ) ;

				*lhs.mutable_best_to_genome() = rhs.best_to_genome() ;
				//! \todo diff to chromosome, chromosome class? dunno...
			}
		}
		else
		{
			// no hit at left side --> just assign
			*lhs.mutable_best_to_genome() = rhs.best_to_genome() ;
			if( rhs.has_diff_to_next() )
				lhs.set_diff_to_next( rhs.diff_to_next() ) ;

			if( rhs.has_diff_to_next_chromosome() )
				lhs.set_diff_to_next_chromosome( rhs.diff_to_next_chromosome() ) ;

			if( rhs.has_diff_to_next_chromosome_class() )
				lhs.set_diff_to_next_chromosome_class( rhs.diff_to_next_chromosome_class() ) ;
		}
	}
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

unsigned len_from_bin_cigar( const string& cigar )
{
	unsigned l = 0 ;
	for( size_t i = 0 ; i != cigar.size() ; ++i )
	{
		if( (uint8_t)cigar[i] < 128 ) l += (uint8_t)cigar[i] ;
		else if( (uint8_t)cigar[i] < 192 ) l += (unsigned)(uint8_t)cigar[i] - 128 ;
	}
	return l ;
}

bool ScoreFilter::xform( Result& r ) {
	if( genome_ && r.has_best_to_genome() && r.best_to_genome().genome_name() == genome_ &&
			r.best_to_genome().score() > slope_ * (len_from_bin_cigar( r.best_to_genome().cigar() ) - intercept_) )
	{
		r.set_reason( bad_alignment ) ;
		r.clear_best_to_genome() ;
	}
	else if( !genome_ && r.has_best_hit() &&
			r.best_hit().score() > slope_ * (len_from_bin_cigar( r.best_hit().cigar() ) - intercept_) )
	{
		r.set_reason( bad_alignment ) ;
		r.clear_best_hit() ;
	}
	return true ;
}

bool LengthFilter::xform( Result& r ) {
	int len = ( r.has_trim_right() ? r.trim_right() : r.sequence().size() )
		    - ( r.has_trim_left() ? r.trim_left() : 0 ) ;
	if( r.has_best_to_genome() && len < minlength_ )
	{
		r.set_reason( no_policy ) ;
		r.clear_best_to_genome() ;
	}
	if( r.has_best_hit() && len < minlength_ )
	{
		r.set_reason( no_policy ) ;
		r.clear_best_hit() ;
	}
	return true ;
}

bool HitFilter::xform( Result& r ) 
{
	if( g_ && r.has_best_to_genome() && r.best_to_genome().genome_name() == g_ ) return true ; // XXX which genome?
	if( !g_ && r.has_best_hit() ) return true ;
	return false ;
}

//! \todo configure which genome we're interested in
bool RmdupStream::is_duplicate( const output::Result& lhs, const output::Result& rhs ) 
{
	if( !lhs.has_best_to_genome() || !rhs.has_best_to_genome()
			|| lhs.sequence().size() != rhs.sequence().size() ) return false ;
	const output::Hit &l = lhs.best_to_genome(), &r = rhs.best_to_genome() ;

	return l.genome_name() == r.genome_name() && l.sequence() == r.sequence()
		&& l.start_pos() == r.start_pos() && l.aln_length() == r.aln_length() ;
}

//! \todo How do we deal with ambiguity codes?  What's the meaning of
//!       their quality scores anyway?
void RmdupStream::add_read( const Result& rhs ) 
{
	Read *rd = cur_.add_member() ;
	if( rhs.has_seqid() ) rd->set_seqid( rhs.seqid() ) ;
	if( rhs.has_description() ) rd->set_description( rhs.description() ) ;
	if( rhs.has_sequence() ) rd->set_sequence( rhs.sequence() ) ;
	if( rhs.has_quality() ) rd->set_quality( rhs.quality() ) ;

	for( size_t i = 0 ; i != rhs.sequence().size() ; ++i )
	{
		int base = -1 ;
		switch( rhs.sequence()[i] ) {
			case 'a': case 'A': base = 0 ; break ;
			case 'c': case 'C': base = 1 ; break ;
			case 't': case 'T':
			case 'u': case 'U': base = 2 ; break ;
			case 'g': case 'G': base = 3 ; break ;
		}

		Logdom qual = Logdom::from_phred( rhs.has_quality() ? rhs.quality()[i] : 30 ) ;
		for( int j = 0 ; j != 4 ; ++j )
			// XXX distribute errors sensibly
			quals_[j].at(i) *= j != base 
				? qual / Logdom::from_float( 3 ) 
				: Logdom::from_float(1) - qual ;
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

//! \brief signals end of input stream
//! After the footer no more input is possible.  We were waiting for
//! input, so no output was available.  cur_ might contain valid data,
//! and if so, we call a consensus and offer it as output.  Else we
//! signal end of stream.
void RmdupStream::put_footer( const Footer& f ) { 
	foot_ = f ;
	if( cur_.IsInitialized() ) 
	{
		state_ = have_output ;
		call_consensus() ;
		swap( cur_, res_ ) ;
	}
	else state_ = end_of_stream ;
}

//! \brief receives a result record and merges it if appropriate
//! \todo We do not want to merge sequences that have a bad alignment
//!       score, instead they should pass through without disturbing the
//!       merging process.
//!
//! There are the following possibilities what to do here:
//! - A result with a bad alignment is passed through (means it is
//!   stored in res_ and output becomes available). XXX
//! - If cur_ is invalid, the result is stored there and cur_ becomes
//!   valid.
//! - A result with correct coordinates (according to is_duplicate) is
//!   directly merged into cur_, no output becomes available.
//! - Anything else cannot be merged, so a consensus is called, cur_
//!   moves to res_, next moves to cur_, and output becomes available.
void RmdupStream::put_result( const Result& next ) 
{
	if( 0 /* check for bad alignment */ ) 
	{
		res_ = next ;
		state_ = have_output ;
	}
	else if( !cur_.IsInitialized() ) {
		cur_ = next ;
	}
	else if( is_duplicate( cur_, next ) )
	{
		// Merge them.  If cur is a plain result, turn it into a
		// degenerate merged one first...
		if( cur_.member_size() == 0 )
		{
			for( size_t i = 0 ; i != 4 ; ++i )
			{
				quals_[i].clear() ;
				quals_[i].resize( cur_.sequence().size() ) ;
			}

			add_read( cur_ ) ;
			cur_.set_seqid( "C_" + cur_.seqid() ) ;
			cur_.clear_description() ;
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
		state_ = have_output ;
	}
}

void RmdupStream::call_consensus()
{
	if( !cur_.member_size() ) return ;

	cur_.clear_sequence() ;
	cur_.clear_quality() ;
	for( size_t i = 0 ; i != quals_[0].size() ; ++i )
	{
		// select base with highest quality
		size_t m = 0 ;
		for( size_t j = 1 ; j != 4 ; ++j )
			if( quals_[j].at( i ) > quals_[m].at( i ) ) m = j ;

		// but calculate the error probability from _all_others_ to
		// retain precision

		Logdom num = quals_[m?0:1].at(i) ;
		Logdom denom = quals_[0].at(i) ;
		for( size_t j = 1 ; j != 4 ; ++j )
		{
			denom += quals_[j].at( i ) ;
			if( j != m ) num += quals_[j].at( i ) ;
		}

		int qscore = (num/denom).to_phred() ;
		if( qscore > 127 ) qscore = 127 ;

#if 0
		// debug code
		for( size_t j = 0 ; j != 4 ; ++j )
			std::cerr << ' ' << (j==m ? '*' : ' ') << ' ' << quals_[j].at( i ).to_float() << std::endl ;
		std::cerr << num.to_float() << '/' << denom.to_float() 
			<< " == " << (num/denom).to_phred() << std::endl ;
#endif

		cur_.mutable_sequence()->push_back( m["ACTG"] ) ;
		cur_.mutable_quality()->push_back( qscore ) ;
	}
}

void ConcatStream::add_stream( Stream* s )
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
		merge_sensibly( foot_, streams_[0]->fetch_footer() ) ;
		delete s ;
	}
}

Result ConcatStream::fetch_result()
{
	Result r = streams_[0]->fetch_result() ;
	if( streams_[0]->get_state() != have_output )
	{
		merge_sensibly( foot_, streams_[0]->fetch_footer() ) ;
		delete streams_[0] ;
		streams_.pop_front() ;
	}
	if( streams_.empty() ) state_ = end_of_stream ;
	return r ;
}

bool QualFilter::xform( Result& r )
{
	if( r.has_quality() ) 
		for( size_t i = 0 ; i != r.sequence().size() && i != r.quality().size() ; ++i )
			if( r.quality()[i] < q_ ) (*r.mutable_sequence())[i] = '-' ;
	return true ;
}


} ; // namespace
