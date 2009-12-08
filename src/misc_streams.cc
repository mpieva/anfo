#if HAVE_CONFIG_H
#include "../config.h"
#endif

#include "misc_streams.h"

#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/text_format.h>

#include <algorithm>
#include <iostream>

namespace streams {

using namespace google::protobuf ;
using namespace google::protobuf::io ;

Result MergeStream::fetch_result() 
{
	assert( state_ == have_output ) ;
	assert( rs_.size() == streams_.size() ) ;
	assert( streams_.size() != 0 ) ;

	int min_idx = 0 ;
	for( size_t i = 1 ; i != rs_.size() ; ++i ) 
		if( ( mode_ == by_coordinate && by_genome_coordinate( gs_ )( &rs_[ i ], &rs_[ min_idx ] ) )
				|| ( mode_ == by_name && by_seqid()( &rs_[ i ], &rs_[ min_idx ] ) ) )
			min_idx = i ;

	Result res = rs_[ min_idx ] ;
	Stream &sm = *streams_[ min_idx ] ;
	if( sm.get_state() == have_output )
	{
		rs_[ min_idx ] = sm.fetch_result() ; 
	}
	else
	{
		merge_sensibly( foot_, sm.fetch_footer() ) ;
		streams_.erase( streams_.begin() + min_idx ) ;
		rs_.erase( rs_.begin() + min_idx ) ;
		if( rs_.empty() ) state_ = end_of_stream ;
	}
	return res ;
}

//! Checks if at least as many inputs are known as needed.  It assumes
//! only one genome index was used and takes the number of slices
//! declared for it.  (This is mostly useful in MegaMergeStream to clean
//! up the mess left by a heavily sliced grid job.  Try to avoid using
//! this 'feature'.)
bool BestHitStream::enough_inputs() const
{
	if( hdr_.config().policy_size() == 0 ||
			hdr_.config().policy(0).use_compact_index_size() == 0 ) return false ;
	if( !hdr_.config().policy(0).use_compact_index(0).has_number_of_slices() ) return true ;
	return nstreams_ >= hdr_.config().policy(0).use_compact_index(0).number_of_slices() ;
}

//! Enough information is read until one sequence has been reported on
//! by each input stream.  Everything else is buffered if necessary.
//! This works fine as long as the inputs come in in roughly the same
//! order and are complete.
Result BestHitStream::fetch_result()
{
	while( !streams_.empty() ) {
		if( cur_input_ == streams_.end() ) cur_input_ = streams_.begin() ;

		Stream &s = **cur_input_ ;
		assert( s.get_state() == have_output ) ;

		Result r = s.fetch_result() ;
		if( s.get_state() == end_of_stream )
		{
			merge_sensibly( foot_, s.fetch_footer() ) ;
			cur_input_ = streams_.erase( cur_input_ ) ;
		}
		else ++cur_input_ ;

		pair< size_t, Result > &p = buffer_[ r.read().seqid() ] ;
		++nread_ ;
		++p.first ;
		if( p.second.IsInitialized() ) merge_sensibly( p.second, r ) ;
		else p.second = r ;

		if( p.first == nstreams_ )
		{
			++nwritten_ ;
			Result r1 ;
			swap( r1, p.second ) ;
			buffer_.erase( r.read().seqid() ) ;

			if( buffer_.empty() && streams_.empty() )
				state_ = end_of_stream ;

			if( ((nread_ + nwritten_) & 0xFFFF) == 0 )
			{
				stringstream s ;
				s << "BestHitStream: " << nread_ << "in "
				<< nwritten_ << "out "
				<< buffer_.size() << "buf" ;
				progress_( Console::info, s.str() ) ;
			}
			return r1 ;
		}
	}

	Result r1 ;
	swap( r1, buffer_.begin()->second.second ) ;
	buffer_.erase( buffer_.begin()->first ) ;
	if( buffer_.empty() ) state_ = end_of_stream ;
	return r1 ;
}

unsigned SortingStream__ninstances = 0 ;

void MegaMergeStream::add_stream( StreamHolder s ) 
{
	Header h = s->fetch_header() ;
	if( h.has_sge_slicing_stride() )
	{
		Holder< BestHitStream > bhs = stream_per_slice_[ h.sge_slicing_index(0) ] ;
		if( !bhs ) bhs = new BestHitStream ;
		bhs->add_stream( s ) ;

		if( bhs->enough_inputs() ) {
			std::stringstream s ;
			s << "MegaMergeStream: got everything for slice " << h.sge_slicing_index(0) ;
			console.output( Console::info, s.str() ) ;
			stream_per_slice_.erase( h.sge_slicing_index(0) ) ;
			ConcatStream::add_stream( bhs ) ;
		}
	}
	else ConcatStream::add_stream( s ) ;
}

Header MegaMergeStream::fetch_header()
{
	if( !stream_per_slice_.empty() ) 
		console.output( Console::warning, "MegaMergeStream: input appears to be incomplete" ) ;

	for( map< int, Holder< BestHitStream > >::iterator
			l = stream_per_slice_.begin(),
			r = stream_per_slice_.end() ; l != r ; ++l )
		ConcatStream::add_stream( l->second ) ;

	return ConcatStream::fetch_header() ;
}


void RepairHeaderStream::put_header( const Header& h ) 
{
	char tmpname[] = "/tmp/anfo_header_XXXXXX" ;
	int fd = mkstemp( tmpname ) ;
	{
		FileOutputStream fos( fd ) ;
		TextFormat::Print( h, &fos ) ;
	}

	string cmd = (editor_.empty() ? getenv("EDITOR") : editor_) + " " + tmpname ;
	for(;;) {
		if( system( cmd.c_str() ) ) break ;;
		lseek( fd, 0, SEEK_SET ) ;
		FileInputStream fis( fd ) ;
		if( TextFormat::Parse( &fis, &hdr_ ) ) break ;
	} 
	throw_errno_if_minus1( unlink( tmpname ), "unlinking", tmpname ) ;
	state_ = need_input ;
}

void FanOut::put_header( const Header& h )
{
	for( citer i = streams_.begin() ; i != streams_.end() ; ++i )
		(*i)->put_header( h ) ;
	state_ = streams_.front()->get_state() ;
}

void FanOut::put_result( const Result& r )
{
	for( citer i = streams_.begin() ; i != streams_.end() ; ++i )
		(*i)->put_result( r ) ;
	state_ = streams_.front()->get_state() ;
}

void FanOut::put_footer( const Footer& f )
{
	Stream::put_footer( f ) ;
	for( citer i = streams_.begin() ; i != streams_.end() ; ++i )
		(*i)->put_footer( f ) ;
}

Footer FanOut::fetch_footer()
{
	for( citer i = streams_.begin() ; i != streams_.end() ; ++i )
		merge_sensibly( foot_, (*i)->fetch_footer() ) ;
	return Stream::fetch_footer() ;
}
	
void Compose::put_header( const Header& h_ )
{
	Header h = h_ ;
	citer i = streams_.begin(), e = streams_.end() ;
	for( e-- ; i != e ; ++i )
	{
		(*i)->put_header( h ) ;
		h = (*i)->fetch_header() ;
	}
	(*i)->put_header( h ) ;
}

Header Compose::fetch_header() 
{
	citer i = streams_.begin(), e = streams_.end() ;
	Header h = (*i)->fetch_header() ;
	for( ++i ; i != e ; ++i )
	{
		(*i)->put_header( h ) ;
		h = (*i)->fetch_header() ;
	}
	return h ;
}

// note weird calls to base(): makes a weird compiler happy...
Stream::state Compose::get_state()
{
	// look at a stream at a time, starting from the end
	for( criter i = streams_.rbegin() ;; )
	{
		// if we fell of the far end, we need more input
		if( i.base() == streams_.rend().base() ) { return need_input ; }
		// else consider the state
		state s = (*i)->get_state() ;
		if( s == need_input ) {
			// input comes from previous stream
			++i ;
		}
		else if( s == have_output ) {
			// output's available, either to the outside or to
			// downstream filters
			if( i.base() == streams_.rbegin().base() ) { return have_output ; }
			Result r = (*i)->fetch_result() ;
			--i ;
			(*i)->put_result( r ) ;
		}
		else if( s == end_of_stream ) {
			// nothing left, pass the footer to see if some data is
			// buffered
			if( i.base() == streams_.rbegin().base() ) { return end_of_stream ; }
			Footer f = (*i)->fetch_footer() ;
			--i ;
			(*i)->put_footer( f ) ;
		}
		else {
			// easy: something's broken
			return invalid ;
		}
	}
}

//! \brief prints statistics
//! The fields are:
//! -# some read name (as an audit trail)
//! -# total number of reads
//! -# number of mapped reads
//! -# number of uniquely mapped reads
//! -# number of unique sequences mapped uniquely
//! -# GC content in raw data
//! -# GC content in mapped data
//! -# total number of bases (== first moment of length dist.)
//! -# sum of squared lengths (== second moment of length dist.)
//! -# number of mapped bases (== first moment of mapped length dist.)
//! -# sum of squared lengths of mapped reads (== second moment of
//!    mapped length dist.)
//! 
//! Note that given first (m1) and second moment (m2) of a sample of size
//! n, the average is simply
//! \f[ m1 / n \f] 
//! and the variance is
//! \f[ (m2 - m1^2 / n) / (n-1) \f]
//! Also, the mean of the weighted contig length (the median of this
//! value would be the N50) is 
//! \f[ m2 / m1 \f]

void StatStream::printout( ostream& s, bool h )
{
	s << (h?"name\t#total\t#mapped\t#mapuniq\t#distinct\t"
		    "GCraw\tGCmapped\trawbases\trawbases^2\t"
			"mapbases\tmapbases^2\n":"")
	  << name_ << '\t' << total_ << '\t' 							// arbitrary name, reads
	  << mapped_ << '\t' << mapped_u_ << '\t' 						// mapped reads, uniquely mapped reads
	  << different_ << '\t'
      << 100*(float)bases_gc_ / (float)bases_ << '\t'				// GC content
	  << 100*(float)bases_gc_m_ / (float)bases_m_ << '\t'			// GC content, mapped only
	  << bases_  << '\t' << bases_squared_ << '\t'					// number of bases & 2nd moment
	  << bases_m_ << '\t' << bases_m_squared_ << endl ;             // number of mapped bases & 2nd moment
} 

void StatStream::put_result( const Result& r )
{
	unsigned bases = r.read().sequence().size() ;
	unsigned gc = count( r.read().sequence().begin(), r.read().sequence().end(), 'G' )
		        + count( r.read().sequence().begin(), r.read().sequence().end(), 'C' ) ;
	unsigned count = std::max( r.member_size(), 1 ) ;
	total_ += count ;
	bases_ += bases ;
	bases_gc_ += gc ;
	bases_squared_ += bases*bases ;
	if( const Hit* h = hit_to( r ) )
	{
		mapped_ += count ;
		bases_m_ += bases ;
		bases_m_squared_ += bases*bases ;
		bases_gc_m_ += gc ;
		if( !h->has_diff_to_next() || h->diff_to_next() >= 60 )
		{
			mapped_u_ += count ;
			++different_ ;
		}
	}
	if( name_.empty() ) name_ = r.read().seqid() ;
}

void StatStream::put_footer( const Footer& f )
{
	Stream::put_footer( f ) ;
	if( fn_.empty() || fn_ == "-" ) printout( cout, true ) ;
	else if( fn_ == "+-" ) printout( cout, false ) ;
	else if( fn_[0] == '+' ) 
	{
		ofstream s( fn_.substr(1).c_str(), ios_base::app ) ;
		printout( s, false ) ;
	}
	else 
	{
		ofstream s( fn_.c_str(), ios_base::trunc ) ;
		printout( s, true ) ;
	}
}

RegionFilter::Regions3 RegionFilter::all_regions ;

RegionFilter::RegionFilter( const pair< istream*, string >& p )
{
	auto_ptr< istream > deffile( p.first ) ;
	const string &fn = p.second ;

	Regions3::iterator i = all_regions.find( fn ) ;
	if( i != all_regions.end() ) my_regions = &i->second ;
	else
	{
		my_regions = &all_regions[ fn ] ;

		console.output( Console::notice, "Loading annotation from " + fn ) ;
		std::string junk, chrom ;
		unsigned start, end ;

		std::getline( *deffile, junk ) ;
		while( std::getline( *deffile >> junk >> chrom >> start >> end, junk ) )
			(*my_regions)[ chrom ][ end ] = start ;
	}
}

void Sanitizer::put_header( const Header& h )
{
	Filter::put_header( h ) ;
	hdr_.clear_sge_slicing_index() ;
	hdr_.clear_sge_slicing_stride() ;
	hdr_.clear_command_line() ;
	hdr_.clear_sge_job_id() ;
	hdr_.clear_sge_task_id() ;
	hdr_.mutable_config()->clear_genome_path() ;
}

bool Sanitizer::xform( Result& r )
{
	r.clear_aln_stats() ;
	return true ;
}

} ; // namespace

