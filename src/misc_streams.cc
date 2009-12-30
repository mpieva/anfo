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
	Stream::put_header( h ) ;
	for( citer i = streams_.begin() ; i != streams_.end() ; ++i )
		(*i)->put_header( h ) ;
}

void FanOut::put_result( const Result& r )
{
	for( citer i = streams_.begin() ; i != streams_.end() ; ++i )
		if( (*i)->get_state() == need_input ) (*i)->put_result( r ) ;
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
	if( fn_.empty() ) {}
	else if( fn_ == "-" ) printout( cout, true ) ;
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

#if HAVE_ELK_SCHEME_H
Object StatStream::get_summary() const
{
	GC_Node ;

	Object n0 = Make_String( name_.data(), name_.size() ) ;		GC_Link( n0 ) ;
	Object n = Cons( Intern( "name" ), n0 ) ;					GC_Link( n ) ;
	Object t = Cons( Intern( "total" ), Make_Integer( total_ ) ) ; GC_Link( t ) ;
	Object m = Cons( Intern( "mapped" ), Make_Integer( mapped_ ) ) ; GC_Link( m ) ;
	Object mu = Cons( Intern( "mapuniq" ), Make_Integer( mapped_u_ ) ) ; GC_Link( mu ) ;
	Object d = Cons( Intern( "distinct" ), Make_Integer( different_ ) ) ; GC_Link( d ) ;
	Object gc = Cons( Intern( "gc-raw" ), Make_Flonum( 100*(float)bases_gc_ / (float)bases_ ) ) ; GC_Link( gc ) ;
	Object gcm = Cons( Intern( "gc-mapped" ), Make_Flonum( 100*(float)bases_gc_m_ / (float)bases_m_ ) ) ; GC_Link( gcm ) ;
	Object b = Cons( Intern( "bases" ), Make_Integer( bases_ ) ) ; GC_Link( b ) ;
	Object b2 = Cons( Intern( "bases-m2" ), Make_Integer( bases_squared_ ) ) ; GC_Link( b2 ) ;
	Object bm = Cons( Intern( "mapbases" ), Make_Integer( bases_m_ ) ) ; GC_Link( b ) ;
	Object bm2 = Cons( Intern( "mapbases-m2" ), Make_Integer( bases_m_squared_ ) ) ; GC_Link( b2 ) ;

	Object r = Cons( n, Cons( t, Cons( m, Cons( mu, Cons( d, Cons( gc, Cons( gcm, Cons( b, Cons( b2, Cons( bm, Cons( bm2, Null ))))))))))) ;
	GC_Unlink ;
	return r ;
}
#endif

void MismatchStats::put_result( const Result& r )
{
	for( int i = 0 ; i != r.hit_size() ; ++i )
	{
		const Hit& h = r.hit(i) ;
		assert( h.has_aln_ref() && h.has_aln_qry() ) ;
		for( size_t j = 0 ; j != h.aln_ref().size() && j != h.aln_qry().size() ; ++j )
		{
			char a = h.aln_ref()[j], b = h.aln_qry()[j] ;
			int u = a == 'A' ? 0 : a == 'C' ? 1 : a == 'G' ? 2 : a == 'T' ? 3 : -1 ;
			int v = b == 'A' ? 0 : b == 'C' ? 1 : b == 'G' ? 2 : b == 'T' ? 3 : -1 ;
			if( u >= 0 && v >= 0 ) ++mat_[u][v] ;
		}
	}
}

#if HAVE_ELK_SCHEME_H
Object MismatchStats::get_summary() const
{
	GC_Node ;
	Object r = Make_Vector( 4, False ) ;
	GC_Link( r ) ;
	for( int i = 0 ; i != 4 ; ++i )
	{
		Object s = Make_Vector( 4, False ) ;
		for( int j = 0 ; j != 4 ; ++j )
			VECTOR(s)->data[j] = Make_Integer( mat_[i][j] ) ;
		VECTOR(r)->data[i] = s ;
	}
	GC_Unlink ;
	return r ;
}
#endif

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

inline static bool good( char x ) { return x == 'A' || x == 'C' || x == 'G' || x == 'T' ; }

void DivergenceStream::put_result( const Result& r ) 
{
	const Hit *pri = hit_to( r, primary_genome_ ) ;
	const Hit *sec = hit_to( r, secondary_genome_ ) ;
	if( pri && sec )
	{
		if( !pri->has_aln_ref() || !pri->has_aln_qry() ||
				!sec->has_aln_ref() || !sec->has_aln_qry() )
			throw "Divergence: need reconstructed alignments" ;

		string::const_iterator pri_ref = pri->aln_ref().begin(),
					pri_qry = pri->aln_qry().begin(),
					pri_qry_end = pri->aln_qry().end(),
					sec_ref = sec->aln_ref().begin(),
					sec_qry = sec->aln_qry().begin(),
					sec_qry_end = sec->aln_qry().end() ;

		for( int skip = 0 ; skip != chop_ && pri_qry != pri_qry_end ; )
		{
			if( *pri_qry != '-' ) ++skip ;
			++pri_qry, ++pri_ref ;
			++sec_qry, ++sec_ref ;
		}

		for( int skip = 0 ; skip != chop_ && pri_qry != pri_qry_end ; )
		{
			--pri_qry_end ;
			--sec_qry_end ;
			if( *pri_qry_end != '-' ) ++skip ;
		}

		while( pri_qry != pri_qry_end ) {
			while( pri_qry != pri_qry_end && !good( *pri_qry ) )
				++pri_qry, ++pri_ref ;
			while( sec_qry != sec_qry_end && !good( *sec_qry ) )
				++sec_qry, ++sec_ref ;

			char p = *pri_ref, s = *sec_ref, q = *pri_qry ;
			if( q != *sec_qry ) throw "disagreement in alignments" ;

			if( good(p) && good(s) && good(q) ) {
				if( q == p && q == s ) ++b1 ;
				else if( q == p && q != s ) ++b2 ;
				else if( q != s && p == s ) ++b3 ;
				else if( q == s && q != p ) ++b4 ;
				else ++b5 ;
			}

			++pri_ref ;
			++sec_ref ;
			++pri_qry ;
			++sec_qry ;
		}
	}
}

#if HAVE_ELK_SCHEME_H
static inline Object Make_StringL( const char* c ) { return Make_String( c, strlen(c) ) ; }

Object DivergenceStream::get_summary() const
{
	GC_Node ;
	Object aa = Make_Vector( 5, False ) ;
	GC_Link( aa ) ;
	Object bb = Make_Vector( 5, False ) ;
	GC_Link( bb ) ;
	Object cd = False ;

	Object* b = VECTOR(bb)->data ;
	b[0] = Make_Flonum( b1 ) ;
	b[1] = Make_Flonum( b2 ) ;
	b[2] = Make_Flonum( b3 ) ;
	b[3] = Make_Flonum( b4 ) ;
	b[4] = Make_Flonum( b5 ) ;

	bb = Cons( Make_StringL( "raw-counts" ), bb ) ;
	bb = Cons( bb, Null ) ;
	bb = Cons( Cons( Make_StringL( "raw-div" ), Make_Flonum( 2.0*b4 / (b2+b4) )), bb ) ;
	
	if( int64_t d1 = 3*b1 - b2 + 3*b3 - b4 - b5 ) {
		double e = ( 3 * b3 - 3 * b4 ) / (double)d1 ;
		double d = 4*e - 3 ;

		double a2, a4 ;

		Object* a = VECTOR(aa)->data ;
		a[0] = Make_Flonum(      ( -3*b1 + b1*e + b3*e )/d ) ;
		a[1] = Make_Flonum( a2 = ( -3*b2 + b2*e + b4*e + b5*e )/d ) ;
		a[2] = Make_Flonum(      ( -3*b3 + 3*b1*e + 3*b3*e )/d ) ;
		a[3] = Make_Flonum( a4 = ( -3*b4 + b2*e + b4*e + b5*e )/d ) ;
		a[4] = Make_Flonum(      ( -3*b5 + 2*b2*e + 2*b4*e + 2*b5*e )/d ) ;

		cd = Make_Flonum( 2*a4 / (a2+a4) ) ;
	}
	else aa = False ;

	aa = Cons( Make_StringL( "corrected-counts" ), aa ) ;
	aa = Cons( aa, bb ) ;
	aa = Cons( Cons( Make_StringL( "corrected-div" ), cd ), aa ) ;

	GC_Unlink ;
	return aa ;
}
#endif

} ; // namespace

