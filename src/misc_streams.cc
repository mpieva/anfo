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

Header MergeStream::fetch_header()
{
	state_ = end_of_stream ;
	for( size_t i = 0 ; i != streams_.size() ; )
	{
		Header h = streams_[i]->fetch_header() ;
		merge_sensibly( hdr_, h ) ;

		if( h.is_sorted_by_name() ) {
			if( mode_ == unknown ) mode_ = by_name ;
			else if( mode_ != by_name ) 
				throw "MergeStream: inconsistent sorting of input (by_name)" ;
		}
		else if( h.is_sorted_by_all_genomes() ) {
			if( mode_ == unknown ) {
				mode_ = by_coordinate ;
				gs_.clear() ;
			}
			else if( mode_ != by_coordinate || !gs_.empty() )
				throw "MergeStream: inconsistent sorting of input (by_all_genomes)" ;
		}
		else if( h.is_sorted_by_coordinate_size() ) {
			if( mode_ == unknown ) {
				mode_ = by_coordinate ;
				gs_.assign( h.is_sorted_by_coordinate().begin(), h.is_sorted_by_coordinate().end() ) ;
			}
			else if( mode_ != by_coordinate || (int)gs_.size() != h.is_sorted_by_coordinate_size()
					|| !equal( gs_.begin(), gs_.end(), h.is_sorted_by_coordinate().begin() ) )
				throw "MergeStream: inconsistent sorting of input (by_genome)" ;
		}

		if( streams_[i]->get_state() == have_output )
		{
			rs_.push_back( streams_[i]->fetch_result() ) ;
			state_ = have_output ;
			++i ;
		}
		else
		{
			merge_sensibly( foot_, streams_[i]->fetch_footer() ) ;
			streams_.erase( streams_.begin() + i ) ;
		}
	}
	return hdr_ ;
}

//! \todo uses linear search, but a heap would be more appropriate
Result MergeStream::fetch_result() 
{
	if( state_ != have_output ) throw "logic error: MergeStream::fetch_result called in wrong stream state" ;
	if( rs_.size() != streams_.size() ) throw "logic error: # of records does not equal # of streams in MergeStream::fetch_result" ;
	if( !streams_.size() ) throw "logic error: no streams left in MergeStream::fetch_result" ;

	int min_idx = 0 ;
	for( size_t i = 1 ; i != rs_.size() ; ++i ) 
		if( ( mode_ == by_coordinate && by_genome_coordinate( gs_ )( &rs_[ i ], &rs_[ min_idx ] ) )
				|| ( mode_ == by_name && by_seqid()( &rs_[ i ], &rs_[ min_idx ] ) ) )
			min_idx = i ;

	string min_name = rs_[ min_idx ].read().seqid() ;
	Result res ;
	for( size_t i = 0 ; i != rs_.size() ; )
	{
		if( rs_[i].read().seqid() == min_name )
		{
			if( res.IsInitialized() ) merge_sensibly( res, rs_[i] ) ;
			else res = rs_[i] ;

			Stream &sm = *streams_[i] ;
			if( sm.get_state() == have_output ) rs_[i++] = sm.fetch_result() ; 
			else
			{
				merge_sensibly( foot_, sm.fetch_footer() ) ;
				streams_.erase( streams_.begin() + i ) ;
				rs_.erase( rs_.begin() + i ) ;
				if( rs_.empty() ) state_ = end_of_stream ;
			}
		}
        else ++i ;
	}
	return res ;
}

Header NearSortedJoin::fetch_header()
{
	state_ = end_of_stream ;
	for( size_t i = 0 ; i != streams_.size() ; )
	{
		merge_sensibly( hdr_, streams_[i]->fetch_header() ) ;
		++nstreams_ ;
		if( streams_[i]->get_state() == have_output )
		{
			state_ = have_output ;
			++i ;
		}
		else
		{
			merge_sensibly( foot_, streams_[i]->fetch_footer() ) ;
			streams_.erase( streams_.begin() + i ) ;
		}
	}
	cur_input_ = streams_.begin() ;
	return hdr_ ;
}

//! Enough information is read until one sequence has been reported on
//! by each input stream.  Everything else is buffered if necessary.
//! This works fine as long as the inputs come in in roughly the same
//! order and are complete.
Result NearSortedJoin::fetch_result()
{
	while( !streams_.empty() ) {
		if( cur_input_ == streams_.end() ) cur_input_ = streams_.begin() ;

		Stream &s = **cur_input_ ;
		if( s.get_state() != have_output ) throw "logic error: NearSortedJoin::fetch_result wasn't provided with input" ;

		Result r = s.fetch_result() ;
		if( s.get_state() == end_of_stream )
		{
			merge_sensibly( foot_, s.fetch_footer() ) ;
			cur_input_ = streams_.erase( cur_input_ ) ;
		}
		else ++cur_input_ ;

        if( buffer_.size() >= 100000 ) throw "input must be complete and nearly sorted the same way for join to work" ;

        if( ((nread_ + nwritten_) & 0xFFFF) == 0 )
        {
            stringstream s ;
            s << "NearSortedJoin: " << nread_ << "in "
                << nwritten_ << "out "
                << buffer_.size() << "buf" ;
            progress_( Console::info, s.str() ) ;
        }

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

void StatStream::put_result( const Result& r )
{
	unsigned bases = r.read().sequence().size() ;
	unsigned gc = count( r.read().sequence().begin(), r.read().sequence().end(), 'G' )
		        + count( r.read().sequence().begin(), r.read().sequence().end(), 'C' ) ;
	unsigned count = std::max( r.members_size() + r.nmembers(), 1U ) ;
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
		if( !h->has_map_quality() || h->map_quality() >= 30 )
		{
			mapped_u_ += count ;
			++different_ ;
		}
	}
}

static inline Object Cons2U( uint64_t x, Object l )
{ return Cons( Make_Unsigned( x & ~(~0LL << 31) ), Cons( Make_Unsigned( x >> 31 ), l ) ) ; }

static inline Object ConsU( unsigned x, Object l )
{ return Cons( Make_Unsigned( x ), l ) ; }

Object StatStream::get_summary() const
{
	Object r = Null ;
	r = Cons2U( bases_m_squared_, r ) ;
	r = Cons2U( bases_m_, r ) ;
	r = Cons2U( bases_squared_, r ) ;
	r = Cons2U( bases_, r ) ;
	r = Cons2U( bases_gc_m_, r ) ;
	r = Cons2U( bases_gc_, r ) ;
	r = ConsU( different_, r ) ;
	r = ConsU( mapped_u_, r ) ;
	r = ConsU( mapped_, r ) ;
	r = ConsU( total_, r ) ;
	return r ;
}

void MismatchStats::put_result( const Result& r )
{
	for( int i = 0 ; i != r.hit_size() ; ++i )
	{
		const Hit& h = r.hit(i) ;
		if( !h.has_aln_ref() || !h.has_aln_qry() ) 
			throw "MismatchStats needs reconstructed alignments" ;

		for( size_t j = 0 ; j != h.aln_ref().size() && j != h.aln_qry().size() ; ++j )
		{
			char a = h.aln_ref()[j], b = h.aln_qry()[j] ;
			int u = a == 'A' ? 0 : a == 'C' ? 1 : a == 'G' ? 2 : a == 'T' ? 3 : -1 ;
			int v = b == 'A' ? 0 : b == 'C' ? 1 : b == 'G' ? 2 : b == 'T' ? 3 : -1 ;
			if( u >= 0 && v >= 0 ) ++mat_[u][v] ;
		}
	}
}

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

void DivergenceStream::put_header( const Header& h )
{
    Stream::put_header( h ) ;
    ancient_ = h.config().has_aligner() && h.config().aligner().has_mean_overhang_length() ;
}

// Excludes gaps, Ns, and nonsensical symbols
inline static bool good( char x ) { return x == 'A' || x == 'C' || x == 'G' || x == 'T' ; }

// Excludes transitions that may have been confounded by deamination.
// Idea is the same as in contamination checker: "If you see a C you can
// trust that you see a C."  (Excluding all transitions doesn't change
// the results.)
inline static bool fine( char x, char y, char z )
{
    if( x == 'C' && y == 'T' && z == 'T' ) return false ;
    if( x == 'T' && y == 'T' && z == 'C' ) return false ;
    if( x == 'G' && y == 'A' && z == 'A' ) return false ;
    if( x == 'A' && y == 'A' && z == 'G' ) return false ;
    return true ;
}

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
		}
		for( int skip = 0 ; skip != chop_ && pri_qry != pri_qry_end ; )
		{
			--pri_qry_end ;
			if( *pri_qry_end != '-' ) ++skip ;
		}
		for( int skip = 0 ; skip != chop_ && sec_qry != sec_qry_end ; )
		{
			if( *sec_qry != '-' ) ++skip ;
			++sec_qry, ++sec_ref ;
		}
		for( int skip = 0 ; skip != chop_ && sec_qry != sec_qry_end ; )
		{
			--sec_qry_end ;
			if( *sec_qry_end != '-' ) ++skip ;
		}

		for(;;)
		{
			while( pri_qry != pri_qry_end && *pri_qry == '-' ) ++pri_qry, ++pri_ref ;
			while( sec_qry != sec_qry_end && *sec_qry == '-' ) ++sec_qry, ++sec_ref ;
			if( pri_qry == pri_qry_end ) break ;
			if( sec_qry == sec_qry_end ) break ;

			char p = *pri_ref, s = *sec_ref, q = *pri_qry ;
			if( q != *sec_qry ) throw "disagreement in alignments (or ugly bug)" ;

			if( good(p) && good(s) && good(q) && (!ancient_ || fine(p,q,s)) ) {
				if( q == p && q == s ) ++b1 ;
				else if( q == p && q != s ) ++b2 ;
				else if( q != s && p == s ) ++b3 ;
				else if( q == s && q != p ) ++b4 ;
				else ++b5 ;
			}

			++pri_ref ;
			++pri_qry ;
			++sec_ref ;
			++sec_qry ;
		}
	}
}

Object DivergenceStream::get_summary() const
{
	Object bb = Make_Vector( 5, False ) ;
	Object* b = VECTOR(bb)->data ;
	b[0] = Make_Flonum( b1 ) ;
	b[1] = Make_Flonum( b2 ) ;
	b[2] = Make_Flonum( b3 ) ;
	b[3] = Make_Flonum( b4 ) ;
	b[4] = Make_Flonum( b5 ) ;
	return bb ;
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
	hdr_.clear_command_line() ;
	hdr_.clear_sge_job_id() ;
	hdr_.clear_sge_task_id() ;
	hdr_.DiscardUnknownFields() ;
}

bool Sanitizer::xform( Result& r )
{
	r.clear_aln_stats() ;
	return true ;
}

AgreesWithChain::AgreesWithChain( const string& l, const string& r, const pair<istream*,string>& p ) 
	: left_genome_( l ), right_genome_( r ), map_()
{
}


AgreesWithChain::Chains::iterator AgreesWithChain::find_any_most_specific_overlap( 
		unsigned start, unsigned end, Chains *chains )
{
	Chains::iterator best = chains->end() ;
	for(;;)
	{
		// first chain that starts at or to the right of end (this one
		// doesn't overlap, but the one _before_ it might)
		Chains::iterator i = chains->lower_bound( end ) ;
		if( i == chains->begin() ) return best ;

		Entry* e = &(--i)->second ;
		// no overlap? return parent.
		if( start >= e->left_end || i->first >= end ) return best ;

		// continue with nested chains, making current entry the best
		best = i ;
		chains = &e->nested ;
	}
}

// Idea is to find the most specific chain that overlaps a read, then
// make sure the read is contained in it.
bool AgreesWithChain::xform( Result& r ) 
{
	const Hit* lh = hit_to( r, left_genome_ ) ;
	const Hit* rh = hit_to( r, right_genome_ ) ;
	if( !lh || !rh ) return false ;

	// find chain hierarchy for correct chromosome
	Map1::iterator i1 = map_.find( lh->sequence() ) ;
	if( i1 == map_.end() ) {
        cerr << "chromosome not found" << endl ;
        return false ;
    }

	// find most specific chain that overlaps left hit
	Chains::iterator i2 = find_any_most_specific_overlap( 
			lh->start_pos(),
			lh->start_pos() + abs(lh->aln_length()),
			&i1->second ) ;
	if( i2 == i1->second.end() ) return false ;
	
	// check if found chain contains both hits
	if(
			lh->start_pos() < i2->first ||
			lh->start_pos() + abs( lh->aln_length() ) > i2->second.left_end ||
			rh->start_pos() < i2->second.right_start ||
			rh->start_pos() + abs( rh->aln_length() ) > i2->second.right_end ||
			rh->sequence() != i2->second.right_chr 
	  ) return false ;

	// check strand
	if( ((rh->aln_length() < 0) != (lh->aln_length() < 0)) != i2->second.strand ) return false ;

	// Now both hits are covered by a single, most specific chain.
	// (Note: currently, hitting anywhere in the chain is fine; this
	// should probably be made somewhat more precise)
	return true ;
}


} ; // namespace

