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

#ifndef INCLUDED_ALIGN_H
#define INCLUDED_ALIGN_H

#include "align_fwd.h"
#include "index.h"
#include "stream.h"

#include <cmath>
#include <deque>
#include <ostream>
#include <sstream>

template< int N, typename T > class Array
{
	private:
		T arr_[N] ;

	public:
		Array() {}
		Array( T t ) { for( int i = 0 ; i != N ; ++i ) arr_[i] = t ; }

		T& operator[] ( int i ) { return arr_[i] ; }
		const T& operator[] ( int i ) const { return arr_[i] ; }
} ;

//! \brief an alignment that has been seeded
//! When initializing an alignment, the seed region is traversed and
//! then extended greedily.  We store the final seed region (a pointer
//! into the reference, an offset into the query, the seeded length) and
//! the score over the seed region.
struct SeededAlignment {
	DnaP reference_ ;
	Logdom score_ ;
	int qoffs_ ;
	int size_ ;

	SeededAlignment() : qoffs_(0), size_(0) {}
	SeededAlignment( const adna_parblock& pb, DnaP reference, const QSequence& query, int qoffs, int size )
		: reference_( reference ), score_( Logdom::one() ), qoffs_( qoffs ), size_( size )
	{
		// greedy initialization: run over the seed, accumulating a
		// score.  then extend greedily as long as there are
		// matches, store the resulting score.
		for( int i = 0 ; i <= size_ ; ++i )
			score_ *= pb.subst_penalty( 0, reference_[i], query[qoffs_+i] ) ;

		for( ; reference_[size_] && query[qoffs_+size_].ambicode &&
				reference_[size_] == query[qoffs_+size_].ambicode ; ++size_ )
			score_ *= pb.subst_penalty( 0, reference_[size_], query[qoffs_+size_] ) ;

		for( ; reference_[-1] && query[qoffs_-1].ambicode &&
				reference_[-1] == query[qoffs_-1].ambicode ;
				++size_, --reference_, --qoffs_ )
			score_ *= pb.subst_penalty( 0, reference_[-1], query[qoffs_-1] ) ;
	}
} ;


struct FullCell {
	Logdom score ;
	uint8_t from_state ;
	uint8_t from_x_offset ;
	uint8_t from_y_offset ;

	void clear() { score = Logdom::null() ; }
	void assign( Logdom z, int os, int xo, int yo )
	{ 
		if( score < z )
		{
			score = z ;
			from_state = os ;
			from_x_offset = xo ;
			from_y_offset = yo ;
		}
	} 
} ;

struct SimpleCell {
	Logdom score ;

	void clear() { score = Logdom::null() ; }
	void assign( Logdom z, int, int, int ) { if( score < z ) score = z ; } 
} ;

//! \brief extension of an alignment through DP
//! This starts with two pointers and a cost limit, and it returns the
//! penalty that was incurred.  The limit can be exceeded while
//! producing an alignment (it would be foolish to throw it away), if
//! nothing is found, std::numeric_limits<uint32_t>::max() is returned.
//!
//! We align until we hit a zero in either query or reference.  Starting
//! state is 0, late the state is inly held implicitly.
//!
//! If only there were real lexical closures... *sigh*
template< typename Cell > class ExtendAlignment {
	private:
		int width_ ;
		std::vector< Array< adna_parblock::num_states, Cell > > cells_ ;
		std::vector< int > mins_, maxs_ ;

		Logdom limit_, result_ ;
		int max_s_, max_x_, max_y_, max_tail_ ;

		void extend(const adna_parblock &pb_,Logdom score,int s,int x,int y,DnaP ref,const QSequence::Base *qry );

	public:
		ExtendAlignment() {}
		ExtendAlignment( const adna_parblock& pb, DnaP reference, const QSequence::Base *query, Logdom limit ) ;

		Logdom get_result() const { return result_ ; }
		int max_x() const { return max_x_ ; }

		void backtrace( std::vector<unsigned>& ) const ;

		void swap( ExtendAlignment<Cell>& rhs ) throw()
		{
			std::swap( width_, rhs.width_ ) ;
			std::swap( cells_, rhs.cells_ ) ;
			std::swap( mins_, rhs.mins_ ) ;
			std::swap( maxs_, rhs.maxs_ ) ;
			std::swap( limit_, rhs.limit_ ) ;
			std::swap( result_, rhs.result_ ) ;
			std::swap( max_s_, rhs.max_s_ ) ;
			std::swap( max_x_, rhs.max_x_ ) ;
			std::swap( max_y_, rhs.max_y_ ) ;
			std::swap( max_tail_, rhs.max_tail_ ) ;
		}
} ;

template< typename Cell > struct ExtendBothEnds {
	ExtendAlignment<Cell> forwards_, backwards_ ;
	Logdom score_ ;

	ExtendBothEnds() {}
	ExtendBothEnds( const adna_parblock& pb, const QSequence& query, const SeededAlignment& seed, Logdom limit ) ;
	
	//! \brief backtraces an alignment and return a CIGAR line
	//!
	//! Backtracing works by simply walking the chain of ols states and
	//! x/y offsets stored in the DP matrix.  See output.proto for the
	//! encoding of the produced CIGAR lines.
	//!
	//! \param minpos will be filled by position of first aligned
	//!               reference base
	//! \param maxpos will be filled by position of first non-aligned
	//!               reference base, so that maxpos-minpos gives the
	//!               aligned length
	//! \return binary CIGAR string
	//! \internal
	//! \todo Write mismatches with different code (to allow various
	//!       calculations without the genome being available).
	std::vector<unsigned> backtrace( const SeededAlignment& seed, DnaP &minpos, DnaP &maxpos ) const ;

	void swap( ExtendBothEnds<Cell>& rhs ) throw()
	{
		forwards_.swap( rhs.forwards_ ) ;
		backwards_.swap( rhs.backwards_ ) ;
		score_.swap( rhs.score_ ) ;
	}
} ;


namespace {
	inline int query_length( const QSequence::Base *query )
	{
		int r = 0 ;
		while( query->ambicode ) ++query, ++r ;
		return r ;
	}
} ;

// Alignment proper: the intial greedy matching must have been done,
// here we extend one side of this into a full alignment, as long as it
// doesn't score more than a prescribed limit.
template< typename Cell >
ExtendAlignment<Cell>::ExtendAlignment( const adna_parblock& pb, DnaP reference, const QSequence::Base *query, Logdom limit ) :
	width_( 2*query_length( query )+2 ), cells_( width_*width_ ), mins_( width_ ), maxs_( width_ ), 
	limit_( limit ), result_( Logdom::null() )
{
	if( limit > Logdom::one() ) return ;

	mins_[0] = 0 ;
	maxs_[0] = 1 ;
	cells_[0][0].assign( Logdom::one(), 0, 0, 0 ) ;
	for( int s = 1 ; s != adna_parblock::num_states ; ++s ) cells_[ 0 ][ s ].clear() ; 

	for( int y = 0 ; y != width_-1 && mins_[y] != maxs_[y] ; ++y )
	{
		assert( y <= width_ ) ;
		assert( mins_[y] >= 0 ) ;
		assert( maxs_[y] <= width_ ) ;
		assert( mins_[y] <= maxs_[y] ) ;

		mins_[y+1] = maxs_[y+1] = 0 ;

		// expand the current row for each state in turn... of course,
		// each state is a special case.
		for( int x = mins_[y] ; x != width_-1 && x != maxs_[y] ; ++x )
		{
			for( int s = 0 ; s != adna_parblock::num_states ; ++s )
			{
				assert( width_*y + x < width_ * width_ ) ;
				assert( width_*y + x < (int)cells_.size() ) ;

				Logdom score = cells_[ width_*y + x ][ s ].score ;
				if( score > limit_ ) extend( pb, score, s, x, y, reference+x, query+y ) ;
			}
		}
	}
}

#define PUT(ns,xo,yo,z) 																\
	if( (z) >= limit_ ) { 																\
		if( mins_[ y+(yo) ] == maxs_[ y+(yo) ] ) 										\
			mins_[ y+(yo) ] = maxs_[ y+(yo) ] = x+(xo) ; 								\
																						\
		for( ; maxs_[ y+(yo) ] <= x+(xo) ; ++maxs_[ y+(yo) ] ) 							\
			for( int ss = 0 ; ss != adna_parblock::num_states ; ++ss ) 					\
				cells_[ width_*(y+(yo)) + maxs_[ y+(yo) ] ][ ss ].clear() ; 			\
																						\
		assert( y+(yo) < width_ ) ; 													\
		assert( x+(xo) < width_ ) ; 													\
		assert( width_*(y+(yo)) + x+(xo) < width_ * width_ ) ; 							\
																						\
		cells_[ width_*(y+(yo)) + x+(xo) ][ ns ].assign( z, s, xo, yo ) ; 				\
	} else {}


// what to do?  
// If in matching state, we know there's no immediate match, so we can...
// - mismatch
// - detect deamination and change to SS state (while matching)
// - open ref gap
// - open query gap
//
// If a gap is open, we can...
// - extend it
// - close it
//
// If we hit a gap symbol, we must...
// - start over at second half in initial state

template< typename Cell >
void ExtendAlignment<Cell>::extend( const adna_parblock &pb_, Logdom score, int s, int x, int y, DnaP ref, const QSequence::Base *qry )
{
	// Note the penalties: The appropriate substitution penalty is
	// applied whenever we (mis-)match two codes, the gap open penalties
	// are applied when opening/extending a gap, the
	// overhang_enter_penalty is applied when changing to SS mode and
	// the overhang_ext_penalty is applied whenever moving along the
	// query while single stranded, even when a gap is open!  This gives
	// correct scores for a geometric distribution of overhang lengths.
	if( !*ref || !qry->ambicode ) 
	{
		// We hit a gap in either the reference or the query.  Whatever
		// is left of the query (if any) must be penalized.  To do this,
		// we virtually extend the reference with Ns and align to those.
		// This is a white lie in that it will overestimate the real
		// penalty, but that's okay, because such an alignment isn't all
		// that interesting in reality anyway.  Afterwards we're
		// finished and adjust result and limit accordingly.
		int yy = 0 ;
		for( ; qry[ yy ].ambicode ; ++yy )
		{
			score *= pb_.subst_penalty( s, 15, qry[ yy ] ) ;
			if( s & adna_parblock::mask_ss ) score *= pb_.overhang_ext_penalty ;
		}
		if( score > result_ ) {
			result_ = limit_ = score ;
			max_s_ = s ;
			max_x_ = x ;
			max_y_ = y ;
			max_tail_ = yy ;
		}
	}
	else if( (s & adna_parblock::mask_gaps) == 0 )
	{
		// no gaps open --> mismatch, open either gap, enter SS
		PUT( s, 1, 1, score * pb_.subst_penalty( s, *ref, *qry ) 
				* ( s & adna_parblock::mask_ss ? pb_.overhang_ext_penalty : Logdom::one() ) ) ;
		if( *ref != qry->ambicode ) {	// only on a mismatch try anything fancy
			PUT( s | adna_parblock::mask_gap_qry, 1, 0, score * pb_.gap_open_penalty ) ;
			PUT( s | adna_parblock::mask_gap_ref, 0, 1, score * pb_.gap_open_penalty ) ;
			if( pb_.overhang_enter_penalty.is_finite() && (s & adna_parblock::mask_ss) == 0 ) {
				// To enter single stranded we require that the penalty for
				// doing so is immediately recovered by the better match.
				// This is easily the case for the observed deamination
				// rates in aDNA.
				Logdom p0 = pb_.subst_penalty( s, *ref, *qry ) ;
				Logdom p4 = pb_.subst_penalty( s | adna_parblock::mask_ss, *ref, *qry ) 
					* pb_.overhang_enter_penalty * pb_.overhang_ext_penalty ;
				if( p4 > p0 ) { PUT( s | adna_parblock::mask_ss, 1, 1, score * p4 ) ; }
			}
		}
	}
	else if( (s & adna_parblock::mask_gaps) == adna_parblock::mask_gap_ref )
	{
		PUT( s, 0, 1, score * pb_.gap_ext_penalty *
				( s & adna_parblock::mask_ss ? pb_.overhang_ext_penalty : Logdom::one() ) ) ;
		PUT( s & ~adna_parblock::mask_gap_ref, 1, 1, score * pb_.subst_penalty( s, *ref, *qry ) *
				( s & adna_parblock::mask_ss ? pb_.overhang_ext_penalty : Logdom::one() ) ) ;
	}
	else
	{
		PUT( s, 1, 0, score * pb_.gap_ext_penalty ) ;
		PUT( s & ~adna_parblock::mask_gap_qry, 1, 1, score * pb_.subst_penalty( s, *ref, *qry ) *
				( s & adna_parblock::mask_ss ? pb_.overhang_ext_penalty : Logdom::one() ) ) ;
	}
}

#undef PUT

// Extension of both sides.  We first run a forward extension at half
// the limit.  If this succeeds, we do the backwards extension limited
// to whatever is left.  If forward extension fails, we do backwards
// extension to half the limit, then add forward using up what's left.
//
// Only one alignment is produced, but we make sure it is the cheapest
// one.  The score may exceed the limit, if we happen to finish right
// when stepping over the limit.  If really nothing is found, we return
// Logdom::null().

template< typename Cell>
ExtendBothEnds<Cell>::ExtendBothEnds(
		const adna_parblock& pb,
		const QSequence& query,
		const SeededAlignment& seed,
		Logdom limit ) :
	forwards_(
			pb,
			seed.reference_ + seed.size_,
			query.start() + seed.qoffs_ + seed.size_,
			( limit / seed.score_ ).sqrt() ),
	backwards_(
			pb,
			seed.reference_.reverse() + 1,
			query.start() - seed.qoffs_ + 1,
			limit / seed.score_ / forwards_.get_result() ),
	score_( seed.score_ * forwards_.get_result() * backwards_.get_result() )
{
	if( score_.is_finite() ) return ;

	ExtendAlignment<Cell> backwards2(
			pb,
			seed.reference_.reverse() + 1,
			query.start() - seed.qoffs_ + 1,
			(limit / seed.score_).sqrt() ) ;
	ExtendAlignment<Cell> forwards2(
			pb,
			seed.reference_ + seed.size_,
			query.start() + seed.qoffs_ + seed.size_,
			limit / backwards2.get_result() / seed.score_ ) ;

	Logdom score2 = forwards2.get_result() * backwards2.get_result() * seed.score_ ;
	if( score2.is_finite() )
	{
		score_ = score2 ;
		forwards_.swap( forwards2 ) ;
		backwards_.swap( backwards2 ) ;
	}
}

template<> inline void ExtendAlignment<FullCell>::backtrace( std::vector<unsigned>& out ) const
{
	if( max_tail_ ) streams::push_i( out, max_tail_ ) ;
	for( size_t x = max_x_, y = max_y_, s = max_s_ ; x || y ; )
	{
		const FullCell& c = cells_[ width_*y + x ][ s ] ;
		if( !c.from_x_offset && !c.from_y_offset ) throw "stuck in backtracing" ;
		else if( !c.from_x_offset ) streams::push_i( out, c.from_y_offset ) ;
		else if( !c.from_y_offset ) streams::push_d( out, c.from_x_offset ) ;
		else if( c.from_y_offset == c.from_x_offset ) streams::push_m( out, c.from_x_offset ) ;
		else throw "inconsistency in backtracing" ;

		x -= c.from_x_offset ;
		y -= c.from_y_offset ;
		s = c.from_state ;
	}
}

template<> inline 
std::vector<unsigned> ExtendBothEnds<FullCell>::backtrace( const SeededAlignment& seed, DnaP &minpos, DnaP &maxpos ) const 
{
	minpos = seed.reference_ - backwards_.max_x() ;
	maxpos = seed.reference_ + seed.size_ + forwards_.max_x() ;

	std::vector<unsigned> trace ;
	backwards_.backtrace( trace ) ;

	trace.push_back( 0 ) ;
	streams::push_m( trace, seed.size_ ) ;
	trace.push_back( 0 ) ;

	std::vector<unsigned> rtrace ;
	forwards_.backtrace( rtrace ) ;
	std::copy( rtrace.rbegin(), rtrace.rend(), back_inserter( trace ) ) ;
	return trace ;
}

#endif
