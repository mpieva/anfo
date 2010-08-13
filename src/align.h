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

#include "index.h"
#include "logdom.h"
#include "stream.h"

#include <cmath>
#include <deque>
#include <ostream>
#include <sstream>

//! \page adna_alignment Operations on Ancient DNA Alignments
//! This type of alignment models aDNA as two additional states, one
//! with higher deamination rate, with an asymmetric substitution
//! matrix, and with affine gap costs.  Additional parameters are all
//! static.
//!
//! \todo We actually don't need the parameter set to be static, and it
//!       may even pay to have it configurable so differently
//!       parameterized alignments an be mixed in a single run.
//!
//! @{

//! \brief substitution matrix
//! We prepare a full matrix of 16 ambiguity codes vs. 16 ambiguity
//! codes, this avoid the need for expensive additions in the log-domain
//! should the need to align ambiguity codes arise.  First index is
//! "from" (reference code), second index is "to" (query code).
typedef Logdom subst_mat[16][16] ;

namespace config { class Aligner ; } ;

//! \brief parameters for simple_adna
//! One such parblock is a static variable for the aligner proper,
//! various support tools shunt additional structures around.
struct adna_parblock
{
	enum {
		mask_ss      = 1,
		mask_gap_ref = 2,
		mask_gap_qry = 4,

		mask_gaps = mask_gap_ref | mask_gap_ref,
		num_states = 6
	} ;

	adna_parblock() {} 
	adna_parblock( const config::Aligner& conf ) ;

	//! \brief DS substitution matrix, forward direction
	subst_mat ds_mat ;

	//! \brief SS substitution matrix, forward direction
	//! Deamination shows up as C->T as it is best understood this way.
	//! To process reverse-complemented deamination, we have to do
	//! rev-complemented lookups while actually moving in the forward
	//! (5'->3') direction.  \see simple_adna::subst_penalty()
	subst_mat ss_mat ;

	//! \brief Penalty for extending an overhang.
	//! Having a constant penalty for the overhang length models its
	//! length distribution as geometric.
	Logdom overhang_ext_penalty ;

	//! \brief penalty for entering SS state
	//! This is essentially the probability of having an overhang at
	//! all.
	Logdom overhang_enter_penalty ;

	//! \brief gap open penalty
	Logdom gap_open_penalty ;

	//! \brief gap extension penalty
	Logdom gap_ext_penalty ;


	Logdom subst_penalty( int s, Ambicode r, const QSequence::Base &qry ) const
	{
		// if reference has a gap, pretend it was an N
		r = !r ? 15 : r ;

		Logdom prob = Logdom::null() ;
		for( uint8_t p = 0 ; p != 4 ; ++p )
			prob += ( s & mask_ss ? ss_mat[r][1<<p] : ds_mat[r][1<<p] )
				* Logdom::from_phred( qry.qscores[p] ) ;
		return prob ;
	}
} ;

std::ostream& operator << ( std::ostream&, const adna_parblock& ) ;

#if 0
//! \brief backtraces an alignment and return a CIGAR line
//!
//! This backtraces an alignment after it has been created by
//! find_cheapest() called with a backtracing structure.  Backtracing
//! works by looking at two intermediate states, the \em current one and
//! its \em predecessor.  In between those two states, exactly one
//! invocation of forward() and one of greedy() have happened.  To
//! backtrace, we first check in which direction we moved (depends on
//! which half of the alignment we were in).  Next, if both the reference
//! and the query offsets differ between states, we copy symbols from
//! both (this corresponds to the greedy extension or a mismatch).  If
//! only one differs, we copy one symbol and fill it up with a gap.
//! Depending on the direction we moved in, the pair is added at the
//! front or the end of the trace.  If the internal state changed, we
//! don't trace at all (we just jump).  See output.proto for the
//! encoding of these CIGAR lines.
//!
//! \param cl the ClosedMap that was used in find_cheapest()
//! \param a final state to start backtracing from
//! \param minpos will be filled by position of smaller end of alignment
//! \param maxpos will be filled by position of greater end of alignment
//! \return binary CIGAR string
//! \internal
//! \todo Write mismatches with different code (to allow various
//!       calculations without the genome being available).

template< typename State > std::vector<unsigned>
backtrace( const typename State::ClosedMap &cl, const State *a, DnaP &minpos, DnaP &maxpos )
{
	std::vector<unsigned> fwd, rev ;

	// When the alignment finished, it pointed to the minimum
	// coordinate (to a gap actually).  Just store it.
	minpos = a->reference + a->ref_offs ;

	// Only trace back second state here, this ends up at the front of
	// the alignment, but we add stuff to the back as we generate it in
	// the wrong order.
	while( const State *b = *lookup( cl, *a ) ) 
	{
		if( (b->state & simple_adna::mask_dir) == 0 ) break ;
		int dr = b->ref_offs - a->ref_offs ;
		int dq = b->query_offs - a-> query_offs ;
		if( int m = std::min( dr, dq ) ) {
			dr -= m ;
			dq -= m ;
			streams::push_m( fwd, m ) ;
		}
		streams::push_i( fwd, dq ) ;
		streams::push_d( fwd, dr ) ;
		a = b ;
	}

	if( a->ref_offs != a->query_offs ) throw  "logic error: first half of backtrace has drifted" ;

	streams::push_m( fwd, -a->ref_offs-1 ) ;
	fwd.push_back( 0 ) ;

	// Skip one state, this is the one aligning the terminal gap.
	a = *lookup( cl, *a ) ; 
		
	// Now at the end of the first phase, we got a pointer to the
	// maximum coordinate (again a gap).  Store it.
	maxpos = a->reference + a->ref_offs ;

	// Trace back first state now, generating a new trace which needs to
	// be reversed in the end.
	while( const State *b = *lookup( cl, *a ) ) 
	{
		int dr = a->ref_offs - b->ref_offs ;
		int dq = a->query_offs - b-> query_offs ;
		if( int m = std::min( dr, dq ) ) {
			dr -= m ;
			dq -= m ;
			streams::push_m( rev, m ) ; 
		}
		streams::push_i( rev, dq ) ;
		streams::push_d( rev, dr ) ;
		a = b ;
	}

	if( a->ref_offs != a->query_offs ) throw "logic error: second half of backtrace has drifted" ;

	// Now (*b) is null, (*a) is the last state ever generated.  Now
	// trace further until we hit the initial state (both offsets
	// vanish).  We can again add this to t2, as we're going in the same
	// direction.
	streams::push_m( rev, a->ref_offs ) ;

	fwd.insert( fwd.end(), rev.rbegin(), rev.rend() ) ;
	return fwd ;
}
#endif

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

	SeededAlignment() {}
	SeededAlignment( const adna_parblock& pb, DnaP reference, const QSequence& query, int qoffs_, int size ) ;
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
//!
//! \todo This forward/backward business is brittle.  We should just
//!       flip the pointers: easy for the genome, for the query we
//!       should simply duplicate the sequence.
class ExtendAlignment {
	// private:
	public:
		struct Cell {
			Logdom score ;
			uint8_t from_state ;
			uint8_t from_x_offset ;
			uint8_t from_y_offset ;

			Cell() : score( Logdom::null() ) {}
			bool operator < ( const Cell& rhs ) const { return score < rhs.score ; }
		} ;

		struct Line {
			std::vector< Array< adna_parblock::num_states, Cell > > cells ;
			int min ;
		} ;

		Logdom limit_, result_ ;
		std::deque< Line > matrix ;

		void put( Line& l, int s, int x, int xo, Logdom z, int yo, int os ) ;
		void extend(const adna_parblock &pb_,Logdom score,int s,int x,int y,DnaP ref,const QSequence::Base *qry );

	public:
		ExtendAlignment() {}
		ExtendAlignment( const adna_parblock& pb, DnaP reference, const QSequence::Base *query, Logdom limit ) ;

		Logdom get_result() const { return result_ ; }

		std::vector<unsigned> backtrace( DnaP &minpos, DnaP &maxpos ) 
		{
			// minpos = reference_ ;
			// maxpos = reference_ + seedsize_ ;

			return std::vector<unsigned>() ;
		}

		void swap( ExtendAlignment& rhs ) 
		{
			std::swap( limit_, rhs.limit_ ) ;
			std::swap( result_, rhs.result_ ) ;
			std::swap( matrix, rhs.matrix ) ;
		}
} ;

struct ExtendBothEnds {
	ExtendAlignment forwards_, backwards_ ;
	Logdom score_ ;

	ExtendBothEnds() {}
	ExtendBothEnds( const adna_parblock& pb, const QSequence& query, const SeededAlignment& seed, Logdom limit ) ;
	
	std::vector<unsigned> backtrace( const SeededAlignment& seed, DnaP &minpos, DnaP &maxpos ) 
	{
		minpos = seed.reference_ ;
		maxpos = seed.reference_ + seed.size_ ;

		return std::vector<unsigned>() ;
	}

	void swap( ExtendBothEnds &rhs )
	{
		forwards_.swap( rhs.forwards_ ) ;
		backwards_.swap( rhs.backwards_ ) ;
		score_.swap( rhs.score_ ) ;
	}
} ;

#endif
