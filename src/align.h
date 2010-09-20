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
class ExtendAlignment {
	// private:
	public:
		struct Cell {
			Logdom score ;
			uint8_t from_state ;
			uint8_t from_x_offset ;
			uint8_t from_y_offset ;

			bool operator < ( const Cell& rhs ) const { return score < rhs.score ; }
		} ;

		int width_ ;
		std::vector< Array< adna_parblock::num_states, Cell > > cells_ ;
		std::vector< int > mins_, maxs_ ;

		Logdom limit_, result_ ;
		int max_s_, max_x_, max_y_, max_tail_ ;

		void put( int s, int os, int x, int xo, int y, int y0, Logdom z ) ;
		void extend(const adna_parblock &pb_,Logdom score,int s,int x,int y,DnaP ref,const QSequence::Base *qry );

	public:
		ExtendAlignment() {}
		ExtendAlignment( const adna_parblock& pb, DnaP reference, const QSequence::Base *query, Logdom limit ) ;

		Logdom get_result() const { return result_ ; }
		int max_x() const { return max_x_ ; }

		void backtrace( std::vector<unsigned>& ) const ;

		void swap( ExtendAlignment& rhs ) 
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

struct ExtendBothEnds {
	ExtendAlignment forwards_, backwards_ ;
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

	void swap( ExtendBothEnds &rhs )
	{
		forwards_.swap( rhs.forwards_ ) ;
		backwards_.swap( rhs.backwards_ ) ;
		score_.swap( rhs.score_ ) ;
	}
} ;

#endif
