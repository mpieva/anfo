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

#ifndef INCLUDED_TRIM_H
#define INCLUDED_TRIM_H

#include <algorithm>
#include <vector>

#include "sequence.h"

namespace {
	template< typename T > T max3( T a, T b, T c )
	{ return std::max( std::max( a, b ), c ) ; }
}

inline bool match( char     a, char     b ) { return to_ambicode(a) & to_ambicode(b) ; }
inline bool match( Ambicode a, char     b ) { return             a  & to_ambicode(b) ; }
inline bool match( char     a, Ambicode b ) { return to_ambicode(a) &             b  ; }
inline bool match( Ambicode a, Ambicode b ) { return             a  &             b  ; }

inline bool match( const QSequence::Base &a, char b ) { return a.ambicode & to_ambicode(b) ; }
inline bool match( char a, const QSequence::Base &b ) { return to_ambicode(a) & b.ambicode ; }

/*! \brief aligns two sequences to find an adapter sequence.
 
	The alignment algorithm following Eugene W. Myers: "An O(ND)
	Difference Algorithm and Its Variations", adapted to 'overlap
	alignment'.

	We want to align two sequences so they overlap with maximum score.
	This is not possible with the plain "minimum distance" method, since
	no overlap always results in the minimum distance of zero.  We do
	need to score the overlapping part, but not the same as an ordinary
	difference.  Now suppose you have an alignment:

		xxxxxx--
		--yyyyzz

    It scores (k * mat/2 - d * (mat-mis)), where k is |x|+|y| (see the
    megablast paper by Zhang for an explanation; we stop keeping score
    once we reach the end of x).  Suppose x is longer:

		xxxxxxx--
		---yyyyzzz

	This scores ((k+1) * mat/2 - (d+D) * (mat-mis)), and we want to set
	D such that this is the same as before.  This gives:

		D = 1/2 * mat / (mat-mis)

	and for typical parameters mat=1 and mis=-3 it equals D=1/8.
	Therefore, end gaps should score 1/8 (but not for local alignments
	(not implemented here), where only an X-drop algorithm would be
	appropriate).

	To implement this, the following changes to the original algorithm
	are necessary:
	- allow a start gap that grows 8x quicker than a normal gap,
	  consider the many more involved diagonals,
    - once finished, prefer alignments on lower numbered diagonals.


	This function will try to overlap the left end of sequence A (query
	sequence) with the right end of sequence B (adapter sequence).
	(Typical adapter trimming therefore requires the sequences to be
	passed in in reverse.)  It returns the score of the resulting
	alignment, which will be positive, or -1 if the alignment would be
	too bad.  The number of characters in sequence A that are overlapped
	is returned in *yout.

	\param seq_a iterator to beginning of sequence A.
	\param seq_a_end iterator to end of sequence A.
	\param seq_b iterator to beginning of sequence B.
	\param seq_b_end iterator to end of sequence B.
	\param yout pointer to location where number of trimmed chars is stores.
	\return alignment score.
*/
	
template< typename A, typename B >
int overlap_align( A seq_a, A seq_a_end, B seq_b, B seq_b_end, int*yout )
{
	static const int discount = 8 ;
	using namespace std ;

	int len_a = seq_a_end - seq_a ;
	int len_b = seq_b_end - seq_b ;
	int maxd = 2 * len_b / discount ;						// makes sure that score>=0
    int dim = discount*maxd + discount-1 ;

	vector<int> v_d( dim+maxd+1 ), v_dm1( dim+maxd+1 ) ;	// v[d] & v[d-1]

	for( int d = 0 ; d != maxd ; ++d, swap( v_d, v_dm1 ) ) {			// D-paths in order of increasing D
		int dm = discount*d + discount-1 ;
		for( int k = -d ; k <= dm ; ++k ) {								// diagonals
			int x = 0 ==  d            ? k                                                                     :
				    1 ==  d && 0 == k  ?                             v_dm1[  k +maxd ]+1					   :
			        k == -d            ?                                                  v_dm1[ k+1 +maxd ]   :
			        k  > dm-discount+1 ? k                                                                     :
			        k == -d+1          ? max(                        v_dm1[  k +maxd ]+1, v_dm1[ k+1 +maxd ] ) :
			        k == dm-discount+1 ? max(  v_dm1[ k-1 +maxd ]+1,                      k )                  :
			        k == dm-discount   ? max3( v_dm1[ k-1 +maxd ]+1, v_dm1[  k +maxd ]+1, k )                  :
			                             max3( v_dm1[ k-1 +maxd ]+1, v_dm1[  k +maxd ]+1, v_dm1[ k+1 +maxd ] ) ;

			int y = x - k ;
			while( x < len_b && y < len_a && match( seq_b[x], seq_a[y] ) ) ++x, ++y ;
			v_d[ k +maxd ] = x ;

			if( x == len_b )
			{
				*yout = y ;
				return x + y - discount * d ; // computes score.
			}
		}
	}
	return -1 ;
}

#endif
