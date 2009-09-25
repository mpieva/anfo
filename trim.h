//    Copyright 2009 Udo Stenzel
//    This file is part of ANFO
//
//    ANFO is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    Foobar is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

#ifndef INCLUDED_TRIM_H
#define INCLUDED_TRIM_H

#include <algorithm>
#include <vector>

enum Mode { equal, a_is_prefix, b_is_prefix, either_prefix, overlap, primer_special } ;

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

/* Alignment algorithm following Eugene W. Myers: "An O(ND) Difference
 * Algorithm and Its Variations".
 * seq_a and seq_b are input sequences, null-terminated.  maxd is the
 * maximum edit distance to consider.  Set gap_a to true to get a free
 * end gap in seq_a, similarly for gap_b and seq_b.  If you want
 * backtracing, point bt to a buffer of at least
 * (strlen(seq_a)+strlen(seq_b)+1)*2 bytes.
 *
 * [*blech*, this feels like FORTRAN.]
 */

template < typename A, typename B >
int align( A seq_a, A seq_a_end, B seq_b, B seq_b_end, int maxd, Mode mode, /* char*bt, */ int*xout, int*yout )
{
    int ntargets = 1 ;
    return alignN( seq_a, seq_a_end, &seq_b, &seq_b_end, maxd, &ntargets, mode, /*bt,*/ xout, yout ) ;
}

/*
	Regarding parallel alignment:
    
    To make this tractable, the monolithic algorithm needs to be split.
    An object is to be used to simulate a closure, most local variables
    end up as members, the loop over d is taken out and the result type
    is converted to return either a result, a failure notice, or nothing
    at all (using a union).  The loop over d is done separately, outside
    the core algorithm.  The d-loop should be able to iterate over many
    such objects and objects drop out of the iteration early if they
    fail.

TODO: To keep memory consumption down, it may be worthwhile to implement
    this without backtracing; a second pass can always be done to get
    the actual alignment.  In fact, this may be sensible in any case,
    parallel alignment or not.  Implementing the O(d) space variant
    (twice the calculation effort, but backtracing becomes simple
    recursion) might actually yield simpler code anyway.
*/

/*
	Regarding overlap alignment:

	We want to align two sequences so they overlap with maximum score.
	This is not possible with the plain "minimum distance" method, since
	no overlap always results in the minimum distance of zero.  We do
	need to score the overlapping part, but not as ordinary difference.
	Now suppose you have an alignment:

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
    (not implemented here), where only an X-drop algorithm is
    appropriate).

	To implement this, the following changes are necessary:
	- allow a start gap that grows 8x quicker than a normal gap,
	  consider the many more involved diagonals,
    - once finished, prefer alignments on lower numbered diagonals.


	XXX: Backtracing of an overlap-alignment doesn't currently work.
*/
	
static const int discount = 8 ;

enum AlignResultTag { art_progressing, art_failed, art_finished } ;

template< typename A, typename B >
struct AlignClosure
{
	A seq_a ;
	B seq_b ;
	int len_a, len_b ;
	int maxd, id ;
	Mode mode ;

	std::vector<int> v_ ;

	AlignClosure( A seq_a_, int len_a_, B seq_b_, int len_b_, int maxd_, Mode mode_, int id_ )
		: seq_a( seq_a_ ), seq_b( seq_b_ )
		, len_a( len_a_ ), len_b( len_b_ )
		, maxd( std::min( maxd_, len_a_ + len_b_ ) )
		, id( id_ ), mode( mode_ ), v_( (1+maxd_) * ((mode==overlap?discount+1:2)*maxd_+1) )
	{ }

	int& v( int k, int d ) { return v_[ k + maxd + ((mode==overlap?discount+1:2)*maxd+1) * d ] ; }

	AlignResultTag operator () ( int d, /*char* bt,*/ int*xout, int*yout ) ;
} ;


	
template < typename A, typename B, typename C >
AlignClosure<A,B> *runMany( C& closures, int& maxd, /*char* bt,*/ int*xout, int*yout )
{
	for( int d = 0 ; d != maxd && !closures.empty() ; ++d )		// D-paths in order of increasing D
	{
		for( typename C::iterator i = closures.begin() ; i != closures.end() ; )
		{
			AlignResultTag tag = (*i)( d, /* bt,*/ xout, yout ) ;
			if(      tag == art_finished ) { maxd = d ; return &(*i) ; }
			else if( tag == art_failed   ) i = closures.erase( i ) ;
			else ++i ;
		}
	}
	return 0 ;
}

template< typename A, typename B >
AlignResultTag AlignClosure<A,B>::operator () ( int d, /*char* bt,*/ int*xout, int*yout )
{
	using std::max ;
    int dm = mode == overlap ? discount*d + discount-1 :
             mode == primer_special ? d + 2 : d ;

	for( int k = -d ; k <= dm ; ++k ) // diagonals
	{
		int x ;
		if( mode != overlap )
		{
			if     ( 0 ==  d   ) x = k ;
			else if( k == -d   ) x =                                        v( k+1, d-1 )   ;
			else if( k == -d+1 ) x = max(                   v(  k, d-1 )+1, v( k+1, d-1 ) ) ;
			else if( k == dm   ) x =       v( k-1, d-1 )+1                                  ;
			else if( k == dm-1 ) x = max(  v( k-1, d-1 )+1, v(  k, d-1 )+1                ) ;
			else                 x = max3( v( k-1, d-1 )+1, v(  k, d-1 )+1, v( k+1, d-1 ) ) ;
		}
		else
		{
			if     ( 0 ==  d   )          x = k ;
			else if( k == -d   )          x =                                        v( k+1, d-1 )   ;
			else if( k == -d+1 )          x = max(                   v(  k, d-1 )+1, v( k+1, d-1 ) ) ;
			else if( k  > dm-discount+1 ) x = k ;
			else if( k == dm-discount+1 ) x = max(  v( k-1, d-1 )+1,                 k ) ;
			else if( k == dm-discount   ) x = max3( v( k-1, d-1 )+1, v(  k, d-1 )+1, k ) ;
			else                          x = max3( v( k-1, d-1 )+1, v(  k, d-1 )+1, v( k+1, d-1 ) ) ;
		}

		int y = x-k ;
		while( x < len_b && y < len_a && match( seq_b[x], seq_a[y] ) ) ++x, ++y ;
		v(k,d) = x ;

		if( 
				(mode == equal && y == len_a && x == len_b) ||
				(mode == a_is_prefix && y == len_a) ||
				(mode == b_is_prefix && x == len_b) ||
				(mode == primer_special && x == len_b) ||
				(mode == either_prefix && (y == len_a || x == len_b)) || 
                (mode == overlap && x == len_b)
		  )
		{
			if( xout ) *xout = x ;
			if( yout ) *yout = y ;
			/* if( bt )
			{
				char *out = bt + 2 * (len_a + len_b +1) ;
				char *out_end = bt + 2 * (len_a + len_b +1) ;
				*--out = 0 ;
				*--out = 0 ;
				for( int dd = d ; dd != 0 && y != 0 && x != 0 ; )
				{
					if( k != -dd && k != dd && x == v( k, dd-1 )+1 )
					{
						--dd ;
						--x ;
						--y ;
						*--out = seq_b[x] ;
						*--out = seq_a[y] ;
					}
					else if( k > -dd+1 && x == v( k-1, dd-1 )+1 )
					{
						--x ;
						--k ;
						--dd ;
						*--out = seq_b[x] ;
						*--out = '-' ;
					}
					else if( k < dd-1 && x == v( k+1, dd-1 ) )
					{
						++k ;
						--y ;
						--dd ;
						*--out = '-' ;
						*--out = seq_a[y] ;
					}
					else
					{
						--x ;
						--y ;
						*--out = seq_b[x] ;
						*--out = seq_a[y] ;
					}
				}
				while( y > 0 && x > 0 )
				{
					--y ;
					--x ;
					*--out = seq_b[x] ;
					*--out = seq_a[y] ;
				}
				while( x > 0 )
				{
					--x ;
					*--out = seq_b[x] ;
					*--out = '-' ;
				}
				while( y > 0 )
				{
					--y ;
					*--out = '-' ;
					*--out = seq_a[y] ;
				}
				memmove( bt, out, out_end - out ) ;
			} */
			return art_finished ;
		}
	}
	return art_progressing ;
}

// same as above, parallelized on many target sequences
// ntarget points to the number of targets initially and contains the
// best aligning target index finally.  seq_b points to an array of
// ntarget pointers to null-terminated sequences.
template< typename A, typename B >
int alignN( A seq_a, A seq_a_end, B *seqs_b, B *seqs_b_end, int maxd, int* ntargets,
		Mode mode, /*char*bt,*/ int*xout, int*yout )
{
	typedef std::vector< AlignClosure< A, B > > C ;
	C closures ;
	for( int i = 0 ; i != *ntargets ; ++i )
		closures.push_back( AlignClosure<A,B>(
					seq_a, seq_a_end - seq_a,
					seqs_b[i], seqs_b_end[i] - seqs_b[i],
					maxd, mode, i ) ) ;

	AlignClosure<A,B> *r = runMany<A,B,C>( closures, maxd, /*bt,*/ xout, yout ) ;
	if( r ) *ntargets = r->id ;
	return maxd ;
}

#endif
