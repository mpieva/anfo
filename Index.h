#ifndef INCLUDED_INDEX_H
#define INCLUDED_INDEX_H

#include <cassert>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <stdint.h>
#include <sys/mman.h>

// Useful typedefs.  These are used mostly for their documentation
// value, C++ unfortunately won't be able to check their differences
// most of the time.

typedef uint64_t Oligo ;
typedef uint8_t  Nucleotide ;	// 0,1,2,3 == A,C,T,G
typedef uint8_t  Ambicode ;		// 1,2,4,8 == A,C,T,G

class CompactGenome
{
	public:
		static const uint32_t signature = 0x30414e44 ; // DNA0 

		CompactGenome( const char* fp ) ;
		~CompactGenome() ;

		// get nucleotide at position ix, relative to the beginning of
		// the file
		Ambicode operator[]( uint32_t ix ) {
			uint8_t w = base[ix >> 1] ;
			if( ix & 1 ) w >>= 4 ;
			return w & 0xf ;
		}

		// scan over finite words of the dna; these may well contain
		// ambiguity codes, but are guaranteed not to contain gaps
		template< typename F, typename G >
			void scan_words( unsigned w, F mk_word, G consume_word, const char* msg = 0 ) ;

	private:
		uint8_t *base ;
		uint32_t length ;
		int fd ;

		static void report( unsigned, unsigned, const char* ) ;
} ;

// What do we need to know about a seed?  A size and two coordinates.
// One is the start on the query sequence, the other we choose to be the
// "diagonal", since seeds on the same or close diagonal will usually be
// combined.  Offset is negative for rc'ed matches.
struct Seed
{
	uint32_t diagonal ;
	uint32_t size ;
	int32_t offset ;
} ;

class FixedIndex 
{
	public:
		enum { signature = 0x30584449u } ; // IDX0 

		FixedIndex( const char* fp, unsigned w ) ;
		~FixedIndex() ;

		// direct lookup of an oligo, results are placed in the vector,
		// the number of results is returned.
		unsigned lookup1( Oligo, std::vector<Seed>&, int32_t offset = 0 ) const ;

		// Lookup of a sequence.   The sequence is split into words as
		// appropriate for the index, then each one of them is looked
		// up.  This method can be implemented for any kind of index.
		unsigned lookup( const std::string& seq, std::vector<Seed>& ) const ;

	private:
		uint32_t *base, *secondary ;
		uint32_t first_level_len ;
		uint64_t length ;
		int fd ;
		unsigned wordsize ;
} ;

// Scan over the dna words in the genome.  Size of the words is limited
// to what fits into a (long long unsigned), typically 16.
//
// Words are encoded as four bits per nucleotide, first nucleotide in
// the MSB(!).  See make_dense_word for why that makes sense.  Unused
// MSBs in words passed to mk_word contain junk.
template< typename F, typename G >
void CompactGenome::scan_words( unsigned w, F mk_word, G consume_word, const char* msg ) {
	assert( (unsigned)std::numeric_limits< Oligo >::digits >= 4 * w ) ;
#ifdef _BSD_SOURCE
	madvise( base, length, MADV_SEQUENTIAL ) ;
#endif
	uint32_t offs = 0 ;
	Oligo dna = 0 ;

	while( (*this)[ offs ] != 0 ) ++offs ;		// first first gap
	for( unsigned i = 0 ; i != w ; ++i )				// fill first word
	{
		dna <<= 4 ;
		dna |= (*this)[ offs ] ;
		++offs ;
	}

	report(offs,length,msg) ;
	while( offs != 2 * length )
	{
		if( (offs & 0xfffff) == 0 ) report(offs,length,msg) ;

		dna <<= 4 ;
		dna |= (*this)[ offs ] ;
		++offs ;

		// throw away words containing gap symbols
		// (This is necessary since we may want to construct
		// discontiguous words, but not if a "don't care" position is a
		// gap.)
		bool clean = true ;
		for( unsigned i = 0 ; i != w ; ++i )
			clean &= ((dna >> (4*i)) & 0xf) != 0 ;

		if( clean ) mk_word( w, offs-w, dna, consume_word ) ;
	}
	if( msg ) std::clog << "\r\e[K" << std::flush ;
}

// Combine short, adjacent seeds into longer ones.  The exact policy for
// that isn't quite clear yet, but what is clear is that we can always
// combine directly adjacent seeds.  It is also guaranteed that no such
// seeds span multiple contigs, so this is absolutely safe.
//
// How to do this?  We sort seeds first by diagonal index, then by
// offset.  Seeds are adjacent iff they have the same diagonal index and
// their offsets differ by no more than the word size.
struct compare_diag_then_offset {
	bool operator()( const Seed& a, const Seed& b ) {
		if( a.diagonal < b.diagonal ) return true ;
		if( b.diagonal < a.diagonal ) return false ;
		return a.offset < b.offset ;
	}
} ;

template < typename C > void combine_seeds( C& v ) 
{
	if( !v.empty() )
	{
		std::sort( v.begin(), v.end(), compare_diag_then_offset() ) ;
		typename C::const_iterator a = v.begin(), e = v.end() ;
		typename C::iterator       d = v.begin() ;
		Seed s = *a ; 
		while( ++a != e )
		{
			if( a->diagonal == s.diagonal )
			{
				if( a->offset - s.offset <= (int32_t)s.size ) {
					uint32_t size_ = a->offset - s.offset + a->size ;
					if( size_ > s.size ) s.size = size_ ;
				}
			}
			else
			{
				*d = s ; ++d ;
				s = *a ;
			}
		}
		*d = s ; ++d ;
		v.erase( d, v.end() ) ;
	}
}

// How do we select seeds?  Seeds that are "close enough" should be
// collapsed into a cluster, then the best seed from any cluster that
// has enough seeds is used.  
//
// XXX: This wants to be configurable.  Seeds on the same diagonal ±d
// and not further apart than ±r are clumped.  A clump is good enough a
// seed if the total match length of the seeds in it reaches m.
template < typename C > void select_seeds( C& v, uint32_t d, int32_t r, uint32_t m )
{
	if( !v.empty() )
	{
		std::sort( v.begin(), v.end(), compare_diag_then_offset() ) ;
		typename C::iterator dest = v.begin(), e = v.end() ;

		for( typename C::const_iterator a = v.begin() ; a != e ; )
		{
			// find range of seeds that form a clump, this will be
			// pointed to by [a,b).  Also calculate total seed length.
			typename C::const_iterator b = a ;
			uint32_t mlen = 0 ;
			while( b != e && b->diagonal - a->diagonal <= d
					&& abs( b->offset - a->offset ) <= r
					&& (b->offset>=0) == (a->offset>=0) ) 
			{
				mlen += b->size ;
				++b ;
			}

			// if good enough, find the longest seed in the clump and
			// keep it
			if( mlen > m ) {
				typename C::const_iterator c = a ;
				for( typename C::const_iterator i = a ; i != b ; ++i )
					if( i->size > c->size ) c = i ;

				*dest = *c ; ++dest ;
			}
			a = b ;
		}
		v.erase( dest, e ) ;
	}
}

#endif

