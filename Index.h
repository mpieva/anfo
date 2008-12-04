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

#ifndef _BSD_SOURCE
inline int madvise( void*, size_t, int ) { return 0 ; }
static const int MADV_SEQUENTIAL = 0 ;
static const int MADV_WILLNEED = 0 ;
#endif

class DnaP
{
	private:
		uint8_t const *base ;

	public:
		explicit DnaP( void const *p = 0 ) : base(static_cast<uint8_t const*>(p)) {}
		Ambicode operator[]( uint32_t ix ) const {
			uint8_t w = base[ix >> 1] ;
			if( ix & 1 ) w >>= 4 ;
			return w & 0xf ;
		}
		operator void const * () const { return base ; }
		void assign( void const *p ) { base = static_cast<uint8_t const*>(p) ; }

		void *unsafe_ptr() const { return const_cast<void*>( static_cast<const void*>( base ) ) ; }
} ;

class CompactGenome
{
	public:
		static const uint32_t signature = 0x30414e44 ; // DNA0 

		CompactGenome( const char* fp ) ;
		~CompactGenome() ;

		DnaP get_base() const { return base ; }

		// scan over finite words of the dna; these may well contain
		// ambiguity codes, but are guaranteed not to contain gaps
		template< typename F, typename G >
			void scan_words( unsigned w, F mk_word, G& consume_word, const char* msg = 0 ) ;

	private:
		DnaP base ;
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

		FixedIndex( const char* fp, unsigned w, unsigned cutoff = std::numeric_limits<unsigned>::max() ) ;
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
		unsigned cutoff ;
} ;

// Scan over the dna words in the genome.  Size of the words is limited
// to what fits into a (long long unsigned), typically 16.
//
// Words are encoded as four bits per nucleotide, first nucleotide in
// the MSB(!).  See make_dense_word for why that makes sense.  Unused
// MSBs in words passed to mk_word contain junk.
template< typename F, typename G >
void CompactGenome::scan_words( unsigned w, F mk_word, G& consume_word, const char* msg ) {
	assert( (unsigned)std::numeric_limits< Oligo >::digits >= 4 * w ) ;
	madvise( base.unsafe_ptr(), length, MADV_SEQUENTIAL ) ;

	uint32_t offs = 0 ;
	Oligo dna = 0 ;

	while( base[ offs ] != 0 ) ++offs ;		// first first gap
	for( unsigned i = 0 ; i != w ; ++i )				// fill first word
	{
		dna <<= 4 ;
		dna |= base[ offs ] ;
		++offs ;
	}

	report(offs,length,msg) ;
	while( offs != 2 * length )
	{
		if( (offs & 0xfffff) == 0 ) report(offs,length,msg) ;

		dna <<= 4 ;
		dna |= base[ offs ] ;
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
// has enough seeds is used.  Seeds remaining after this process will be
// grouped at the beginning of their container, the rest is erased.
//
// In a sense, this wants to be a single, configurable function.  Right
// now configuration is simply done by three additional parameters:
// Seeds on the same diagonal ±d and not further apart than ±r are
// clumped.  A clump is good for an alignment if the total match length
// of the seeds in it reaches m, and if so, its longest seed is actually
// used.  We may need different parameters for input sequences of
// different lengths or characteristics; that decision, however, has to
// be deferrred to some higher abstraction level.
//
// The container type used in the following must be a random access
// container.  Either std::vector or std::deque should be fine.
template < typename C > void select_seeds( C& v, uint32_t d, int32_t r, uint32_t m )
{
	if( !v.empty() )
	{
		std::sort( v.begin(), v.end(), compare_diag_then_offset() ) ;
		typename C::iterator clump_begin = v.begin(),
				 input_end = v.end(), out = v.begin() ;

		// Start building a clump, assuming there's is still something
		// to build from
		while( clump_begin != input_end )
		{
			typename C::iterator clump_end = clump_begin + 1 ;
			typename C::iterator last_touched = clump_end ;

			// Still anything in the clump that may have unrecognized
			// neighbors?
			for( typename C::iterator open_in_clump = clump_begin ; open_in_clump != clump_end ; ++open_in_clump )
			{
				// Decide whether open_in_clump and candidate are
				// actually neighbors.  XXX: They are not if they end up
				// in different contigs :XXX; else they are if their
				// diagonals are dloser than ±d and their offsets are
				// closer than ±r.
				for( typename C::iterator candidate = clump_end ;
						candidate != input_end &&
						candidate->diagonal <= open_in_clump->diagonal + d ;
						++candidate )
				{
					if( abs( candidate->offset - open_in_clump->offset ) <= r
							&& (candidate->offset>=0) == (open_in_clump->offset>=0) ) 
					{
						// Include the candidate by swapping it with the
						// first seed not in our clump and extending the
						// clump.  Remember that we swapped, we may have
						// messed up the sorting.
						std::swap( *clump_end++, *candidate ) ;
						if( candidate < last_touched ) last_touched = candidate ;
					}
				}
			}

			// Okay, we have built our clump, but we left a mess behind
			// it (between clump_end and last_touched).  Clean up the
			// mess first.
			if( clump_end < last_touched ) 
				std::sort( clump_end, last_touched, compare_diag_then_offset() ) ;
			
			// The new clump sits between clump_begin and clump_end.  We
			// will now reduce this to at most one "best" seed.
			typename C::iterator best = clump_begin ;
			uint32_t mlen = 0 ;
			for( typename C::iterator cur = clump_begin ; cur != clump_end ; ++cur )
			{
				mlen += cur->size ;
				if( cur->size > best->size ) best = cur ;
			}

			// The whole clump is good enough if the total match length
			// reaches m.  If so, we keep its biggest seed by moving it
			// to out, else we do nothing.
			if( mlen >= m ) *out++ = *best ;

			// Done with the clump.  Move to the next.
			clump_begin = clump_end ;
		}

		// We updated the input vector from begin() to out, the rest is
		// just junk that's left over.  Get rid of it.
		v.erase( out, input_end ) ;
	}
}

#endif

