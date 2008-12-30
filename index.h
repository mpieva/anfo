#ifndef INCLUDED_INDEX_H
#define INCLUDED_INDEX_H

#include "sequence.h"
#include "util.h"

#include <metaindex.pb.h>

#include <cassert>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

//! \brief a genome as stored in a DNA file
//! This class mmaps a DNA file and wraps it with a sensible interface.
class CompactGenome
{
	public:
		typedef std::map< uint32_t, std::pair< const metaindex::Sequence*, const metaindex::Contig* > > ContigMap ;

	private:
		DnaP base_ ;
		uint32_t length_ ;
		int fd_ ;
		ContigMap contig_map_ ;

	public:
		//! \brief constructs an invalid genome
		//! Genomes constructed in the default fashion are unusable;
		//! however, this makes \c CompactGenome default constructible
		//! for use in standard containers.
		CompactGenome() : base_(0), length_(0), fd_(0) {}

		//! \brief makes accessible a file described by metadata
		//! \param g genome definition as found in a configuration file
		//! \param adv advise passed to \c madvise, if you anticipate
		//!            specific use of the genome
		CompactGenome( const metaindex::Genome &g, int adv = MADV_NORMAL ) ;
		~CompactGenome() ;

		DnaP get_base() const { return base_ ; }

		//! \brief scan over finite words of the dna
		//! Fixed size words in the genome are iterated.  The words may
		//! well contain ambiguity codes, but are guaranteed not to
		//! contain gaps.
		//!
		//! Words are encoded as four bits per nucleotide, first nucleotide in
		//! the MSB(!).  See ::make_dense_word for why that makes sense.  Unused
		//! MSBs in words passed to mk_word contain junk, do not rely on
		//! them.
		//!
		//! The functor objects in the following will be passed by
		//! value.  Make sure they are tiny and that their function call
		//! operators can be inlined.
		//!
		//! \param w word size, 4*w must not be more than there are bits
		//!          in an unsigned long long.
		//! \param mk_word functor object called to transform a
		//!                sequence of ambiguity codes into words of
		//!                nucleotides
		//! \param msg if set, switches on progress reports and is
		//!            included in them
		template< typename F > void scan_words( unsigned w, F mk_word, const char* msg = 0 ) ;

		void swap( CompactGenome& g )
		{
			std::swap( base_, g.base_ ) ;
			std::swap( length_, g.length_ ) ;
			std::swap( fd_, g.fd_ ) ;
			std::swap( contig_map_, g.contig_map_ ) ;
		}

		const ContigMap &get_contig_map() const { return contig_map_ ; }

	private:
		static const uint32_t signature = 0x30414e44 ; // DNA0 

		//! \brief reports a position while scanning
		//! \internal
		static void report( uint32_t, uint32_t, const char* ) ;
} ;

//! \brief Representation of a seed.
//! A seed is described by a size and two coordinates.  One is the start
//! on the query sequence, the other we choose to be the "diagonal"
//! (difference between coordinates), since seeds on the same or close
//! diagonal will usually be combined.  Offset is negative for rc'ed
//! matches, in this case its magnitude is the actual offset from the
//! end of the sequence.
struct Seed
{
	uint32_t diagonal ;
	uint32_t size ;
	int32_t offset ;
} ;

inline std::ostream& operator << ( std::ostream& o, const Seed& s )
{
	return o << '@' << s.offset << '+' << s.diagonal << ':' << s.size ;
}

class FixedIndex 
{
	public:
		enum { signature = 0x30584449u } ; // IDX0 

		FixedIndex() : base(0), secondary(0), first_level_len(0), length(0), fd(0), wordsize(0) {}
		FixedIndex( const char* fp, unsigned w ) ;
		~FixedIndex() ;

		unsigned lookup1( Oligo, std::vector<Seed>&, uint32_t cutoff, int32_t offset = 0 ) const ;
		unsigned lookup( const QSequence& seq, std::vector<Seed>&,
				uint32_t cutoff = std::numeric_limits<uint32_t>::max() ) const ;

		operator const void * () const { return base ; }

		void swap( FixedIndex& i ) {
			std::swap( base, i.base ) ;
			std::swap( secondary, i.secondary ) ;
			std::swap( first_level_len, i.first_level_len ) ;
			std::swap( length, i.length ) ;
			std::swap( fd, i.fd ) ;
			std::swap( wordsize, i.wordsize ) ;
		}

	private:
		uint32_t *base, *secondary ;
		uint32_t first_level_len ;
		uint64_t length ;
		int fd ;
		unsigned wordsize ;
} ;

template< typename F > void CompactGenome::scan_words( unsigned w, F mk_word, const char* msg )
{
	assert( (unsigned)std::numeric_limits< Oligo >::digits >= 4 * w ) ;
	madvise( (void*)base_.unsafe_ptr(), length_, MADV_SEQUENTIAL ) ;

	uint32_t offs = 0 ;
	Oligo dna = 0 ;

	while( base_[ offs ] != 0 ) ++offs ;	// find first gap
	for( unsigned i = 0 ; i != w ; ++i )	// fill first word
	{
		dna <<= 4 ;
		dna |= base_[ offs ] ;
		++offs ;
	}

	report(offs,length_,msg) ;
	while( offs != 2 * length_ )
	{
		if( (offs & 0xfffff) == 0 ) report(offs,length_,msg) ;

		dna <<= 4 ;
		dna |= base_[ offs ] ;
		++offs ;

		// throw away words containing gap symbols
		// (This is necessary since we may want to construct
		// discontiguous words, but not if a "don't care" position is a
		// gap.)
		bool clean = true ;
		for( unsigned i = 0 ; i != w ; ++i )
			clean &= ((dna >> (4*i)) & 0xf) != 0 ;

		if( clean ) mk_word( w, offs-w, dna ) ;
	}
	if( msg ) std::clog << "\r\e[K" << std::flush ;
}

//! \brief compares seeds first by diagonal, then by offset
//! \internal
//! Functor object passed to std::sort in various places.
struct compare_diag_then_offset {
	bool operator()( const Seed& a, const Seed& b ) {
		if( a.diagonal < b.diagonal ) return true ;
		if( b.diagonal < a.diagonal ) return false ;
		return a.offset < b.offset ;
	}
} ;

//! \brief Combines short, adjacent seeds into longer ones.
//! The exact policy for aggregating seeds isn't quite clear yet, but
//! what is clear is that we can always combine directly adjacent or
//! overlapping seeds.  It is also guaranteed that no such seeds span
//! multiple contigs, so this is absolutely safe.
//!
//! How to do this?  We sort seeds first by diagonal index, then by
//! offset.  Seeds are adjacent iff they have the same diagonal index and
//! their offsets differ by no more than the seed size.
//!
//! \param v container of seeds, will be modifed in place.
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
				if( (a->offset >= 0) == (s.offset >= 0) &&
						a->offset - s.offset <= (int32_t)s.size ) {
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

//! \brief aggregates and selects seeds
//! How to select seeds isn't finalized, however, the following seems to
//! be reasonable:  Seeds that are "close enough" should be collapsed
//! into a clump, then the best seed from any clump that has enough
//! seeds is used.  Seeds remaining after this process will be grouped
//! at the beginning of their container, the rest is erased.
//!
//! Configuration of this function is simply done by three additional
//! parameters: Seeds on the same diagonal ±d and not further apart than
//! ±r are clumped.  A clump is good enough for an alignment if the
//! total match length of the seeds in it reaches m, and if so, its
//! longest seed is actually used.  We may need different parameters for
//! input sequences of different lengths or characteristics; that
//! decision, however, has to be deferrred to some higher abstraction
//! level.
//!
//! The container type used in the following must be a random access
//! container.  Either std::vector or std::deque should be fine.
//!
//! \param v container of seeds, will be modified in place.
//! \param d maximum number of gaps between seeds to be clumped
//! \param r maximum unseeded length between seeds to be clumped
//! \param m minimum total seed length in a good clump
//! \param cm contig map from indexed genome
template < typename C > void select_seeds( C& v, uint32_t d, int32_t r, uint32_t m,
		const CompactGenome::ContigMap &cm )
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
			// Decide whether open_in_clump and candidate are actually
			// neighbors.  They are not if they end up in different
			// contigs; else they are if their diagonals are dloser than
			// ±d and their offsets are closer than ±r.
			for( typename C::iterator candidate = clump_end ;
					candidate != input_end &&
					candidate->diagonal <= open_in_clump->diagonal + d ;
					++candidate )
			{
				if( abs( candidate->offset - open_in_clump->offset ) <= r
						&& (candidate->offset>=0) == (open_in_clump->offset>=0) )
				{
					CompactGenome::ContigMap::const_iterator low = cm.lower_bound(
							open_in_clump->offset + open_in_clump->diagonal + open_in_clump->size ) ;
					CompactGenome::ContigMap::const_iterator high = cm.upper_bound(
							candidate->offset + candidate->diagonal ) ;

					// make sure no contig start is in between
					if( low == high ) {
						// Include the candidate by swapping it with the
						// first seed not in our clump and extending the
						// clump.  Remember that we swapped, we may have
						// messed up the sorting.
						std::swap( *clump_end++, *candidate ) ;
						if( candidate < last_touched ) last_touched = candidate ;
					}
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

#endif

