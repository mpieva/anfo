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

#ifndef INCLUDED_INDEX_H
#define INCLUDED_INDEX_H

#include "config.pb.h"
#include "sequence.h"
#include "util.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <glob.h>

//! \brief a genome as stored in a DNA file
//! This class mmaps a DNA file and wraps it with a sensible interface.
class CompactGenome
{
	public:
		// Why are we using indices here instead of pointers?  Because
		// pointers would be internal and make copying of this object a
		// lot harder.

		// the pair is (index of sequence, index of contig)
		typedef std::map< uint32_t, std::pair< int, int > > ContigMap ;
		
		// value is index of contig
		typedef std::map< uint32_t, int > PosnMap1 ;

		// first half of value is index of sequence
		typedef std::map< std::string, std::pair< int, PosnMap1 > > PosnMap ;

	private:
		DnaP base_ ;
		size_t file_size_ ;
		uint32_t length_ ;
		int fd_ ;
		ContigMap contig_map_ ;
		PosnMap posn_map_ ;

		CompactGenome( const CompactGenome& ) ; // not implemented
		void operator = ( const CompactGenome& ) ; // not implemented

	public:
		config::Genome g_ ;

	public:
		//! \brief constructs an invalid genome
		//! Genomes constructed in the default fashion are unusable;
		//! however, this makes \c CompactGenome default constructible
		//! for use in standard containers.
		CompactGenome() : base_(0), file_size_(0), length_(0), fd_(-1), contig_map_(), posn_map_(), g_() {}

		//! \brief makes accessible a genome file
		//! \param name file name of the genome
		//! \param c program configuration, needed for the search path
		//! \param adv advise passed to madvise(), if you anticipate
		//!            specific use of the genome
		CompactGenome( const std::string& name, int adv = MADV_NORMAL ) ;
		~CompactGenome() ;

		std::string name() const { return g_.name() ; }
		std::string describe() const 
		{
			std::string d = g_.name() ;
			if( g_.has_description() ) d += " (" + g_.description() + ")" ;
			return d ;
		}

		uint32_t total_size() const { return g_.total_size() ; }
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
			std::swap( file_size_, g.file_size_ ) ;
			std::swap( length_, g.length_ ) ;
			std::swap( fd_, g.fd_ ) ;
			std::swap( contig_map_, g.contig_map_ ) ;
			std::swap( posn_map_, g.posn_map_ ) ;
			std::swap( g_, g.g_ ) ;
		}

		const ContigMap &get_contig_map() const { return contig_map_ ; }

		//! \brief translates a DNA pointer back to sequence coordinates
		//! If the DNA pointer points into this genome, it is translated
		//! to the name of the sequence and the offset into it.  The
		//! strand is disregarded, the result is the same regardless of
		//! the direction the pointer would be moving in.
		//! \param pos position to be translated
		//! \param offset is assigned the offset after translation
		//! \return pointer to sequence that includes the hit, else null
		const config::Sequence *translate_back( DnaP pos, uint32_t& offset ) const ;

		//! \brief translates "human" cooordinates into a pointer
		//! Return 0 if anything goes wrong.
		DnaP find_pos( const std::string& seq, uint32_t pos ) const ;

		enum { signature = 0x31414e44u } ; // DNA1 

	private:
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
	int64_t diagonal ;
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
		enum { old_signature = 0x31584449u     // IDX1 
		     , signature     = 0x32584449u } ; // IDX2

		FixedIndex() : p_(0), base(0), secondary(0), first_level_len(0), length(0), fd_(0), ci_() {}

		//! \brief loads an index from a file
		//! \param name filename
		//! \param c I'll be damned if I remember what this was for
		//! \param adv advice passed to mmap()
		FixedIndex( const std::string &name, const config::Config &c, int adv = MADV_NORMAL ) ;
		~FixedIndex() ;

		unsigned lookupS( const std::string&  seq, std::vector<Seed>&,
				bool near_perfect = false, int *num_useless = 0,
				uint32_t cutoff = std::numeric_limits<uint32_t>::max() ) const ;
		unsigned lookup1( Oligo, std::vector<Seed>&, uint32_t cutoff, int32_t offs, int *num_useless ) const ;
		unsigned lookup1m( Oligo o, std::vector<Seed>& v, uint32_t cutoff, int32_t offs, int *num_useless ) const ;

		operator const void * () const { return base ; }

		void swap( FixedIndex& i ) {
			std::swap( p_, i.p_ ) ;
			std::swap( base, i.base ) ;
			std::swap( secondary, i.secondary ) ;
			std::swap( first_level_len, i.first_level_len ) ;
			std::swap( length, i.length ) ;
			std::swap( fd_, i.fd_ ) ;
			std::swap( ci_, i.ci_ ) ;
		}

	private:
		const void* p_ ;
		const uint32_t *base, *secondary ;
		uint32_t first_level_len ;
		uint64_t length ;
		int fd_ ;

	public:
		config::CompactIndex ci_ ;
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
			if( a->diagonal == s.diagonal &&
					(a->offset >= 0) == (s.offset >= 0) &&
					a->offset - s.offset <= (int32_t)s.size )
			{
				uint32_t size2 = a->offset - s.offset + a->size ;
				if( size2 > s.size ) s.size = size2 ;
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
//! parameters: Seeds on the same diagonal �d and not further apart than
//! �r are clumped.  A clump is good enough for an alignment if the
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
			// �d and their offsets are closer than �r.
			for( typename C::iterator candidate = clump_end ;
					candidate != input_end &&
					candidate->diagonal <= open_in_clump->diagonal + d ;
					++candidate )
			{
				if( abs( candidate->offset - open_in_clump->offset ) <= r
						&& (candidate->offset>=0) == (open_in_clump->offset>=0) )
				{
					// make sure both parts belong to the same contig
					CompactGenome::ContigMap::const_iterator left = cm.upper_bound(
							open_in_clump->offset + open_in_clump->diagonal ) ;
					CompactGenome::ContigMap::const_iterator right = cm.upper_bound(
							candidate->offset + candidate->diagonal ) ;

					if( left == right ) {
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

typedef std::map< std::string, CompactGenome > Genomes ;
typedef std::map< std::string, FixedIndex > Indices ;

class Metagenome
{
	public:
		enum Persistence { persistent = 0, ephemeral = 1 } ;
		static bool nommap ;
		static void make_room() ;

	private:
		// maps file name to genome object
		typedef std::map< std::string, std::pair< Persistence, CompactGenome* > > Genomes ;

		// maps sequence id to genome object
		typedef std::map< std::string, CompactGenome* > SeqMap1 ;

		// maps genome id to sequence map
		typedef std::map< std::string, SeqMap1 > SeqMap ;

		Genomes genomes ;
		SeqMap seq_map ;
		std::list< std::string > path ;

		static Metagenome the_metagenome ;

	public:

		Metagenome( const char* p ) ;
		~Metagenome() { for( Genomes::iterator i = genomes.begin() ; i != genomes.end() ; ++i ) delete i->second.second ; }

		static void add_path( const std::string& s ) { the_metagenome.path.push_front( s ) ; }

		//! \brief finds a named sequence
		//! If a genome is given, only genome files whose name starts with
		//! the genome name are considered.  If genome is empty, all files
		//! are searched.
		static CompactGenome &find_sequence( const std::string& genome, const std::string& seq, Persistence ) ;

		static CompactGenome &find_genome( const std::string& genome, Persistence ) ;

		static glob_t glob_path( const std::string& genome ) ;

		static bool translate_to_genome_coords( DnaP pos, uint32_t &xpos, const config::Sequence** s_out = 0, const config::Genome** g_out = 0 ) ;

		//! \brief replacement for mmap() that knows how to free memory
		//! This is a wrapper wround mmap().  If the system mmap()
		//! fails, it will try to free up some memory by forgetting
		//! about an ephemeral genome and then call mmap again.
		//! Parameters are passed to mmap() unchanged.  If
		//! Metagenome::nommap is set, doesn't mmap() the file
		//! descriptor, but mmap()s /dev/zero instead and read()s the
		//! requested region from the file descriptor (intended for file
		//! systems where mmap() is agonizingly slow, e.g.  GCFS).
		static void *mmap( void *start, size_t length, int prot, int flags, int fd, off_t offset ) ;
} ;


#endif

