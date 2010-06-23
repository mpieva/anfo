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
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <glob.h>
#include <sys/mman.h>

//! \brief a genome as stored in a DNA file
//! This class mmaps a DNA file and wraps it with a sensible interface.
class CompactGenome
{
	public:
		// Why are we using indices here instead of pointers?  Because
		// pointers would be internal and make copying of this object a
		// lot harder.

		// key is offset into genome, value is (index of sequence, index
		// of contig)
		typedef std::map< uint32_t, std::pair< int, int > > ContigMap ;
		
		// value is index of contig
		typedef std::map< uint32_t, int > PosnMap1 ;

		// first half of value is index of sequence
		typedef std::map< std::string, std::pair< int, PosnMap1 > > PosnMap ;

		static void cleanup( const CompactGenome* ) {}

	private:
		DnaP base_ ;
		size_t file_size_ ;
		uint32_t length_ ;
		int fd_ ;
		ContigMap contig_map_ ;
		PosnMap posn_map_ ;

		CompactGenome( const CompactGenome& ) ; // not implemented
		void operator = ( const CompactGenome& ) ; // not implemented

		~CompactGenome() ; // must control life cycle
		friend class Metagenome ;

	public:
		config::Genome g_ ;
		mutable int refcount_ ;

	public:
		//! \brief makes accessible a genome file
		//! \param name file name of the genome
		//! \param c program configuration, needed for the search path
		CompactGenome( const std::string& name ) ;

		void add_ref() const { ++refcount_ ; }

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
		//! \param mk_gap functor that's called at every gap
		//! \param slice scan this slice of a big genome
		//! \param slices total number of slices
		//! \param msg if set, switches on progress reports and is
		//!            included in them
		template< typename F, typename G > void scan_words( unsigned w, F mk_word, G mk_gap, int slice = 0, int slices = 1, const char* msg = 0 ) const ;

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
// struct Seed
// {
	// int64_t diagonal ;
	// uint32_t size ;
	// int32_t offset ;
// } ;
struct Matches
{
	const uint32_t *begin, *end ;
	int32_t offs ;
	uint16_t wordsize ;
	uint16_t stride ;
	int64_t diag ;
} ;


// match lists are sorted backwards by reference coordinate and
// therefore backwards by diagonal; we also choose to sort
// backwards on the offset.
struct compare_match_lists {
	bool operator()( const Matches &a, const Matches &b ) {
		if( b.diag < a.diag ) return true ;
		if( a.diag < b.diag ) return false ;
		return b.offs < a.offs ;
	}
} ;

struct compare_matches_for_heap {
	// note the `wrong' order; this is intended for a heap!
	bool operator()( const Matches &a, const Matches &b ) {
		return compare_match_lists()( b, a ) ;
	}
} ;


//! \brief Collection of short matches.
//! We actually collect the ranges of sorted(!) matches in the index
//! data structure, then merge-sort them on demand.
class PreSeeds
{
	private:
		struct equal_match_lists {
			bool operator()( const Matches &a, const Matches &b ) {
				return a.begin == b.begin && a.offs == b.offs ;
			}
		} ;
		std::vector< Matches > heap_ ;

		bool normalize() 
		{
			while( heap_.back().begin != heap_.back().end &&
					*heap_.back().begin % heap_.back().stride != 0 ) 
				++heap_.back().begin ;

			if( heap_.back().begin == heap_.back().end ) {
				heap_.pop_back() ;
				return false ;
			}
			else {
				heap_.back().diag = (int64_t)(*heap_.back().begin) - (int64_t)(heap_.back().offs) ;
				return true ;
			}
		}

	public:
		PreSeeds() {}

		void post( const uint32_t *begin, const uint32_t *end, int32_t offs, int wordsize, int stride )
		{
			if( begin != end ) {
				heap_.push_back( Matches() ) ;
				heap_.back().begin = begin ;
				heap_.back().end = end ;
				heap_.back().offs = offs ;
				heap_.back().wordsize = wordsize ;
				heap_.back().stride = stride ;
				normalize() ;
			}
		}

		bool empty() const { return heap_.empty() ; }

		const Matches& get() const { return heap_.front() ; }

		void start_traversal() {
			std::sort( heap_.begin(), heap_.end(), compare_match_lists() ) ;
			heap_.erase( std::unique( heap_.begin(), heap_.end(), equal_match_lists() ), heap_.end() ) ;
			std::make_heap( heap_.begin(), heap_.end(), compare_matches_for_heap() ) ;
		}

		void take() 
		{
			std::pop_heap( heap_.begin(), heap_.end(), compare_matches_for_heap() ) ;

			++heap_.back().begin ;
			assert( heap_.back().begin == heap_.back().end ||
				heap_.back().begin[-1] > heap_.back().begin[0] ) ;

			if( normalize() ) std::push_heap( heap_.begin(), heap_.end(), compare_matches_for_heap() ) ;
		}
} ;


// inline std::ostream& operator << ( std::ostream& o, const Seed& s )
// {
	// return o << '@' << s.offset << '+' << s.diagonal << ':' << s.size ;
// }


class FixedIndex 
{
	public:
		struct LookupParams {
			uint32_t cutoff ; 				// more common seeds will be ignored
			uint32_t allow_mismatches ;		// number of changed positions per seed word
			uint32_t wordsize ;				// nucleotides per seed word
			uint32_t stride ;				// tiling stepsize
		} ;

		enum { signature = 0x33584449u } ; // "IDX3"

		FixedIndex() : p_(0), base(0), secondary(0), first_level_len(0), length(0), fd_(0) {}

		//! \brief loads an index from a file
		//! \param name filename
		FixedIndex( const std::string &name ) ;
		~FixedIndex() { if( p_ ) munmap( (void*)p_, length ) ; if( fd_ != -1 ) close( fd_ ) ; }

		unsigned lookupS( const std::string&  seq, PreSeeds&, const LookupParams &p, int *num_useless ) const ;
		unsigned lookup1(  Oligo, PreSeeds&, const LookupParams &p, int32_t offs, int *num_useless ) const ;
		unsigned lookup1m( Oligo, PreSeeds&, const LookupParams &p, int32_t offs, int *num_useless ) const ;
		unsigned lookup2m( Oligo, PreSeeds&, const LookupParams &p, int32_t offs, int *num_useless ) const ;

		operator const void * () const { return base ; }
		const config::CompactIndex& metadata() const { return meta_ ; }

		void swap( FixedIndex& i ) {
			std::swap( p_, i.p_ ) ;
			std::swap( base, i.base ) ;
			std::swap( secondary, i.secondary ) ;
			std::swap( first_level_len, i.first_level_len ) ;
			std::swap( length, i.length ) ;
			std::swap( fd_, i.fd_ ) ;
		}

	private:
		const void* p_ ;
		const uint32_t *base, *secondary ;
		uint32_t first_level_len ;
		uint64_t length ;
		int fd_ ;
		config::CompactIndex meta_ ;

		unsigned lookup1m( Oligo, PreSeeds&, const LookupParams &p, int32_t offs, int *num_useless, size_t ) const ;
} ;

template< typename F, typename G > void CompactGenome::scan_words(
		unsigned w, F mk_word, G mk_gap, int slice, int slices, const char* msg ) const
{
	if( (unsigned)std::numeric_limits< Oligo >::digits < 4 * w ) 
		throw "cannot build index: oligo doesn't fit" ;

	// start here, in case we want slices
	uint32_t offs = 2 * slice * (int64_t)length_ / slices ;
	uint32_t eoffs = 2 * (slice+1) * (int64_t)length_ / slices ;
	Oligo dna = 0 ;

    // do not start before first contig (that's the header region)
    assert( g_.sequence_size() && g_.sequence(0).contig_size() ) ;
    offs = std::max( offs, g_.sequence(0).contig(0).offset()-1 ) ;
	while( base_[ offs ] != 0 ) ++offs ;	// find first gap

	for( unsigned i = 0 ; i != w ; ++i )	// fill first word
	{
		dna <<= 4 ;
		dna |= base_[ offs ] ;
		++offs ;
	}

	report(offs,length_,msg) ;
	// run to end of slice, but stop at gaps only
	while( offs < eoffs || base_[ offs-1 ] )
	{
		if( (offs & 0xffffff) == 0 ) report(offs,length_,msg) ;

		dna <<= 4 ;
		dna |= base_[ offs ] ;
		++offs ;

		if( !base_[offs] ) mk_gap( offs ) ;

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
/*
struct compare_diag_then_offset {
	bool operator()( const Seed& a, const Seed& b ) {
		if( a.diagonal > b.diagonal ) return true ;
		if( b.diagonal > a.diagonal ) return false ;
		return a.offset > b.offset ;
	}
} ;
*/

//! \brief stores a new seed
//! We immediately combine adjacent, overlapping and adjacent seeds: if
//! seeds have the same diagonal and overlap or are adjacent, they can
//! always be combined.
//inline void add_seed( PreSeeds& v, int64_t diag, int32_t offs, uint32_t size )
//{
	//v.push_back( Seed() ) ;
	//v.back().diagonal = diag ;
	//v.back().offset = offs ;
	//v.back().size = size ;
//}

//! \brief Combines short, adjacent seeds into longer ones.
//! This is the cheap method to combine seeds: only overlapping and
//! adjacent seeds are combined, neighboring diagonals are not
//! considered.  The code is short and direct, and works well even for
//! imperfect seeds.
//!
//! How to do this?  We sort seeds first by diagonal index, then by
//! offset.  Seeds are adjacent iff they have the same diagonal index
//! and their offsets differ by no more than the seed size.
//!
//! \note Formerly we tried to somehow deal with neighboring diagonals.
//!       This has been declared as not worth the hassle, so it was
//!       dropped.  If we get seeds for the same region on neighboring
//!       diagonals, they are useless repeats anyway and excluded from
//!       the index to begin with.
//!
//! \param v container of seeds, will be modifed in place.
//! \param m minimum length of a good seed
//! \param ss output container for seed positions
//! \return number of good seeds produced
inline int combine_seeds( PreSeeds& v, uint32_t m, output::Seeds *ss )
{
	int out = 0 ;
	if( !v.empty() )
	{
		v.start_traversal() ;
		// std::sort( v.begin(), v.end(), compare_diag_then_offset() ) ;

		// combine overlapping and adjacent seeds into larger ones
		// PreSeeds::const_iterator a = v.begin(), e = v.end() ;
		Matches s = v.get() ;
		v.take() ;
		for( ; !v.empty() ; v.take() )
		{
			const Matches& a = v.get() ;
			assert( a.diag <= s.diag ) ;
			if( a.diag == s.diag ) assert( a.offs <= s.offs ) ;

			if( a.diag == s.diag &&
					(a.offs >= 0) == (s.offs >= 0) &&
					a.offs + (int32_t)a.wordsize >= s.offs )
			{
				s.offs = a.offs ;
				s.wordsize += s.offs - a.offs ;
			}
			else
			{
				if( s.wordsize >= m ) 
				{
					ss->mutable_ref_positions()->Add( s.diag + s.offs ) ;
					ss->mutable_query_positions()->Add( s.offs ) ;
					ss->mutable_seed_sizes()->Add( s.wordsize ) ;
					++out ;
				}
				s = a ;
			}
		}
		if( s.wordsize >= m ) 
		{
			ss->mutable_ref_positions()->Add( s.diag + s.offs ) ;
			ss->mutable_query_positions()->Add( s.offs ) ;
			ss->mutable_seed_sizes()->Add( s.wordsize ) ;
			++out ;
		}
	}
	return out ;
}

typedef std::map< std::string, CompactGenome > Genomes ;
typedef std::map< std::string, FixedIndex > Indices ;
typedef Holder< const CompactGenome > GenomeHolder ;

class Metagenome
{
	public:
		static int nommap ;
		static void make_room() ;

	private:
		// maps file name to genome object
		typedef std::map< std::string, CompactGenome* > Genomes ;

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
		~Metagenome() { for( Genomes::iterator i = genomes.begin() ; i != genomes.end() ; ++i ) delete i->second ; }

		static void add_path( const std::string& s ) { the_metagenome.path.push_front( s ) ; }

		//! \brief finds a named sequence
		//! If a genome is given, only genome files whose name starts with
		//! the genome name are considered.  If genome is empty, all files
		//! are searched.
		static GenomeHolder find_sequence( const std::string& genome, const std::string& seq ) ;

		static GenomeHolder find_genome( const std::string& genome ) ;

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
		static void *mmap( void *start, size_t length, int *fd, off_t offset ) ;
} ;

#endif

