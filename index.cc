#include "index.h"
#include "util.h"

#include <cmath>

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>

using namespace config ;
using namespace std ; 

CompactGenome::CompactGenome( const std::string &name, const config::Config &c, int adv )
	: base_(), file_size_(0), length_(0), fd_(-1), contig_map_(), g_()
{
	try
	{
		fd_ = path_open( name, "dna", "ANFO_PATH", c.genome_path().begin(), c.genome_path().end() ) ;

		struct stat the_stat ;
		throw_errno_if_minus1( fstat( fd_, &the_stat ), "statting", name.c_str() ) ;
		file_size_ = the_stat.st_size ;
		void *p = mmap( 0, file_size_, PROT_READ, MAP_SHARED, fd_, 0 ) ;
		throw_errno_if_minus1( p, "mmapping", name.c_str() ) ;
		base_.assign( (uint8_t*)p ) ;

		if( ((uint32_t const*)p)[0] != signature ) 
			throw name + string(" does not have 'DNA1' signature") ;

		uint32_t meta_off = ((uint32_t const*)p)[1] ;
		uint32_t meta_len = ((uint32_t const*)p)[2] ;

		if( !g_.ParseFromArray( (const char*)p + meta_off, meta_len ) )
			throw "error parsing meta information" ;

		length_  = meta_off ;
		if( adv ) madvise( p, length_, adv ) ;

		for( int i = 0 ; i != g_.sequence_size() ; ++i )
		{
			const Sequence &s = g_.sequence(i) ;
			for( int j = 0 ; j != s.contig_size() ; ++j )
			{
				const Contig &c = s.contig(j) ;
				contig_map_[ c.offset() ].first = i ;
				contig_map_[ c.offset() ].second = j ;
			}
		}
		contig_map_[ g_.total_size() ] = make_pair( -1, -1 ) ;
	}
	catch(...) 
	{
		if( base_ ) munmap( (void*)base_.unsafe_ptr(), file_size_ ) ;
		if( fd_ != -1 ) close( fd_ ) ;
		throw ;
	}
}

CompactGenome::~CompactGenome()
{
	if( base_ ) munmap( (void*)base_.unsafe_ptr(), file_size_ ) ;
	if( fd_ != -1 ) close( fd_ ) ;
}


void CompactGenome::report( uint32_t o, uint32_t l, const char* msg ) {
	if( msg )
		clog << '\r' << msg << ": at offset " << o << " (of " << 2*l
			<< ", " << round((50.0*o)/l) << "%)\e[K" << flush ;
}

FixedIndex::FixedIndex( const std::string& name, const config::Config& c, int adv )
	: base(0), secondary(0), first_level_len(0), length(0), fd_(-1)
{
	void *p = 0 ;
	try 
	{
		fd_ = path_open( name, "idx", "ANFO_PATH", c.genome_path().begin(), c.genome_path().end() ) ;

		struct stat the_stat ;
		throw_errno_if_minus1( fstat( fd_, &the_stat ), "statting", name.c_str() ) ;
		length = the_stat.st_size ;
		p = mmap( 0, length, PROT_READ, MAP_SHARED, fd_, 0 ) ;
		throw_errno_if_minus1( p, "mmapping", name.c_str() ) ;
		p_ = p ;

		if( adv ) madvise( p, length, adv ) ;

		if( *(const uint32_t*)p != signature ) 
			throw name + string(" does not have 'IDX1' signature") ;

		uint32_t meta_len = ((const uint32_t*)p)[1] ;
		if( !ci_.ParseFromArray( (const char*)p + 8, meta_len ) )
			throw "error parsing meta information" ;

		// XXX: this line is actually correct, but indexes need to be
		// rebuilt
		// first_level_len = (1 << (2 * ci_.wordsize())) + 1 ;
		first_level_len = 1 << ((2 * ci_.wordsize()) + 1) ;
		base = (const uint32_t*)( (const char*)p + 8 + ((3+meta_len) & ~3) ) ;
		secondary = base + first_level_len ; 
	}
	catch(...)
	{
		if( p_ ) munmap( (void*)p_, length ) ;
		if( fd_ != -1 ) close( fd_ ) ;
		throw ;
	}
}

FixedIndex::~FixedIndex() 
{
	if( p_ ) munmap( (void*)p_, length ) ;
	if( fd_ != -1 ) close( fd_ ) ;
}

//! \brief directly looks up an oligo
//! A single oligo is looked up and results are appended to a vector.
//! Seeds that occur too often can be ignored, but for statistical
//! purposes, they are always counted in the result.
//!
//! \param o oligo to be looked up.
//! \param v receiver for the resulting seeds
//! \param cutoff discard oligos more frequent than this
//! \param offs offset value to be placed in the seeds
//! \param num_useless filled in with the number of seeds ignored due to
//!                    being too frequent
//! \return number of seeds found, including repetitive ones
unsigned FixedIndex::lookup1( Oligo o, vector<Seed>& v, uint32_t cutoff, int32_t offs, int *num_useless ) const 
{
	assert( o < first_level_len ) ;

	Seed seed ;
	seed.size = ci_.wordsize() ;
	seed.offset = offs ;

	if( base[o+1] - base[o] < cutoff )
	{
		for( uint32_t p = base[o] ; p != base[o+1] ; ++p )
		{
			seed.diagonal = (int64_t)secondary[p] - (int64_t)offs ;
			v.push_back( seed ) ;
		}
	}
	else if( num_useless ) *num_useless += base[o+1] - base[o] ;
	return base[o+1] - base[o] ;
} 

//! \brief looks up an oligo with up to one mismatch
//! The oligo itself will first be looked for, then all its
//! one-substituion variants are generated (by successively xor'ing
//! every position with 01, 10 and 11) and looked up, too.
//!
//! \todo The way this is implemented, the sorting of seeds is
//!       completely overwhelmed, negating the performance gain from
//!       seeding at all.  Needs to be fixed by a better data structure.
//!
//! \param o oligo to be looked up.
//! \param v receiver for the resulting seeds
//! \param cutoff discard oligos more frequent than this
//! \param offs offset value to be placed in the seeds
//! \param num_useless filled in with the number of seeds ignored due to
//!                    being too frequent
//! \return number of seeds found, including repetitive ones

unsigned FixedIndex::lookup1m( Oligo o, vector<Seed>& v, uint32_t cutoff, int32_t offs, int *num_useless ) const 
{
	unsigned total = lookup1( o, v, cutoff, offs, num_useless ) ;
	Oligo m1 = 1, m2 = 2, m3 = 3 ;
	for( size_t i = 0 ; i != ci_.wordsize() ; ++i )
	{
		total += lookup1( o ^ m1, v, cutoff, offs, num_useless ) ;
		total += lookup1( o ^ m2, v, cutoff, offs, num_useless ) ;
		total += lookup1( o ^ m3, v, cutoff, offs, num_useless ) ;
		m1 <<= 2 ;
		m2 <<= 2 ;
		m3 <<= 2 ;
	}
	return total ;
} 

//! \brief looks up a whole sequence
//! The sequence is split into words as appropriate for the index, then
//! each one of them is looked up.  This method can be implemented for
//! any kind of index, whether based on fixed words or not.
//!
//! \param dna sequence to search
//! \param v receiver for resulting seeds
//! \param near_perfect set to one to search for seeds with one mismatch
//! \param num_useless is filled in with the number of seeds thrown away
//!                    due to being too frequent
//! \param cutoff disregard oligos more frequent than this
//! \return number of seeds found

unsigned FixedIndex::lookupS( const QSequence& dna, std::vector<Seed>& v,
		bool near_perfect, int *num_useless, uint32_t cutoff ) const
{
	Oligo o_f = 0, o_r = 0 ;
	Oligo mask = ~( ~0 << (ci_.wordsize() * 2) ) ;

	int s = 2 * ci_.wordsize() - 2 ;
	unsigned filled = 0 ;
	unsigned total = 0 ;
	int32_t fraglen = dna.length() ;

	for( int32_t offset = 0 ; offset != fraglen ; ++offset )
	{
		o_f >>= 2 ;
		o_r = (o_r << 2) & mask ;
		++filled ;

		switch( dna[offset].ambicode )
		{
			case 1: o_f |= 0 << s ; o_r |= 2 ; break ;
			case 2: o_f |= 1 << s ; o_r |= 3 ; break ;
			case 4: o_f |= 2 << s ; o_r |= 0 ; break ;
			case 8: o_f |= 3 << s ; o_r |= 1 ; break ;
			default: filled = 0 ; break ;
		}
		if( filled >= ci_.wordsize() ) 
		{
			if( near_perfect )
				total += lookup1m( o_f, v, cutoff,   offset - ci_.wordsize() + 1, num_useless ) 
					   + lookup1m( o_r, v, cutoff, - offset                     , num_useless ) ;
			else
				total += lookup1( o_f, v, cutoff,   offset - ci_.wordsize() + 1, num_useless ) 
					   + lookup1( o_r, v, cutoff, - offset                     , num_useless ) ;
		}
	}
	combine_seeds( v ) ;
	return total ;
}

const config::Sequence *CompactGenome::translate_back( DnaP pos, uint32_t& offset ) const 
{
	ContigMap::const_iterator high = contig_map_.upper_bound( pos.abs() - base_ ) ;
	if( pos.abs() < base_ ) return 0 ; // before start
	if( high == contig_map_.end() ) return 0 ; // after end
	--high ;

	const config::Sequence &sequ = g_.sequence( high->second.first ) ;
	const config::Contig &contig = sequ.contig( high->second.second ) ;

	offset = pos.abs() - base_ - contig.offset() + contig.range_start() ;
	return &sequ ;
}

DnaP CompactGenome::find_pos( const std::string& seq, uint32_t pos ) const
{
	for( int i = 0 ; i != g_.sequence_size() ; ++i )
		if( g_.sequence(i).name() == seq )
			for( int j = 0 ; j != g_.sequence(i).contig_size() ; ++j )
				if( g_.sequence(i).contig(j).range_start() <= pos &&
				    g_.sequence(i).contig(j).range_end() >= pos )
					return base_ + g_.sequence(i).contig(j).offset() + (pos - g_.sequence(i).contig(j).range_start()) ;
	return DnaP(0) ;
}


