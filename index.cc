#include "index.h"
#include "util.h"

#include <cmath>

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>

CompactGenome::CompactGenome( const metaindex::Genome &g, int adv )
	: base_(), length_(0), fd_(-1)
{
	try
	{
		const char *fp = g.filename().c_str() ;
		fd_ = open( fp, O_RDONLY ) ;
		throw_errno_if_minus1( fd_, "opening", fp ) ;
		struct stat the_stat ;
		throw_errno_if_minus1( fstat( fd_, &the_stat ), "statting", fp ) ;
		length_ = the_stat.st_size ;
		void *p = mmap( 0, length_, PROT_READ, MAP_SHARED, fd_, 0 ) ;
		throw_errno_if_minus1( p, "mmapping", fp ) ;
		if( adv ) madvise( p, length_, adv ) ;
		base_.assign( (uint8_t*)p ) ;

		if( *((uint32_t const*)p) != signature ) 
			throw fp + std::string(" does not have 'DNA0' signature") ;
	}
	catch(...) 
	{
		if( base_ ) munmap( (void*)base_.unsafe_ptr(), length_ ) ;
		if( fd_ != -1 ) close( fd_ ) ;
		throw ;
	}
}

CompactGenome::~CompactGenome()
{
	if( base_ ) munmap( (void*)base_.unsafe_ptr(), length_ ) ;
	if( fd_ != -1 ) close( fd_ ) ;
}


void CompactGenome::report( uint32_t o, uint32_t l, const char* msg ) {
	if( msg )
		std::clog << '\r' << msg << ": at offset " << o << " (of " << 2*l
			<< ", " << round((50.0*o)/l) << "%)\e[K" << std::flush ;
}

FixedIndex::FixedIndex( const char* fp, unsigned w )
	: base(0), secondary(0), first_level_len(0), length(0), fd(-1), wordsize(w)
{
	try 
	{
		fd = open( fp, O_RDONLY ) ;
		throw_errno_if_minus1( fd, "opening", fp ) ;
		struct stat the_stat ;
		throw_errno_if_minus1( fstat( fd, &the_stat ), "statting", fp ) ;
		length = the_stat.st_size ;
		void *p = mmap( 0, length, PROT_READ, MAP_SHARED, fd, 0 ) ;
		throw_errno_if_minus1( p, "mmapping", fp ) ;
		base = (uint32_t*)p ;

		// std::cerr << "index base: " << p << std::endl ;

		if( base[0] != signature ) 
			throw fp + std::string(" does not have 'IDX0' signature") ;

		first_level_len = base[1] ;
		base += 2 ;
		secondary = base + first_level_len + 1 ;
	}
	catch(...)
	{
		if( base ) munmap( base-2, length ) ;
		if( fd != -1 ) close( fd ) ;
		throw ;
	}
}

FixedIndex::~FixedIndex() 
{
	if( base ) munmap( base-2, length ) ;
	if( fd != -1 ) close( fd ) ;
}

unsigned FixedIndex::lookup1( Oligo o, std::vector<Seed>& v, uint32_t cutoff, int32_t offs ) const 
{
	assert( o < first_level_len ) ;

	Seed seed ;
	seed.size = wordsize ;
	seed.offset = offs ;

	if( base[o+1] - base[o] >= cutoff ) return 0 ;

	for( uint32_t p = base[o] ; p != base[o+1] ; ++p )
	{
		seed.diagonal = secondary[p] - (uint32_t)offs ;
		v.push_back( seed ) ;
	}
	return base[o+1] - base[o] ;
} 

unsigned FixedIndex::lookup( const std::string& dna, std::vector<Seed>& v, uint32_t cutoff ) const
{
	Oligo o_f = 0, o_r = 0 ;
	Oligo mask = ~( ~0 << (wordsize * 2) ) ;

	int s = 2 * wordsize - 2 ;
	unsigned filled = 0 ;
	unsigned total = 0 ;
	int32_t fraglen = dna.size() ;

	for( int32_t offset = 0 ; offset != fraglen ; ++offset )
	{
		o_f >>= 2 ;
		o_r = (o_r << 2) & mask ;
		++filled ;

		switch( dna[offset] )
		{
			case 'a': case 'A': o_f |= 0 << s ; o_r |= 2 ; break ;
			case 'c': case 'C': o_f |= 1 << s ; o_r |= 3 ; break ;
			case 'u': case 'U':
			case 't': case 'T': o_f |= 2 << s ; o_r |= 0 ; break ;
			case 'g': case 'G': o_f |= 3 << s ; o_r |= 1 ; break ;
			default: filled = 0 ; break ;
		}
		if( filled >= wordsize ) 
			total += lookup1( o_f, v, cutoff,   offset - wordsize + 1 ) 
				   + lookup1( o_r, v, cutoff, - offset                ) ;
	}
	combine_seeds( v ) ;
	return total ;
}

