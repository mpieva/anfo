#include "Index.h"
#include "util.h"

#include <cmath>

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>

CompactGenome::CompactGenome( const char *fp )
	: base(), length(0), fd(-1)
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
		base.assign( p ) ;

		if( ((uint32_t const*)(void const*)base)[0] != signature ) 
			throw fp + std::string(" does not have 'DNA0' signature") ;
	}
	catch(...) 
	{
		if( base ) munmap( base.unsafe_ptr(), length ) ;
		if( fd != -1 ) close( fd ) ;
		throw ;
	}
}

CompactGenome::~CompactGenome()
{
	if( base ) munmap( base.unsafe_ptr(), length ) ;
	if( fd != -1 ) close( fd ) ;
}


void CompactGenome::report( uint32_t o, uint32_t l, const char* msg ) {
	if( msg )
		std::clog << '\r' << msg << ": at offset " << o << " (of " << 2*l
			<< ", " << round((50.0*o)/l) << "%)\e[K" << std::flush ;
}

FixedIndex::FixedIndex( const char* fp, unsigned w, unsigned c ) 
	: base(0), secondary(0), first_level_len(0), length(0), fd(-1), wordsize(w), cutoff(c)
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

		std::cerr << "index base: " << p << std::endl ;

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

unsigned FixedIndex::lookup1( Oligo o, std::vector<Seed>& v, int32_t offs ) const 
{
	assert( o < first_level_len ) ;

	Seed seed ;
	seed.size = wordsize ;
	seed.offset = offs ;

	if( base[o+1] - base[o] >= cutoff ) return 0 ;

	for( uint32_t p = base[o] ; p != base[o+1] ; ++p )
	{
		seed.diagonal = secondary[p] - offs ;
		v.push_back( seed ) ;
	}
	return base[o+1] - base[o] ;
} 

unsigned FixedIndex::lookup( const std::string& dna, std::vector<Seed>& v ) const
{
	Oligo o_f = 0, o_r = 0 ;
	int s = 2 * wordsize - 2 ;
	unsigned filled = 0 ;
	unsigned total = 0 ;
	int32_t fraglen = dna.size() ;

	for( int32_t offset = 0 ; offset != fraglen ; ++offset )
	{
		o_f >>= 2 ;
		o_r = (o_r << 2) & ((1<<s)-1) ;
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
			total += lookup1( o_f, v, offset - wordsize + 1 ) 
				   + lookup1( o_r, v, offset - wordsize + 1 - fraglen ) ;
	}
	return total ;
}

