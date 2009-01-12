#include "index.h"
#include "util.h"

#include <cmath>

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>

using namespace config ;
using namespace std ; 

static inline int open_( std::string fp, const char* ext )
{
	int fd = open( fp.c_str(), O_RDONLY ) ;
	if( fd == -1 && errno != ENOENT )
		throw_errno_if_minus1( fd, "opening", fp.c_str() ) ;
	if( fd != -1 ) return fd ;

	fp.append( ext ) ;
	fd = open( fp.c_str(), O_RDONLY ) ;
	if( fd == -1 && errno != ENOENT )
		throw_errno_if_minus1( fd, "opening", fp.c_str() ) ;
	return fd ;
}

CompactGenome::CompactGenome( const std::string &name, const config::Config &c, int adv )
	: base_(), file_size_(0), length_(0), fd_(-1), contig_map_(), g_()
{
	try
	{
		std::string fp = name ;
		fd_ = open_( fp, ".dna" ) ;
		for( int i = 0 ; fd_ == -1 && i != c.genome_path_size() ; ++i )
		{
			fp = c.genome_path(i) + '/' + name ;
			fd_ = open_( fp, ".dna" ) ;
		}
		throw_errno_if_minus1( fd_, "opening", name.c_str() ) ;

		struct stat the_stat ;
		throw_errno_if_minus1( fstat( fd_, &the_stat ), "statting", fp.c_str() ) ;
		file_size_ = the_stat.st_size ;
		void *p = mmap( 0, file_size_, PROT_READ, MAP_SHARED, fd_, 0 ) ;
		throw_errno_if_minus1( p, "mmapping", fp.c_str() ) ;
		base_.assign( (uint8_t*)p ) ;

		if( ((uint32_t const*)p)[0] != signature ) 
			throw fp + string(" does not have 'DNA1' signature") ;

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
				contig_map_[ c.offset() ].first = s ;
				contig_map_[ c.offset() ].second = c ;
			}
		}
		contig_map_[ g_.total_size() ] = make_pair( Sequence(), Contig() ) ;
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

FixedIndex::FixedIndex( const std::string& name, const config::Config& c )
	: base(0), secondary(0), first_level_len(0), length(0), fd_(-1)
{
	try 
	{
		std::string fp = name ;
		fd_ = open_( fp, ".idx" ) ;

		for( int i = 0 ; fd_ == -1 && i != c.index_path_size() ; ++i )
		{
			fp = c.index_path(i) + '/' + name ;
			fd_ = open_( fp, ".idx" ) ;
		}
		throw_errno_if_minus1( fd_, "opening", name.c_str() ) ;

		struct stat the_stat ;
		throw_errno_if_minus1( fstat( fd_, &the_stat ), "statting", fp.c_str() ) ;
		length = the_stat.st_size ;
		p_ = mmap( 0, length, PROT_READ, MAP_SHARED, fd_, 0 ) ;
		throw_errno_if_minus1( p_, "mmapping", fp.c_str() ) ;

		if( *(const uint32_t*)p_ != signature ) 
			throw fp + string(" does not have 'IDX1' signature") ;

		uint32_t meta_len = ((const uint32_t*)p_)[1] ;
		if( !ci_.ParseFromArray( (const char*)p_ + 8, meta_len ) )
			throw "error parsing meta information" ;

		first_level_len = 1 << (2 * ci_.wordsize()) + 1 ;
		base = (const uint32_t*)( (const char*)p_ + 8 + ((3+meta_len) & ~3) ) ;
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
//! \return number of seeds found, including repetitive ones
unsigned FixedIndex::lookup1( Oligo o, vector<Seed>& v, uint32_t cutoff, int32_t offs ) const 
{
	assert( o < first_level_len ) ;

	Seed seed ;
	seed.size = ci_.wordsize() ;
	seed.offset = offs ;

	if( base[o+1] - base[o] < cutoff )
	{
		for( uint32_t p = base[o] ; p != base[o+1] ; ++p )
		{
			seed.diagonal = secondary[p] - (uint32_t)offs ;
			v.push_back( seed ) ;
		}
	}
	return base[o+1] - base[o] ;
} 


//! \brief looks up a whole sequence
//! The sequence is split into words as
//! appropriate for the index, then each one of them is looked
//! up.  This method can be implemented for any kind of index, whether
//! based on fixed words or not
//!
//! \param dna sequence to search
//! \param v receiver for resulting seeds
//! \param cutoff disregard oligos more frequent than this
//! \return number of seeds found
//!
//! \todo move cutoff parameter somewhere else to improve modularity

unsigned FixedIndex::lookup( const QSequence& dna, std::vector<Seed>& v, uint32_t cutoff ) const
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

		switch( dna[offset] )
		{
			case 1: o_f |= 0 << s ; o_r |= 2 ; break ;
			case 2: o_f |= 1 << s ; o_r |= 3 ; break ;
			case 4: o_f |= 2 << s ; o_r |= 0 ; break ;
			case 8: o_f |= 3 << s ; o_r |= 1 ; break ;
			default: filled = 0 ; break ;
		}
		if( filled >= ci_.wordsize() ) 
			total += lookup1( o_f, v, cutoff,   offset - ci_.wordsize() + 1 ) 
				   + lookup1( o_r, v, cutoff, - offset                      ) ;
	}
	combine_seeds( v ) ;
	return total ;
}

bool CompactGenome::translate_back( DnaP pos, std::string& sequ_id, uint32_t& offset ) const 
{
	ContigMap::const_iterator high = contig_map_.upper_bound( pos.abs() - base_ ) ;
	if( pos.abs() < base_ ) return false ; // before start
	if( high == contig_map_.end() ) return false ; // after end
	--high ;

	const config::Sequence &sequ = high->second.first ;
	const config::Contig &contig = high->second.second ;

	sequ_id = sequ.name() ;
	offset = pos.abs() - base_ - contig.offset() + contig.range_start() ;
	return true ;
}

