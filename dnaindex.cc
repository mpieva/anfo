#include "util.h"

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>

/*
 * Scans a dna file, creating an index of fixed-size words. 
 * We need a two-level index and want to keep it simple.  To this end,
 * we'll store the configuration of the index externally in text format
 * and the index files themselves are little more than arrays.  So, this
 * is what an index file looks like:
 *
 * - The signature "IDX0" in ASCII.
 * - The size of the first level array in words.
 * - The first level array itself.  Every DNA-word is converted to an
 *   index by assigning 0,1,2,3 to A,C,T,G and writing nucleotides LSB
 *   to MSB.  The entry at that position is an offset into the second
 *   level array.
 * - The size of the second level array in words.
 * - The second level array.  This is a dense array of offsets into a
 *   compact genome.  Entries in the first level array point to the
 *   begin of the subarray of indices where their key can be found.
 *
 * All words are 32bit unsigned integers written LSB first.  No such
 * word can overflow, since we limit ourselves to genomes no larger than
 * 4G bases.  The index file can be larger than 4G (and often will be).
 * If so, they are only useful on 64bit platforms.
 */

/*
 * This code probably makes little to no sense on a 32bit machine.
 * Anyway, we take care that a 32bit version will still work, though
 * limiting the index size to comfortably fit into the address space
 * will make it impossible to index more than roughly a quarter billion
 * bases.
 *
 * We use (unsigned long long int) to store words cut out from the
 * genome.  Assuming it has at least 64 bits of storage, this limits the
 * word size to 16 bases.  (No problem, more than that is probably
 * impractical anyway.)  
 */

/* Note to self: creating and mmapping a sparse file, then filling it in
 * random order is an exceptionally bad idea...
 */

typedef long long unsigned int Oligo ;
typedef unsigned char uint8_t ;			// ???

class CompactGenome
{
	private:
		uint8_t *base ;
		unsigned length ;
		int fd ;

	public:
		CompactGenome( const char* fp ) ;
		~CompactGenome() ;

		// get nucleotide at position ix, relative to the beginning of
		// the file
		uint8_t operator[]( unsigned ix ) {
			uint8_t w = base[ix >> 1] ;
			if( ix & 1 ) w >>= 4 ;
			return w & 0xf ;
		}

		// scan over finite words of the dna; these may well contain
		// ambiguity codes, but are guaranteed not to contain gaps
		template< typename F, typename G >
			void scan_words( int w, F mk_word, G consume_word, const char* msg = 0 ) ;

		uint8_t deref( unsigned offset ) const
		{ return (base[offset >> 1] >> (4 * (offset & 1))) & 0xf ; }

} ;

CompactGenome::CompactGenome( const char *fp )
	: base(0), length(0), fd(-1)
{
	int fd = open( fp, O_RDONLY ) ;
	throw_errno_if_minus1( fd, "opening", fp ) ;
	struct stat the_stat ;
	throw_errno_if_minus1( fstat( fd, &the_stat ), "statting", fp ) ;
	length = the_stat.st_size ;
	void *p = mmap( 0, length, PROT_READ, MAP_SHARED, fd, 0 ) ;
	throw_errno_if_minus1( p, "mmapping", fp ) ;
	base = (uint8_t*)p ;
}

CompactGenome::~CompactGenome()
{
	if( base ) munmap( base, length ) ;
	if( fd != -1 ) close( fd ) ;
}

// Scan over the dna words in the genome.  Size of the words is limited
// to what fits into a (long long unsigned), typically 16.
//
// Words are encodes as four bits per nucleotide, first nucleotide in
// the MSB(!).  See make_dense_word for why that makes sense.  Unused
// MSBs in words passed to mk_word contain junk.
inline void report( unsigned o, unsigned l, const char* msg ) {
	if( msg )
		std::clog << '\r' << msg << ": at offset " << o << " (of " << 2*l
			<< ", " << round((50.0*o)/l) << "%)\e[K" << std::flush ;
}

template< typename F, typename G >
void CompactGenome::scan_words( int w, F mk_word, G consume_word, const char* msg ) {
	assert( std::numeric_limits< Oligo >::digits >= 4 * w ) ;
#ifdef _BSD_SOURCE
	madvise( base, length, MADV_SEQUENTIAL ) ;
#endif
	unsigned offs = 0 ;
	Oligo dna = 0 ;

	while( deref( offs ) != 0 ) ++offs ;		// first first gap
	for( int i = 0 ; i != w ; ++i )				// fill first word
	{
		dna <<= 4 ;
		dna |= deref( offs ) ;
		++offs ;
	}

	report(offs,length,msg) ;
	while( offs != 2 * length )
	{
		if( (offs & 0xfffff) == 0 ) report(offs,length,msg) ;

		dna <<= 4 ;
		dna |= deref( offs ) ;
		++offs ;

		// throw away words containing gap symbols
		// (This is necessary since we may want to construct
		// discontiguous words, but not if a "don't care" position is a
		// gap.)
		bool clean = true ;
		for( int i = 0 ; i != w ; ++i )
			clean &= ((dna >> (4*i)) & 0xf) != 0 ;

		if( clean ) mk_word( w, offs-w, dna, consume_word ) ;
	}
	std::clog << "\r\e[K" << std::flush ;
}

template< typename G >
void make_dense_word1( int w, unsigned offs, Oligo dna, unsigned acc, G consume_word )
{
	// std::cerr << __PRETTY_FUNCTION__ << std::endl ;
	// std::cerr << "w = " << std::dec << w << ", offs = " << std::dec << offs 
		// << ", dna = " << std::hex << dna << ", acc = " << std::hex << acc << std::endl ;

	while( w ) {
		// Loop over set bits in ambiguity codes.  Special cased and
		// unrolled for speed.  Note that we need recursion, but we can
		// avoid it in the most common case of not having any ambiguity
		// code.
		int nt = dna & 0xf ;
		dna >>= 4 ;
		acc <<= 2 ;
		--w ;

		switch ( nt )
		{
			case  0: break ; // shouldn't actually occur
			case  1: acc |= 0 ; break ;
			case  2: acc |= 1 ; break ;
			case  3: make_dense_word1( w, offs, dna, acc | 0, consume_word ) ;
					 acc |= 1 ; break ;
			case  4: acc |= 2 ; break ;
			case  5: make_dense_word1( w, offs, dna, acc | 0, consume_word ) ;
					 acc |= 2 ; break ;
			case  6: make_dense_word1( w, offs, dna, acc | 1, consume_word ) ;
					 acc |= 2 ; break ;
			case  7: make_dense_word1( w, offs, dna, acc | 0, consume_word ) ;
					 make_dense_word1( w, offs, dna, acc | 1, consume_word ) ;
					 acc |= 2 ; break ;
			case  8: acc |= 3 ; break ;
			case  9: make_dense_word1( w, offs, dna, acc | 0, consume_word ) ;
					 acc |= 3 ; break ;
			case 10: make_dense_word1( w, offs, dna, acc | 1, consume_word ) ;
					 acc |= 3 ; break ;
			case 11: make_dense_word1( w, offs, dna, acc | 0, consume_word ) ;
					 make_dense_word1( w, offs, dna, acc | 1, consume_word ) ;
					 acc |= 3 ; break ;
			case 12: make_dense_word1( w, offs, dna, acc | 2, consume_word ) ;
					 acc |= 3 ; break ;
			case 13: make_dense_word1( w, offs, dna, acc | 0, consume_word ) ;
					 make_dense_word1( w, offs, dna, acc | 2, consume_word ) ;
					 acc |= 3 ; break ;
			case 14: make_dense_word1( w, offs, dna, acc | 1, consume_word ) ;
					 make_dense_word1( w, offs, dna, acc | 2, consume_word ) ;
					 acc |= 3 ; break ;
			case 15: make_dense_word1( w, offs, dna, acc | 0, consume_word ) ;
					 make_dense_word1( w, offs, dna, acc | 1, consume_word ) ;
					 make_dense_word1( w, offs, dna, acc | 2, consume_word ) ;
					 acc |= 3 ; break ;
		}
	}

	// No more bases to add, we can finally do something with this word.
	consume_word( offs, acc ) ;
}

template< typename G >
void make_dense_word( int w, unsigned offs, Oligo dna, G consume_word )
{
	make_dense_word1<G>( w, offs, dna, 0, consume_word ) ;
}

class count_word {
	private:
		uint32_t *base ;

	public:
		count_word( uint32_t *b ) : base(b) {}
		void operator()( unsigned off, unsigned ix ) { ++base[ix] ; }
} ;

class store_word {
	private:
		uint32_t *index_1l ;
		uint32_t *index_2l ;

	public:
		store_word( uint32_t *i1, uint32_t *i2 ) : index_1l(i1), index_2l(i2) {}
		void operator()( unsigned off, unsigned ix )
		{ 
			unsigned& p = index_1l[ ix ] ;
			if( p ) index_2l[ --p ] = off ;
		}
} ;

int main_( int argc, const char * const argv[] )
{
	if( argc != 3 ) return 1 ;
	CompactGenome genome( argv[1] ) ;

	unsigned word_size =   10 ; // XXX
	unsigned cutoff    = 1280 ; // XXX

	unsigned first_level_len = 1 << (2 * word_size) + 1 ;
	assert( std::numeric_limits<size_t>::max() / 4 > first_level_len ) ;

	uint32_t *base = (uint32_t*)malloc( 4 * first_level_len ) ;
	throw_errno_if_null( base, "allocating first level index" ) ;
#ifdef _BSD_SOURCE
	madvise( base, 4*first_level_len, MADV_WILLNEED ) ;
#endif

	// we will scan over the dna twice: once to count occurences of
	// words, once to actually store pointers.
	
	// First scan: only count words.  We'll have to go over the whole
	// table again to convert counts into offsets.
	genome.scan_words( word_size, make_dense_word<count_word>, count_word( base ), "Counting" ) ;

#if 0
	// This is temporary... make a histogram.
	std::map< int, int > hist ;
	for( unsigned *p = base+3 ; p != base+3+first_level_len ; ++p ) ++hist[ *p ] ;

	for( std::map< int,int >::const_iterator i = hist.begin() ; i != hist.end() ; ++i )
		std::cout << i->first << '\t' << i->second << std::endl ;
#endif

	// sum up occurences and replace counts by pointers into array.
	// For practical reasons, pointers (actually indices) will be set to
	// the *end* of the allocated range and count down when writing
	// second-level-pointers.  
	//
	// Whenever we decide to ignore some oligo because it occurs too
	// often, we store the index 0.  An additional pass needs to replace
	// the zeros by sensible values.
	unsigned total = 0 ;
	for( uint32_t *p = base ; p != base + first_level_len ; ++p )
	{
		if( *p < cutoff ) {
			unsigned total_ = total ;
			total += *p ;
			assert( total_ <= total ) ;		// this blows up on overflow
			*p = total ;
		}
		else
		{
			*p = 0 ;
		}
	}
	std::clog << "Need to store " << total << " pointers." << std::endl ;
	assert( std::numeric_limits<size_t>::max() / 4 > first_level_len  + total ) ;

	uint32_t *lists = (uint32_t*)malloc( 4 * total ) ;
	throw_errno_if_null( lists, "allocated second level index" ) ;
#ifdef _BSD_SOURCE
	madvise( lists, 4 * total, MADV_WILLNEED ) ;
#endif

	// Second scan: we actually store the offsets now.
	genome.scan_words( word_size, make_dense_word<store_word>, store_word( base, lists ), "Indexing" ) ;

	// need to fix 0-entries in 1L index
	unsigned last = total ;
	for( uint32_t *p = base + first_level_len ; p != base ; )
	{
		--p ;
		if( !*p ) *p = last ;
		else last = *p ;
	}

	std::clog << "Writing " << argv[2] << "..." << std::endl ;
	int fd = open( argv[2], O_RDWR | O_TRUNC | O_CREAT, 0644 ) ;
	throw_errno_if_minus1( fd, "opening", argv[2] ) ;

	unsigned signature = 0x30584449 ; // IDX0 
	write( fd, &signature, 4 ) ;
	write( fd, &first_level_len, 4 ) ;
	write( fd, base, 4 * first_level_len ) ;
	write( fd, &total, 4 ) ;
	write( fd, lists, 4 * total ) ;
	close( fd ) ;

	std::cout << "index \"" << argv[2] << "\" {\n"
		<< "\tmethod compact_word_list ;\n"
		<< "\tgenome \"" << argv[1] << "\" ;\n"
		<< "\twordsize " << word_size << " ;\n"
		<< "\tcutoff " << cutoff << " ;\n"
		<< "\tindexsize " << total << " ;\n"
		<< "} ;\n" ;
	return 0 ;
}

