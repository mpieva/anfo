#include "Index.h"
#include "util.h"
#include "metaindex.pb.h"

#include <google/protobuf/text_format.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>

#include <cassert>
#include <iostream>
#include <iomanip>

#include <fcntl.h>
#include <unistd.h>

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

template< typename G >
void make_dense_word1( int w, uint32_t offs, Oligo dna, uint32_t acc, G consume_word )
{
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
void make_dense_word( int w, uint32_t offs, Oligo dna, G consume_word )
{
	make_dense_word1<G>( w, offs, dna, 0, consume_word ) ;
}

class count_word {
	private:
		uint32_t *base ;

	public:
		count_word( uint32_t *b ) : base(b) {}
		void operator()( uint32_t off, uint32_t ix ) { ++base[ix] ; }
} ;

class store_word {
	private:
		uint32_t *index_1l ;
		uint32_t *index_2l ;

	public:
		store_word( uint32_t *i1, uint32_t *i2 ) : index_1l(i1), index_2l(i2) {}
		void operator()( uint32_t off, uint32_t ix )
		{ 
			uint32_t& p = index_1l[ (uint64_t)ix ] ;
			if( p ) index_2l[ (uint64_t)(--p) ] = off ;
		}
} ;

int main_( int argc, const char * const argv[] )
{
	if( argc != 3 ) return 1 ;
	CompactGenome genome( argv[1] ) ;

	unsigned word_size =     8 ; // XXX
	unsigned cutoff    =     std::numeric_limits<unsigned>::max() ; // XXX

	uint64_t first_level_len = 1 << (2 * word_size) + 1 ;
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
	uint64_t total = 0 ;
	for( uint32_t *p = base ; p != base + first_level_len ; ++p )
	{
		if( *p < cutoff ) {
			total += *p ;
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
	uint32_t last = total ;
	for( uint32_t *p = base + first_level_len ; p != base ; )
	{
		--p ;
		if( !*p ) *p = last ;
		else last = *p ;
	}

	std::clog << "Writing " << argv[2] << "..." << std::endl ;
	int fd = open( argv[2], O_RDWR | O_TRUNC | O_CREAT, 0644 ) ;
	throw_errno_if_minus1( fd, "opening", argv[2] ) ;

	uint32_t sig = FixedIndex::signature ;
	mywrite( fd, &sig, 4 ) ;
	mywrite( fd, &first_level_len, 4 ) ;
	mywrite( fd, base, 4 * first_level_len ) ;
	mywrite( fd, &total, 4 ) ;
	mywrite( fd, lists, 4 * total ) ;
	close( fd ) ;

	metaindex::MetaIndex mi ;
	metaindex::CompactIndex *ci = mi.add_compact_index() ;
	ci->set_filename( argv[2] ) ;
	ci->set_genome_name( argv[1] ) ;
	ci->set_wordsize( word_size ) ;
	if( cutoff == std::numeric_limits<unsigned>::max() ) ci->clear_cutoff() ;
	else ci->set_cutoff( cutoff ) ;
	ci->set_indexsize( total ) ;

	// std::cout << "index \"" << argv[2] << "\" {\n"
		// << "\tmethod compact_word_list ;\n"
		// << "\tgenome \"" << argv[1] << "\" ;\n"
		// << "\twordsize " << word_size << " ;\n"
		// << "\tcutoff " << cutoff << " ;\n"
		// << "\tindexsize " << total << " ;\n"
		// << "} ;\n" ;

	google::protobuf::io::OstreamOutputStream oos( &std::cout ) ;
	google::protobuf::TextFormat::Print( mi, &oos ) ;
	return 0 ;
}


