/*! \page DNA indexer
 *
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
 *
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
 *
 * \todo This is crap... splitting everything over two files.  DNA and
 *       index files should include their metainformation.  Maybe put an
 *       appropriate message at the end of the file, store an offset at
 *       the start?
 *
 * \note to self: creating and mmapping a sparse file, then filling it in
 *       random order is an exceptionally bad idea...
 */

#include "index.h"
#include "util.h"
#include "metaindex.pb.h"
#include "conffile.h"

#include <google/protobuf/text_format.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>

#include <popt.h>

#include <cassert>
#include <iostream>
#include <iomanip>

#include <fcntl.h>
#include <unistd.h>

template< typename G > struct mk_dense_word {
	public:
		mk_dense_word( G consume_word ) : consume_( consume_word ) {}

		void operator()( int w, uint32_t offs, Oligo dna ) { go( w, offs, dna, 0 ) ; }

	private:
		G consume_ ;

		void go( int w, uint32_t offs, Oligo dna, uint32_t acc ) 
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
					case  0: break ; // not supposed actually occur
					case  1: acc |= 0 ; break ;
					case  2: acc |= 1 ; break ;
					case  3: go( w, offs, dna, acc | 0 ) ;
							 acc |= 1 ; break ;
					case  4: acc |= 2 ; break ;
					case  5: go( w, offs, dna, acc | 0 ) ;
							 acc |= 2 ; break ;
					case  6: go( w, offs, dna, acc | 1 ) ;
							 acc |= 2 ; break ;
					case  7: go( w, offs, dna, acc | 0 ) ;
							 go( w, offs, dna, acc | 1 ) ;
							 acc |= 2 ; break ;
					case  8: acc |= 3 ; break ;
					case  9: go( w, offs, dna, acc | 0 ) ;
							 acc |= 3 ; break ;
					case 10: go( w, offs, dna, acc | 1 ) ;
							 acc |= 3 ; break ;
					case 11: go( w, offs, dna, acc | 0 ) ;
							 go( w, offs, dna, acc | 1 ) ;
							 acc |= 3 ; break ;
					case 12: go( w, offs, dna, acc | 2 ) ;
							 acc |= 3 ; break ;
					case 13: go( w, offs, dna, acc | 0 ) ;
							 go( w, offs, dna, acc | 2 ) ;
							 acc |= 3 ; break ;
					case 14: go( w, offs, dna, acc | 1 ) ;
							 go( w, offs, dna, acc | 2 ) ;
							 acc |= 3 ; break ;
					case 15: go( w, offs, dna, acc | 0 ) ;
							 go( w, offs, dna, acc | 1 ) ;
							 go( w, offs, dna, acc | 2 ) ;
							 acc |= 3 ; break ;
				}
			}

			// No more bases to add, we can finally do something with this word.
			consume_( offs, acc ) ;
		}
} ;

//! \brief makes C++ happy
//! Type inference only works for functions, not objects.  Don't you
//! just love C++?
//! \internal
template< typename G > mk_dense_word<G> make_dense_word( G g ) { return mk_dense_word<G>( g ) ; }

//! \brief counts words in a genome
//! This is passed as a worker to CompactGenome::scan_words, it only
//! counts the words and stores the counts in what later becomes the
//! first level index.
//! \internal
class count_word {
	private:
	public:
		uint32_t *base_ ;
		uint32_t cutoff_ ;
		uint64_t &total_ ;

		count_word( uint32_t *b, uint32_t c, uint64_t &t ) : base_(b), cutoff_(c), total_(t) {}
		void operator()( uint32_t off, uint32_t ix ) const {
			uint32_t x = ++base_[ix] ; 
			if( x <= cutoff_ ) {
				++total_ ;
				if( x == cutoff_ ) total_ -= cutoff_ ;
			}

			if( std::numeric_limits<size_t>::max() / 4 <= total_ ) 
				throw "Sorry, but a bigger machine is needed for this kind of index." ;
		}
} ;

//! \brief stores words in an index
//! This is passed as a worker to CompactGenome::scan_words, it stores
//! the words in the second level index, which must have been correctly
//! allocated beforehand.
//! \internal
class store_word {
	private:
		uint32_t *index_1l_ ;
		uint32_t *index_2l_ ;

	public:
		store_word( uint32_t *i1, uint32_t *i2 ) : index_1l_(i1), index_2l_(i2) {}
		void operator()( uint32_t off, uint32_t ix )
		{ 
			uint32_t& p = index_1l_[ (uint64_t)ix ] ;
			if( p ) index_2l_[ (uint64_t)(--p) ] = off ;
		}
} ;

int main_( int argc, const char * argv[] )
{
	enum option_tags { opt_none, opt_version } ;

	const char* output_file = 0 ;
	const char* config_file = 0 ;
	const char* description = 0 ;
	const char* genome_conf = 0 ;

	unsigned wordsize = 10 ;
	unsigned cutoff   = std::numeric_limits<unsigned>::max() ;
	int      verbose  = 0 ;

	struct poptOption options[] = {
		{ "version",     'V', POPT_ARG_NONE,   0,            opt_version, "Print version number and exit", 0 },
		{ "output",      'o', POPT_ARG_STRING, &output_file, opt_none,    "Output index to FILE", "FILE" },
		{ "config",      'c', POPT_ARG_STRING, &config_file, opt_none,    "Write configuration to FILE", "FILE" },
		{ "genome",      'g', POPT_ARG_STRING, &genome_conf, opt_none,    "Read genome config from FILE", "FILE" },
		{ "description", 'd', POPT_ARG_STRING, &description, opt_none,    "Add TEXT as description to index", "TEXT" },
		{ "wordsize",    's', POPT_ARG_INT,    &wordsize,    opt_none,    "Index words of length SIZE", "SIZE" },
		{ "limit",       'l', POPT_ARG_INT,    &cutoff,      opt_none,    "Do not index words more frequent than LIM", "LIM" },
		{ "verbose",     'v', POPT_ARG_NONE,   &verbose,     opt_none,    "Make more noise while working", 0 },
		POPT_AUTOHELP POPT_TABLEEND
	} ;

	poptContext pc = poptGetContext( "dnaindex", argc, argv, options, 0 ) ;
	poptSetOtherOptionHelp( pc, "[OPTION...]" ) ;
	poptReadDefaultConfig( pc, 0 ) ;

	if( argc <= 1 ) { poptPrintHelp( pc, stderr, 0 ) ; return 1 ; }
	for( int rc = poptGetNextOpt( pc ) ; rc > 0 ; rc = poptGetNextOpt(pc) ) switch( rc )
	{
		case opt_version:
			std::cout << poptGetInvocationName(pc) << ", revision " << VERSION << std::endl ;
			return 0 ;

		default:
			std::cerr << poptGetInvocationName(pc) << ": " << poptStrerror( rc ) 
				<< ' ' << poptBadOption( pc, 0 ) << std::endl ;
			return 1 ; 
	}

	if( !config_file ) throw "missing --config option" ;
	if( !output_file ) throw "missing --output option" ;
	if( !genome_conf ) throw "missing --genome option" ;

	metaindex::Config mi ;
	merge_binary_config( genome_conf, mi ) ;

	if( !mi.genome_size() ) throw "no genome found" ;
	CompactGenome genome( mi.genome(0), MADV_SEQUENTIAL ) ;

	uint64_t first_level_len = 1 << (2 * wordsize) + 1 ;
	assert( std::numeric_limits<size_t>::max() / 4 > first_level_len ) ;

	uint32_t *base = (uint32_t*)malloc( 4 * first_level_len ) ;
	throw_errno_if_null( base, "allocating first level index" ) ;
	madvise( base, 4*first_level_len, MADV_WILLNEED ) ;

	// we will scan over the dna twice: once to count occurences of
	// words, once to actually store pointers.
	
	// First scan: only count words.  We'll have to go over the whole
	// table again to convert counts into offsets.
	uint64_t total0 = 0;
	genome.scan_words( wordsize, make_dense_word( count_word( base, cutoff, total0 ) ), "Counting" ) ;
	std::clog << "Need to store " << total0 << " pointers." << std::endl ;

#if 0
	// XXX This is temporary... make a histogram.
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

	assert( total0 == total ) ;
	uint32_t *lists = (uint32_t*)malloc( 4 * total ) ;
	throw_errno_if_null( lists, "allocating second level index" ) ;
	madvise( lists, 4 * total, MADV_WILLNEED ) ;

	// Second scan: we actually store the offsets now.
	genome.scan_words( wordsize, make_dense_word( store_word( base, lists ) ), "Indexing" ) ;

	// need to fix 0-entries in 1L index
	uint32_t last = total ;
	for( uint32_t *p = base + first_level_len ; p != base ; )
	{
		--p ;
		if( !*p ) *p = last ;
		else last = *p ;
	}

	std::clog << "Writing " << output_file << "..." << std::endl ;
	int fd = throw_errno_if_minus1( creat( output_file, 0644 ), "opening", output_file ) ;

	uint32_t sig = FixedIndex::signature ;
	mywrite( fd, &sig, 4 ) ;
	mywrite( fd, &first_level_len, 4 ) ;
	mywrite( fd, base, 4 * first_level_len ) ;
	mywrite( fd, &total, 4 ) ;
	mywrite( fd, lists, 4 * total ) ;
	close( fd ) ;

	metaindex::Config c ;
	metaindex::CompactIndex *ci = c.add_compact_index() ;
	ci->set_filename( output_file ) ;
	ci->set_genome_name( mi.genome(0).name() ) ;
	ci->set_wordsize( wordsize ) ;
	if( cutoff == std::numeric_limits<unsigned>::max() ) ci->clear_cutoff() ;
	else ci->set_cutoff( cutoff ) ;
	ci->set_indexsize( total ) ;

	write_binary_config( config_file, c ) ;
	return 0 ;
}


