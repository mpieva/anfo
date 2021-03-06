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

/*! \page DNA indexer
 *
 * Scans a dna file, creating an index of fixed-size words.  We need a
 * two-level index and want to keep it simple.  To this end, we'll store
 * plain arrays in the index files so they can be directly mmaped.  That
 * makes them inherently unportable between different architectures, but
 * we'll have to live with that.  So, this is what an index file looks
 * like:
 *
 * - The signature "IDX1" in ASCII.
 * - The length of the meta information.
 * - The metainformation of type config::CompactIndex, serialized in
 *   binary protocol buffer format.
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
 * word size to no less than 16 bases.  (No problem, more than that is
 * probably impractical anyway.)  
 */

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include "index.h"
#include "util.h"
#include "config.pb.h"

#include <google/protobuf/text_format.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>

#include <popt.h>

#include <iostream>
#include <iomanip>
#include <sstream>

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

//! \brief checks for self-overlapping seeds
//! \internal
//! We define a seed as useless if for any n in the range 1..o dropping
//! the first n bases from the seed results in the same value as
//! dropping the last n bases.
//! \param w seed to test
//! \param l length of the seed 
//! \param o max. overlap length to check
//! \return true iff the seed is good (i.e. doesn't self-overlap)
inline bool check_repeat( uint32_t w, uint32_t l, uint32_t o )
{
	uint32_t mask = (1 << (2*l))-1, v = w ;
	for( ; o ; --o )
	{
		mask >>= 2 ;
		v    >>= 2 ;
		if( v == (w & mask) ) return false ;
	}
	return true ;
}

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
		uint32_t stride_ ;
		uint32_t size_ ;
		uint32_t repeat_ ;
		uint64_t &total_ ;

		count_word( uint32_t *b, uint32_t c, uint32_t s, uint32_t sz, uint32_t r, uint64_t &t )
			: base_(b), cutoff_(c), stride_(s), size_(sz), repeat_(r), total_(t) {}
		void operator()( uint32_t off, uint32_t ix ) const {
			if( off % stride_ == 0 && check_repeat( ix, size_, repeat_ ) ) {
				uint32_t x = ++base_[ix] ; 
				if( x <= cutoff_ ) {
					++total_ ;
					if( x == cutoff_ ) total_ -= cutoff_ ;
				}

				if( std::numeric_limits<size_t>::max() / 4 <= total_ ) 
					throw "Sorry, but a bigger machine is needed for this kind of index." ;
			}
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
		uint32_t stride_ ;
		uint32_t size_ ;
		uint32_t repeat_ ;

	public:
		store_word( uint32_t s, uint32_t *i1, uint32_t *i2, uint32_t sz, uint32_t r )
			: index_1l_(i1), index_2l_(i2), stride_(s), size_(sz), repeat_(r) {}
		void operator()( uint32_t off, uint32_t ix )
		{ 
			if( off % stride_ == 0 && check_repeat( ix, size_, repeat_ ) ) {
				uint32_t& p = index_1l_[ (uint64_t)ix ] ;
				if( p ) index_2l_[ (uint64_t)(--p) ] = off ;
			}
		}
} ;

struct noop { void operator()( uint32_t ) const {} } ;
struct put_gap { 
	google::protobuf::RepeatedField<uint32_t> *gaps_ ;
	put_gap( google::protobuf::RepeatedField<uint32_t> *gaps ) : gaps_(gaps) {}
	void operator()( uint32_t g ) const { gaps_->Add(g) ; } 
} ;

WRAPPED_MAIN
{
	GOOGLE_PROTOBUF_VERIFY_VERSION ;
	enum option_tags { opt_none, opt_version } ;

	const char* output_file = 0 ;
	const char* description = 0 ;
	const char* genome_name = 0 ;

	unsigned wordsize  = 12 ;
	unsigned cutoff    = std::numeric_limits<unsigned>::max() ;
	unsigned repeat    = std::numeric_limits<unsigned>::max() ;
	unsigned stride    = 4 ;
	int      verbose   = 0 ;
	int 	 histogram = 0 ;
	int      slice     = 0 ;
	int      slices    = 1 ;

	config::Config mi ;
	struct poptOption options[] = {
		{ "version",     'V', POPT_ARG_NONE,   0,            opt_version, "Print version number and exit", 0 },
		{ "output",      'o', POPT_ARG_STRING, &output_file, opt_none,    "Output index to FILE", "FILE" },
		{ "genome",      'g', POPT_ARG_STRING, &genome_name, opt_none,    "Read genome GENOME", "GENOME" },
		{ "description", 'd', POPT_ARG_STRING, &description, opt_none,    "Add TEXT as description to index", "TEXT" },
		{ "wordsize",    's', POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT,    &wordsize,    opt_none,    "Index words of length SIZE", "SIZE" },
		{ "stride",      'S', POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT,    &stride,      opt_none,    "Index every Nth word", "N" },
		{ "limit",       'l', POPT_ARG_INT,    &cutoff,      opt_none,    "Do not index words more frequent than LIM", "LIM" },
		{ "repetitive",  'r', POPT_ARG_INT,    &repeat,      opt_none,    "Do not index words that self-overlap N bases", "N" },
		{ "slice",       'x', POPT_ARG_INT,    &slice,       opt_none,    "Process slice N if creating a split index", "N" },
		{ "numslice",       'x', POPT_ARG_INT,    &slice,       opt_none,    "Process slice N if creating a split index", "N" },
		{ "histogram",   'h', POPT_ARG_NONE,   &histogram,   opt_none,    "Produce histogram of word frequencies", 0 },
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
			std::cout << poptGetInvocationName(pc) << ", revision " << PACKAGE_VERSION << std::endl ;
			return 0 ;

		default:
			std::clog << poptGetInvocationName(pc) << ": " << poptStrerror( rc ) 
				<< ' ' << poptBadOption( pc, 0 ) << std::endl ;
			return 1 ; 
	}

	if( !genome_name ) throw "missing --genome option" ;
	if( poptGetArg( pc ) ) throw "unexpected non-option argument" ;

	GenomeHolder genome = Metagenome::find_genome( genome_name ) ;

	if( repeat == std::numeric_limits<uint32_t>::max() ) repeat = wordsize >> 1 ;
	uint64_t first_level_len = (1 << (2 * wordsize)) + 1 ;
	if( std::numeric_limits<uint32_t>::max() <= first_level_len 
			|| std::numeric_limits<size_t>::max() / 4 <= first_level_len ) 
		throw "integer overflow in size calculation; cannot handle this word size" ;

	uint32_t *base = (uint32_t*)malloc( 4 * first_level_len ) ;
	throw_errno_if_null( base, "allocating first level index" ) ;

	// we will scan over the dna twice: once to count occurences of
	// words, once to actually store pointers.
	
	// First scan: only count words.  We'll have to go over the whole
	// table again to convert counts into offsets.
	uint64_t total0 = 0;
	genome->scan_words( wordsize, make_dense_word( count_word( base, cutoff, stride, wordsize, repeat, total0 ) ), noop(), slice, slices, "Counting" ) ;
	std::clog << "Need to store " << total0 << " pointers." << std::endl ;

	// Histogram of word frequencies.
	if( histogram )
	{
		std::map< int, int > hist ;
		for( unsigned *p = base+3 ; p != base+3+first_level_len ; ++p )
			++hist[ *p ] ;

		int max = 0 ;
		for( std::map< int,int >::const_iterator i = hist.begin() ; i != hist.end() ; ++i )
			if( max < i->second ) max = i->second ;

		for( std::map< int,int >::const_iterator i = hist.begin() ; i != hist.end() ; ++i )
			std::cout << i->first << '\t' << i->second << '\t'
				<< std::string( 50*i->second / max, '*' ) << std::endl ;
	}

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

	if( total0 != total ) throw "logic error: totals are not equal" ;

	uint32_t *lists = (uint32_t*)malloc( 4 * total ) ;
	throw_errno_if_null( lists, "allocating second level index" ) ;
	madvise( lists, 4 * total, MADV_WILLNEED ) ;

	config::CompactIndex ci ;

	// Second scan: we actually store the offsets now.
	genome->scan_words(
			wordsize,
			make_dense_word( store_word( stride, base, lists, wordsize, repeat ) ),
			put_gap( ci.mutable_gap() ),
			slice, slices,
			"Indexing" ) ;
	std::clog << "Found " << ci.gap_size() << " gap symbols" << std::endl ;

	// need to fix 0-entries in 1L index
	uint32_t last = total ;
	for( uint32_t *p = base + first_level_len ; p != base ; )
	{
		--p ;
		if( !*p ) *p = last ;
		else last = *p ;
	}

	ci.set_genome_name( genome_name ) ;
	ci.set_wordsize( wordsize ) ;
	ci.set_stride( stride ) ;
	ci.set_indexsize( total ) ;
	if( cutoff != std::numeric_limits<unsigned>::max() ) ci.set_cutoff( cutoff ) ;
	std::string metainfo ;
	if( !ci.SerializeToString( &metainfo ) ) throw "error when serializing metainfo" ;

	// Note to self: creating and mmapping a sparse file, then filling
	// it in random order is an exceptionally bad idea...
 
	std::clog << "Writing " << genome->name() << "..." << std::endl ;
	std::stringstream output_path ;
	if( output_file ) output_path << output_file ;
	else output_path << genome->name() << '_' << wordsize << ".idx" ;
	int fd = throw_errno_if_minus1( creat( output_path.str().c_str(), 0644 ), "opening", output_path.str().c_str() ) ;

	uint32_t sig = FixedIndex::signature, len = metainfo.size() ; 
	metainfo.append( "\0\0\0" ) ;
	mywrite( fd, &sig, 4 ) ;
	mywrite( fd, &len, 4 ) ;
	mywrite( fd, metainfo.data(), (len+3)&~3 ) ;
	mywrite( fd, base, 4 * first_level_len ) ;
	mywrite( fd, lists, 4 * total ) ;

	close( fd ) ;
	return 0 ;
}


