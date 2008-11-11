/*
 * FASTA to DNA file converter
 *
 * Reads a FASTA file with one ore more DNA or RNA sequences possibly
 * containing ambiguity codes from stdin, writes it in a compact
 * representation under the filename given as first argument, and writes
 * a human readable index to stdout, which should be stored.
 *
 * We chop the genome into contigs.  Every sequence start also starts a
 * new contig, stretches of more than 2 Ns separate contigs and are not
 * stored.  Ambiguity codes are recognized.
 *
 * A "compact genome" is a file encoding DNA as 4 bits per nucleotide.
 * The bits encode for A,C,T,G from LSB to MSB.  Combinations of bits
 * form ambiguity codes (0 is a gap, used as terminator, 15 is an N),
 * two left rotations reverse-complement a single nucleotide.  In a
 * byte, the first nucleotide occupies the 4 LSBs.
 *
 * The "compact genome" starts with the ACSII-signature "DNA0", followed
 * by the raw DNA.  Each contig is prepended with a single gap, then the
 * contigs are concatenated, then a final single gap is appended.  The
 * DNA is written out with two nucleotides per byte.  The final byte is
 * padded with a gap if necessary.  Note that every contig is always
 * preceded and succeeded by a gap; this is done so we can mmap(2) a
 * genome into memory and operate on pointers into it in the knowledge
 * that every contig is terminated by a gap in either direction.
 *
 * An index of the contigs and a description of the genome is packed
 * into a message of type "Genome" (see metaindex.proto) and written to
 * stdout in protobuf text format.
 */

#include "metaindex.pb.h"
#include "util.h"

#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/text_format.h>

#include <popt.h>

#include <cctype>
#include <fstream>
#include <iostream>

#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>

uint8_t decode_fna( char c ) {
	switch( c ) {
		case 'a': case 'A': return 1 ;
		case 'b': case 'B': return 1 ^ 0xf ;
		case 'c': case 'C': return 2 ;
		case 'd': case 'D': return 2 ^ 0xf ;
		case 'g': case 'G': return 8 ;
		case 'h': case 'H': return 8 ^ 0xf ;
		case 't': case 'T': return 4 ;
		case 'u': case 'U': return 4 ;
		case 'v': case 'V': return 4 ^ 0xf ;

		case 'n': case 'N': return 15 ;
		case 'm': case 'M': return  3 ;
		case 'k': case 'K': return 12 ;
		case 'w': case 'W': return  5 ;
		case 's': case 'S': return 10 ;
		case 'r': case 'R': return  9 ;
		case 'y': case 'Y': return  6 ;
		default:            return  0 ;
	}
}

inline bool encodes_nuc( char c ) { return decode_fna(c) != 0 ; }
inline uint8_t to_mask( char c ) { return decode_fna(c) & 0xf ; }

class FastaDecoder
{
	private:
		void (FastaDecoder::*state)( char ) ;

		unsigned position ;
		unsigned source_position ;
		unsigned num_n  ;
		uint8_t one_nt ;
		const unsigned max_num_n ;
		bool verbose ;

		std::ofstream& dna ;
		metaindex::Genome *cur_genome ;
		metaindex::Sequence *cur_sequence ;
		metaindex::Contig *cur_contig ;
		
		
	public:
		FastaDecoder( unsigned p, std::ofstream& dna_, metaindex::Genome *g_, unsigned maxn, bool v_ = false )
			: state( &FastaDecoder::s_idle ), position( p ), num_n( 0 ), one_nt( 0 )
			, max_num_n( maxn ), verbose( v_ ), dna( dna_ )
			, cur_genome( g_ ), cur_sequence( 0 ), cur_contig( 0 )
		{
			dna.write( "DNA0", 4 ) ;
		}

		void step( char c ) { (this->*state)( c ) ; }

		void begin_new_sequence()
		{
			source_position = 0 ;
			cur_sequence = cur_genome->add_sequence() ;
		}
		void finish_sequence() { cur_sequence = 0 ; }
		void begin_new_contig() 
		{ 
			put_nt( '-' ) ;
			num_n = 0 ;
			cur_contig = cur_sequence->add_contig() ;
			cur_contig->set_offset( position ) ;
			cur_contig->set_range_start( source_position ) ;
			if( verbose )
				std::clog << cur_sequence->name() << " @ " << source_position << std::endl ;
		}
		void finish_contig()
		{ 
			cur_contig->set_range_end( source_position ) ;
		}
		void finish_genome()
		{
			put_nt( '-' ) ;
			cur_genome->set_total_size( position ) ;
			if( position & 1 ) put_nt( '-' ) ;
		}

		void put_nt( char c )
		{
			if( c == 'n' || c == 'N' ) 
			{
				++num_n ;
			}
			else
			{
				for( ; num_n ; --num_n ) put_nt_( '-' ) ;
				put_nt_( c ) ;
			}
		}

		void put_nt_( char c )
		{
			if( position & 1 )
			{
				dna.put( (to_mask(c) << 4) | one_nt ) ;
			}
			else
			{
				one_nt = to_mask(c) ;
			}
			++position ;
		}

		void s_idle( char c )
		{
			if( c == '>' ) {
				begin_new_sequence() ;
				state = &FastaDecoder::s_name ;
			}
		}

		void s_name( char c )
		{
			if( !c ) finish_genome() ;
			else if( std::isspace(c) || c == '\n' )
			{
				if( c == '\n' ) state = &FastaDecoder::s_seq_begin ;
				else state = &FastaDecoder::s_info ;
			}
			else cur_sequence->mutable_name()->push_back( c ) ;
		}

		void s_info( char c )
		{
			if( !c ) finish_genome() ; 
			else if( c == '\n' ) state = &FastaDecoder::s_seq_begin ;
		}

		void s_seq_begin( char c )
		{
			if( !c ) finish_genome() ;
			else if( c == ';' ) state = &FastaDecoder::s_info ;
			else if( c == '>' ) {
				begin_new_sequence() ;
				state = &FastaDecoder::s_name ;
			}
			else if( c == 'n' || c == 'N' ) {
				if( num_n == max_num_n ) {
					source_position += num_n+1 ;
					num_n = 0 ;
					state = &FastaDecoder::s_contig_gap ;
				}
				else
				{
					++source_position ;
					put_nt( c ) ;
				}
			}
			else if( encodes_nuc( c ) ) {
				begin_new_contig() ;
				++source_position ;
				put_nt( c ) ;
				state = &FastaDecoder::s_sequence ;
			}
		}

		void s_sequence( char c )
		{
			if( !c ) {
				finish_contig() ;
				finish_sequence() ;
				finish_genome() ;
			}
			else if( c == 'N' || c == 'n' && num_n == max_num_n ) {
				finish_contig() ;
				source_position += num_n+1 ;
				num_n = 0 ;
				state = &FastaDecoder::s_contig_gap ;
			}
			else if( encodes_nuc( c ) ) {
				++source_position ;
				put_nt( c ) ;
			}
			else if( c == '>' ) {
				finish_contig() ;
				finish_sequence() ;
				begin_new_sequence() ;
				state = &FastaDecoder::s_name ;
			}
		}

		void s_contig_gap( char c ) 
		{
			if( !c ) {
				finish_sequence() ;
				finish_genome() ;
			}
			else if( c == '>' ) {
				finish_sequence() ;
				begin_new_sequence() ;
				state = &FastaDecoder::s_name ;
			}
			else if( encodes_nuc( c ) && c != 'n' && c != 'N' ) {
				begin_new_contig() ;
				++source_position ;
				put_nt( c ) ;
				state = &FastaDecoder::s_sequence ;
			}
			else if( c == 'n' || c == 'N' ) {
				++source_position ;
			}
		}
} ;




int main_( int argc, const char * argv[] )
{
	enum option_tags { opt_none, opt_version } ;

	const char* output_file = 0 ;
	const char* config_file = 0 ;
	const char* description = 0 ;
	const char* genome_name = 0 ;
	int max_num_n = 2 ;
	int verbose = 0 ;

	struct poptOption options[] = {
		{ "version", 'V', POPT_ARG_NONE,   0, opt_version, "print version number and exit", 0 },
		{ "output",  'o', POPT_ARG_STRING, &output_file, opt_none, "write output to FILE", "FILE" },
		{ "maxn",    'm', POPT_ARG_INT,    &max_num_n, opt_none, "treat N consecutive Ns as separator", "N" },
		{ "config",  'c', POPT_ARG_STRING, &config_file, opt_none, "use FILE as configuration", "FILE" },
		{ "genome",  'g', POPT_ARG_STRING, &genome_name, opt_none, "set genome name to NAME", "NAME" },
		{ "description", 'd', POPT_ARG_STRING, &description, opt_none, "add TEXT as description to genome", "TEXT" },
		{ "verbose", 'v', POPT_ARG_NONE,   &verbose, opt_none, "make more noise while working", 0 },
		POPT_AUTOHELP POPT_TABLEEND
	} ;

	poptContext pc = poptGetContext( "fa2dna", argc, argv, options, 0 ) ;
	poptSetOtherOptionHelp( pc, "[OPTION...] [FASTA-file...]" ) ;
	poptReadDefaultConfig( pc, 0 ) ;

	if( argc <= 1 ) { poptPrintHelp( pc, stderr, 0 ) ; return 1 ; }
	for( int rc = poptGetNextOpt( pc ) ; rc > 0 ; rc = poptGetNextOpt(pc) ) switch( rc )
	{
		case opt_version:
			std::cout << poptGetInvocationName(pc) << ", revision " << VERSION << std::endl ;
			return 0 ;

		default:
			std::cerr << argv[0] << ": " << poptStrerror( rc ) 
				<< ' ' << poptBadOption( pc, 0 ) << std::endl ;
			return 1 ; 
	}

	if( !config_file ) throw "missing --config option" ;
	if( !output_file ) throw "missing --output option" ;
	if( !genome_name ) throw "missing --genome option" ;

	metaindex::MetaIndex mi ;
	int config_in = open( config_file, O_RDONLY ) ;
	if( config_in != -1 || errno != ENOENT )
	{
		throw_errno_if_minus1( config_in, "reading config" ) ;
		google::protobuf::io::FileInputStream fis( config_in ) ;
		google::protobuf::TextFormat::Parse( &fis, &mi ) ;
		close( config_in ) ;
	}

	metaindex::Genome *g = 0 ;
	for( int i = 0 ; i != mi.genome_size() && !g ; ++i )
		if( mi.genome(i).name() == genome_name )
			g = mi.mutable_genome(i) ;
	if( g ) g->clear_sequence() ; else g = mi.add_genome() ;

	g->set_name( genome_name ) ;
	g->set_filename( output_file ) ;
	if( description ) g->set_description( description ) ;

	std::ofstream output_stream( output_file ) ;
	FastaDecoder fd( 8, output_stream, g, max_num_n, verbose ) ;
	if( !poptPeekArg( pc ) ) 
		for( int c ; (c = std::cin.get()) != std::istream::traits_type::eof() ; )
			fd.step( c ) ;
	
	while( const char* arg = poptGetArg( pc ) )
	{
		if( strcmp(arg,"-") ) 
		{
			std::ifstream str( arg ) ;
			for( int c ; (c = str.get()) != std::istream::traits_type::eof() ; )
				fd.step( c ) ;
		}
		else
			for( int c ; (c = std::cin.get()) != std::istream::traits_type::eof() ; )
				fd.step( c ) ;
	}

	fd.step( 0 ) ;

	std::string new_config_file = config_file ;
	new_config_file.push_back('#') ;

	int config_out = throw_errno_if_minus1(
			creat( new_config_file.c_str(), 0644 ), "writing config" ) ;
	google::protobuf::io::FileOutputStream fos( config_out ) ;
	fos.SetCloseOnDelete( true ) ;
	google::protobuf::TextFormat::Print( mi, &fos ) ;

	throw_errno_if_minus1( 
			rename( new_config_file.c_str(), config_file ),
			"renaming config" ) ;
	return 0 ;
}


