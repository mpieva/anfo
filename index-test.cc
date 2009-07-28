//! \page index_test Testing the indexing phase
//! This binary just runs the indexer on whatever input it is given,
//! printing the resulting seeds in various detail.  If the input looks
//! as if it contains a coordinate (simulated data does), this program
//! will check whether a seed in the correct region is found at all.

#include "align.h"
#include "anfo_common.h"
#include "compress_stream.h"
#include "conffile.h"
#include "index.h"
#include "util.h"

#include <popt.h>

#include <limits>
#include <string>
#include <iostream>
#include <iomanip>

using namespace config ;
using namespace std ;
using namespace google::protobuf::io ;


int main_( int argc, const char * argv[] )
{
	GOOGLE_PROTOBUF_VERIFY_VERSION ;
	enum option_tags { opt_none, opt_version, opt_quiet } ;

	const char* config_file = 0 ;
	int outputlevel = 0 ;
	int simulation = 0 ;

	struct poptOption options[] = {
		{ "version",     'V', POPT_ARG_NONE,   0,            opt_version, "Print version number and exit", 0 },
		{ "config",      'c', POPT_ARG_STRING, &config_file, opt_none,    "Read config from FILE", "FILE" },
		{ "simulation",  's', POPT_ARG_NONE,   &simulation,  opt_none,    "Verify simulated data", 0 },
		{ "debug",       'd', POPT_ARG_INT,    &outputlevel, opt_none,    "Set debug level to L", "L" },
		POPT_AUTOHELP POPT_TABLEEND
	} ;

	poptContext pc = poptGetContext( "anfo", argc, argv, options, 0 ) ;
	poptSetOtherOptionHelp( pc, "[OPTION...] [sequence-file...]" ) ;
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

	Mapper mapper( get_default_config( config_file ) ) ;

	unsigned total_seeded = 0, total_failures = 0 ;
	while( const char* arg = poptGetArg( pc ) ) 
	{
		int inp_fd = !strcmp( arg, "-" ) ? 0 :
			throw_errno_if_minus1( open( arg, O_RDONLY ), "opening ", arg ) ;

		FileInputStream raw_inp( inp_fd ) ;
		std::auto_ptr<ZeroCopyInputStream> inp( decompress( &raw_inp ) ) ;

		QSequence ps ;
		while( read_fastq( inp.get(), ps ) ) 
		{
			output::Result r ;
			std::deque< alignment_type > ol ;
			mapper.index_sequence( ps, r, ol ) ;

			uint32_t closest = UINT_MAX ;
			uint32_t seq_len = ps.length() ;
			if( simulation )
			{
				size_t p0 = ps.get_name().find( '-', 4 ) ;
				size_t p1 = ps.get_name().find( '_', p0 ) ;

				std::string seq_name = ps.get_name().substr( 4, p0-4 ) ;
				uint32_t orig_pos = atoi( ps.get_name().substr( p0+1, p1-p0-1 ).c_str() ) ;
				if( ps.get_name()[p1+1] == '+' ) orig_pos += seq_len / 2 ;
				else orig_pos += seq_len / 2 ;

				for( size_t i = 0 ; i != ol.size() ; ++i )
				{
					uint32_t offset ;
					const Sequence *sequ ;
					mapper.translate_to_genome_coords( ol[i].reference, offset, &sequ ) ;
					if( sequ->name() == seq_name && orig_pos - offset < closest )
						closest = orig_pos - offset ;
					if( sequ->name() == seq_name && offset - orig_pos < closest )
						closest = offset - orig_pos ;
				}
			}

			if( r.has_aln_stats() )
			{
				++total_seeded ;
				cout << setw(27) << ps.get_name()
					<< setw(5) << seq_len
					<< setw(10) << r.aln_stats().num_raw_seeds()
					<< setw(10) << r.aln_stats().num_useless()
					<< setw(10) << r.aln_stats().num_grown_seeds()
					<< setw(10) << r.aln_stats().num_clumps() ;

				if( simulation )
				{
					if( closest <= seq_len / 2 ) std::cout << "  got correct seed" ;
					else {
						++total_failures ;
						if( closest == UINT_MAX ) std::cout << "  no seed" ;
						else std::cout << "  closest seed: " << closest ;
					}
				}
				std::cout << std::endl ;

				for( size_t i = 0 ; i != ol.size() ; ++i )
				{
					uint32_t offset ;
					const Sequence *sequ ;
					mapper.translate_to_genome_coords( ol[i].reference, offset, &sequ ) ;

					if( outputlevel >= 1 )
						std::cout << sequ->name() << "+" << offset << ", " ;
				}
				if( outputlevel >= 1 )
					std::cout << std::dec << std::endl ;
			}
		}
		if( inp_fd ) close( inp_fd ) ;
	}
	if( simulation )
		cout << total_failures << " failures in " << total_seeded << " attempts." << endl ;

	poptFreeContext( pc ) ;
	return 0 ;
}
