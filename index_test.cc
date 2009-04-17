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

using namespace config ;
using namespace std ;
using namespace google::protobuf::io ;


int main_( int argc, const char * argv[] )
{
	GOOGLE_PROTOBUF_VERIFY_VERSION ;
	enum option_tags { opt_none, opt_version, opt_quiet } ;

	const char* config_file = 0 ;
	int outputlevel = 0 ;

	struct poptOption options[] = {
		{ "version",     'V', POPT_ARG_NONE,   0,            opt_version, "Print version number and exit", 0 },
		{ "config",      'c', POPT_ARG_STRING, &config_file, opt_none,    "Read config from FILE", "FILE" },
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
			// XXX

			cout << ps.get_name() << ": " << r.num_raw_seeds() << " seeds, "
				<< r.num_grown_seeds() << " superseeds, "
				<< r.num_clumps() << " aggregates." << endl ;
			if( outputlevel >= 1 )
			{
				for( size_t i = 0 ; i != ol.size() ; ++i )
					std::cout << ol[i].reference << ", " ;
				std::cout << std::dec << std::endl ;
			}
		}
		if( inp_fd ) close( inp_fd ) ;
	}

	poptFreeContext( pc ) ;
	return 0 ;
}
