#include "config.h"

//! \page dump_params Parameter Extraction
//! dump-params expects exactly one file argument, which must be either
//! an ANFO output file or a textual configuration file.  The file is
//! opened, a configuration extracted, and if it contains the parameters
//! for an ancient DNA aligner, the aligner is configured and its
//! internal parameters (substitution matrices and associated penalties)
//! are dumped to stdout.

#include "align.h"
#include "conffile.h"
#include "stream.h"
#include "util.h"

#include "output.pb.h"
#include <iostream>
#include <memory>

int main_( int argc, const char** argv )
{
	if( argc != 2 ) throw "expected exactly one argument" ;

	config::Config conf ;
	try {
		std::auto_ptr< streams::Stream > af( streams::make_input_stream( argv[1] ) ) ;
		conf = af->fetch_header().config() ;
	} 
	catch( const streams::AnfoReader::ParseError& )
	{
		conf = parse_text_config( argv[1] ) ;
	}
	
	if( !conf.has_aligner() ) throw "couldn't find aligner configuration" ;
	simple_adna::configure( conf.aligner(), &std::cout ) ;
	return 0 ;
}

