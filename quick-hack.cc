//! \page anfo_for_grid_engine ANFO executable for the SGE
//! This is the version of ANFO hacked to run smoothly on hostile
//! environments, e.g. our broken setup of the Sun Grid Engine.  These
//! are the concessions to the environment:
//!
//! * This version is single threaded (doesn't even link pthreads).
//! * It can slice the input according to the value of SGE_TASK_ID.
//! * It installs a signal handler to report punishment by the SGE.
//! * It writes a new file and renames it at the end so that improper
//!   shutdown is detectable.
//! * It never writes anything to stderr or stdout, unless requested.

#include "config.h"

#include "anfo_common.h"
#include "compress_stream.h"
#include "conffile.h"
#include "stream.h"
#include "util.h"

#include "config.pb.h"
#include "output.pb.h"

#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/text_format.h>

#include <cstring>
#include <csignal>
#include <limits>
#include <memory>
#include <sstream>
#include <string>

using namespace config ;
using namespace output ;
using namespace std ;
using namespace google::protobuf::io ;
using namespace streams ;

int main_( int argc, const char * argv[] )
{
	GOOGLE_PROTOBUF_VERIFY_VERSION ;
	console.loglevel = Console::info ;

	if( argc != 4 ) return 1 ;
	const char* config_file = argv[1] ;
	const char* input_file  = argv[2] ;
	const char* output_file = argv[3] ;

	Config conf = get_default_config( config_file ) ;
	Mapper mapper( conf ) ;

	std::string output_file_ = output_file ;
	output_file_ += ".#new#" ;

	google::protobuf::io::FileOutputStream fos( strcmp( output_file, "-" ) ?
			throw_errno_if_minus1( open( output_file_.c_str(), O_WRONLY | O_CREAT, 0666 ),
				                   "opening", output_file_.c_str() ) : 1 ) ;
	if( strcmp( output_file, "-" ) ) fos.SetCloseOnDelete( true ) ;
	std::auto_ptr< google::protobuf::io::ZeroCopyOutputStream > zos( compress_fast( &fos ) ) ;

	Header ohdr ;
	CodedOutputStream cos( zos.get() ) ;
	cos.WriteRaw( "ANFO", 4 ) ; // signature
	*ohdr.mutable_config() = conf ;
	ohdr.set_version( PACKAGE_VERSION ) ;

	for( const char **arg = argv ; arg != argv+argc ; ++arg ) *ohdr.add_command_line() = *arg ;
	streams::write_delimited_message( cos, 1, ohdr ) ;

	for( auto_ptr<Stream> inp( make_input_stream( input_file ) ) ; inp->get_state() == Stream::have_output ; )
	{
		Result r = inp->fetch_result() ;
		std::deque< alignment_type > ol ;
		QSequence ps ;
		int pmax = mapper.index_sequence( r, ps, ol ) ;
		if( pmax != INT_MAX ) mapper.process_sequence( ps, pmax, ol, r ) ;
		streams::write_delimited_message( cos, 4, r ) ;
	}

	Footer ofoot ;
	ofoot.set_exit_code( exit_with ) ;
	streams::write_delimited_message( cos, 3, ofoot ) ;

	if( !exit_with && strcmp( output_file, "-" ) )
		throw_errno_if_minus1(
				rename( output_file_.c_str(), output_file ),
				"renaming output" ) ;
	return 0 ;
}

