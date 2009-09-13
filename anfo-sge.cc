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

#include <popt.h>

#include <cstring>
#include <csignal>
#include <limits>
#include <sstream>
#include <string>

#if HAVE_UNISTD_H
#include <unistd.h>
#endif

#if HAVE_SYS_STAT_H
#include <sys/stat.h>
#endif

#if HAVE_FCNTL_H
#include <fcntl.h>
#endif

using namespace config ;
using namespace std ;
using namespace google::protobuf::io ;

extern "C" RETSIGTYPE sig_handler( int sig ) { exit_with = sig + 128 ; }
	
void expand_placeholcer( string &s, int x )
{
	stringstream ss ; ss << x ;
	for( size_t p ; string::npos != (p = s.find( "$$" )) ; )
		s.replace( p, 2, ss.str() ) ;
}

int main_( int argc, const char * argv[] )
{
	GOOGLE_PROTOBUF_VERIFY_VERSION ;
	enum option_tags { opt_none, opt_version, opt_quiet } ;

	const char* config_file = 0 ;
	const char* output_file = 0 ; 
	int solexa_scale = 0 ;
	int stride = 1 ;
	int fastq_origin = 33 ;

	Metagenome::nommap = true ;

	struct poptOption options[] = {
		{ "version",     'V', POPT_ARG_NONE,   0,            opt_version, "Print version number and exit", 0 },
		{ "config",      'c', POPT_ARG_STRING, &config_file, opt_none,    "Read config from FILE", "FILE" },
		{ "output",      'o', POPT_ARG_STRING, &output_file, opt_none,    "Write output to FILE", "FILE" },
		{ "quiet",       'q', POPT_ARG_NONE,   0,            opt_quiet,   "Don't show progress reports", 0 },
		{ "solexa-scale", 0 , POPT_ARG_NONE,   &solexa_scale,opt_none,    "Quality scores use Solexa formula", 0 },
		{ "fastq-origin", 0 , POPT_ARG_INT,    &fastq_origin,opt_none,    "Quality 0 encodes as ORI, not 33", "ORI" },
		{ "sge-task-last",0 , POPT_ARG_INT,    &stride,      opt_none,    "Override SGE_TASK_LAST env var", "N" },
		POPT_AUTOHELP POPT_TABLEEND
	} ;

	poptContext pc = poptGetContext( "anfo", argc, argv, options, 0 ) ;
	poptSetOtherOptionHelp( pc, "[OPTION...] [sequence-file...]" ) ;
	poptReadDefaultConfig( pc, 0 ) ;

	if( argc <= 1 ) { poptPrintHelp( pc, stderr, 0 ) ; return 1 ; }
	for( int rc = poptGetNextOpt( pc ) ; rc > 0 ; rc = poptGetNextOpt(pc) ) switch( rc )
	{
		case opt_quiet:
			std::clog.rdbuf( 0 ) ;
			break ;

		case opt_version:
			std::cout << poptGetInvocationName(pc) << ", revision " << PACKAGE_VERSION << std::endl ;
			return 0 ;

		default:
			std::clog << poptGetInvocationName(pc) << ": " << poptStrerror( rc ) 
				<< ' ' << poptBadOption( pc, 0 ) << std::endl ;
			return 1 ; 
	}

	if( !output_file ) throw "no output file" ;

	unsigned slicenum = 0, total_slices = 1 ;
	if( const char* tid = getenv("SGE_TASK_ID") ) slicenum = atoi(tid) ? atoi(tid)-1 : 0 ;
	if( const char* tid = getenv("SGE_TASK_LAST") ) if( stride == 1 ) stride = atoi(tid) ? atoi(tid) : 1 ;

	Config conf = get_default_config( config_file ) ;
	for( int i = 0 ; i != conf.policy_size() ; ++i )
	{
		for( int j = 0 ; j != conf.policy(i).use_compact_index_size() ; ++j )
		{
			CompactIndexSpec &ixs = *conf.mutable_policy(i)->mutable_use_compact_index(j) ;
			if( ixs.has_number_of_slices() ) 
			{
				if( ixs.number_of_slices() != total_slices )
				{
					if( total_slices == 1 ) total_slices = ixs.number_of_slices() ;
					else throw "multiple differently sliced indices won't work" ;
				}
				expand_placeholcer( *ixs.mutable_name(), slicenum % total_slices ) ;
			}
		}
	}
	slicenum /= total_slices ;
	stride /= total_slices ;
	Mapper mapper( conf ) ;

	deque<string> files ;
	while( const char* arg = poptGetArg( pc ) ) files.push_back( arg ) ;
	poptFreeContext( pc ) ;
	if( files.empty() ) files.push_back( "-" ) ; 

	streams::AnfoWriter os( output_file ) ;
	output::Header ohdr ;
	*ohdr.mutable_config() = conf ;
	ohdr.set_version( PACKAGE_VERSION ) ;
	if( stride > 1 ) 
	{
		ohdr.set_sge_slicing_stride( stride ) ;
		ohdr.add_sge_slicing_index( slicenum ) ;
	}
	for( const char **arg = argv ; arg != argv+argc ; ++arg ) *ohdr.add_command_line() = *arg ;
	if( const char *jobid = getenv( "SGE_JOB_ID" ) ) ohdr.set_sge_job_id( atoi( jobid ) ) ;
	if( const char *taskid = getenv( "SGE_TASK_ID" ) ) ohdr.set_sge_task_id( atoi( taskid ) ) ;
	os.put_header( ohdr ) ; 

	signal( SIGUSR1, sig_handler ) ;
	signal( SIGUSR2, sig_handler ) ;
	signal( SIGTERM, sig_handler ) ;
	signal( SIGINT, sig_handler ) ;

	for( size_t total_count = 0 ; !exit_with && !files.empty() ; files.pop_front() )
	{
		std::auto_ptr< streams::Stream > inp(
			streams::make_input_stream( files.front().c_str(), solexa_scale, fastq_origin ) ) ;

		for( ; !exit_with && inp->get_state() == streams::Stream::have_output ; ++total_count )
		{
			output::Result r = inp->fetch_result() ;
			if( total_count % stride == slicenum ) {
				QSequence ps ;
				std::deque< alignment_type > ol ;
				int pmax = mapper.index_sequence( r, ps, ol ) ;
				if( pmax != INT_MAX ) mapper.process_sequence( ps, pmax, ol, r ) ;
				os.put_result( r ) ;
			}
		}
	}

	output::Footer ofoot ;
	ofoot.set_exit_code( exit_with ) ;
	os.put_footer( ofoot ) ;
	return 0 ;
}

