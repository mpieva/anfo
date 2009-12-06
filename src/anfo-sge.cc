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

//! \page anfo_for_grid_engine ANFO executable for the SGE
//! This is the version of ANFO hacked to run smoothly on hostile
//! environments, e.g. our broken setup of the Sun Grid Engine.  These
//! are the concessions to the environment:
//!
//! * This version is single threaded (doesn't even link pthreads).
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

#include <cstdio>
#include <csignal>
#include <cstring>
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
	
WRAPPED_MAIN
{
	GOOGLE_PROTOBUF_VERIFY_VERSION ;
	enum option_tags { opt_none, opt_version } ;

	const char* config_file = 0 ;
	const char* output_file = 0 ; 
	int solexa_scale = 0 ;
	int clobber = 0 ;
	int fastq_origin = 33 ;
    int task_id = 0 ;
	if( const char *t = getenv( "SGE_TASK_ID" ) ) task_id = atoi( t ) -1 ; 

	struct poptOption options[] = {
		{ "version",     'V', POPT_ARG_NONE,   0,            opt_version, "Print version number and exit", 0 },
		{ "config",      'c', POPT_ARG_STRING, &config_file, opt_none,    "Read config from FILE", "FILE" },
		{ "output",      'o', POPT_ARG_STRING, &output_file, opt_none,    "Write output to FILE", "FILE" },
		{ "clobber",     'C', POPT_ARG_NONE,   &clobber,     opt_none,    "Overwrite existing output file", 0 },
		{ "nommap",       0 , POPT_ARG_NONE,&Metagenome::nommap,opt_none, "Don't use mmap(), read() indexes instead", 0 },
		{ "solexa-scale", 0 , POPT_ARG_NONE,   &solexa_scale,opt_none,    "Quality scores use Solexa formula", 0 },
		{ "fastq-origin", 0 , POPT_ARG_INT,    &fastq_origin,opt_none,    "Quality 0 encodes as ORI, not 33", "ORI" },
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

	if( !output_file ) throw "no output file" ;

    std::string of( expand( output_file, task_id ) ) ;

	// no-op if output exists and overwriting wasn't asked for
    if( !clobber && 0 == access( of.c_str(), F_OK ) ) return 0 ;

    console.set_quiet() ;
    of.append( ".#new#" ) ;
	streams::ChunkedWriter os( of.c_str(), 25 ) ; // prefer speed over compression

	Config conf = get_default_config( config_file ) ;
	Mapper mapper( conf ) ;

	deque<string> files ;
	while( const char* arg = poptGetArg( pc ) ) files.push_back( expand( arg, task_id ) ) ;
	poptFreeContext( pc ) ;
	if( files.empty() ) files.push_back( "-" ) ; 

	output::Header ohdr ;
	*ohdr.mutable_config() = conf ;
	ohdr.set_version( PACKAGE_VERSION ) ;
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
		streams::StreamHolder inp(
				streams::make_input_stream( files.front().c_str(), solexa_scale, fastq_origin ) ) ;

		for( ; !exit_with && inp->get_state() == streams::Stream::have_output ; ++total_count )
		{
			output::Result r = inp->fetch_result() ;
			QSequence ps ;
			std::deque< alignment_type > ol ;
			int pmax = mapper.index_sequence( r, ps, ol ) ;
			if( pmax != INT_MAX ) mapper.process_sequence( ps, pmax, ol, r ) ;
			os.put_result( r ) ;
		}
	}

	output::Footer ofoot ;
	ofoot.set_exit_code( exit_with ) ;
	os.put_footer( ofoot ) ;

	if( !exit_with ) std::rename( of.c_str(), expand( output_file, task_id ).c_str() ) ;
	return 0 ;
}

