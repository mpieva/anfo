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

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include "anfo_common.h"
#include "compress_stream.h"
#include "concurrent_stream.h"
#include "conffile.h"
#include "index.h"
#include "misc_streams.h"
#include "stream.h"
#include "util.h"

#include "output.pb.h"

#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/text_format.h>

#include <popt.h>

#include <algorithm>
#include <cstring>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <string>

#include <fnmatch.h>

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
using namespace streams ;

//! \page anfo_standalone Standalone ANFO executable
//! What's not obvious from the command line help:  ANFO can run in
//! multiple threads.  If you have an SMP machine, that is definitely
//! what you want, but even a single CPU machine can benefit from
//! multithreading while the index is being (implicitly) fetched to
//! memory.  The -p option can request any number of worker
//! threads for aligning and indexing, and you get a main (control and
//! IO) thread in addition.
//!
//! \todo We want an E-value...
//! \todo We want more than just the best match.  Think about a sensible
//!       way to configure this.

//! \page finding_alns How to find everything we need
//! We look for best hits globally and specifically on one genome.  They
//! will be discovered in the order of decreasing score.  After setup,
//! we can operate in a loop of cleaning out the stuff we don't need
//! anymore and finding more alignments.
//!
//! Cleanup:
//! - If we don't have a best hit, everything is needed.
//! - If we don't have a hit to a different species, alignments to any
//!   different species are needed.
//! - If we don't have a hit to a different order, alignments to any
//!   different order are needed.
//! - If we don't have a hit to the human genome, alignments to the
//!   human genome are needed.
//! - If we don't have the second best hit, non-overlapping alignments
//!   to the human genome are needed.
//! - If we didn't hit a different chromosome, alignments to different
//!   chromosomes are needed.
//! - If we didn't hit a different class (autosomes, sex chromosomes,
//!   organelles), those are needed.
//!
//! We're done if nothing is left to align or no alignment is found (or
//! if we got everything, which means everything is thrown out nothing
//! is left).
//!
//! Assignment:
//! For every alignment, just check if it fits anywhere, then store it
//! appropriately.  Expand the two we were interested in.
//!
//! \todo Actually implement search for multiple alignments in its full generality...

WRAPPED_MAIN
{
	GOOGLE_PROTOBUF_VERIFY_VERSION ;
	enum option_tags { opt_none, opt_version, opt_housekeep, opt_index, opt_genome } ;

	std::vector< pair< int, string > > more_opts ;

	const char* config_file = 0 ;
	const char* output_file = 0 ; 
	int nthreads = 1 ;
	int solexa_scale = 0 ;
	int clobber = 0 ;
	int fastq_origin = 33 ;
	int task_id = 0 ;
	if( const char *t = getenv( "SGE_TASK_ID" ) ) task_id = atoi( t ) -1 ; 
	if( const char *t = getenv( "NSLOTS" ) ) nthreads = atoi( t ) ;

	struct poptOption options[] = {
		{ "version",     'V', POPT_ARG_NONE,   0,                   opt_version,  "Print version number and exit", 0 },
		{ "config",      'c', POPT_ARG_STRING, &config_file,        opt_none,     "Read config from FILE", "FILE" },
		{ "output",      'o', POPT_ARG_STRING, &output_file,        opt_none,     "Write output to FILE", "FILE" },
		{ "housekeeping",'H', POPT_ARG_NONE,   0,                   opt_housekeep,"Perform housekeeping (e.g. trimming)", 0 },
		{ "index",       'i', POPT_ARG_STRING, 0,			        opt_index,    "Use FILE as index according to config", "FILE" },
		{ "genome",      'g', POPT_ARG_STRING, 0,			        opt_genome,   "Use FILE as genome according to config", "FILE" },
		{ "threads",     'p', POPT_ARG_INT,    &nthreads,           opt_none,     "Run next step in N parallel worker threads", "N" },
		{ "clobber",     'C', POPT_ARG_NONE,   &clobber,            opt_none,     "Overwrite existing output file", 0 },
		{ "nommap",       0 , POPT_ARG_NONE,   &Metagenome::nommap, opt_none,     "Don't use mmap(), read() indexes instead", 0 },
		{ "solexa-scale", 0 , POPT_ARG_NONE,   &solexa_scale,       opt_none,     "Quality scores use Solexa formula", 0 },
		{ "fastq-origin", 0 , POPT_ARG_INT,    &fastq_origin,       opt_none,     "Quality 0 encodes as ORI, not 33", "ORI" },
		{ "quiet",       'q', POPT_ARG_VAL,    &console.loglevel, Console::error, "Suppress most output", 0 },
		{ "verbose",     'v', POPT_ARG_VAL,    &console.loglevel, Console::info,  "Produce more output", 0 },
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

		case opt_housekeep:
			more_opts.push_back( make_pair( rc + (nthreads << 4), string() ) ) ;
			break ;

		case opt_index:
		case opt_genome:
			more_opts.push_back( make_pair( rc + (nthreads << 4), poptGetOptArg( pc ) ) ) ;
			break ;

		default:
			std::clog << poptGetInvocationName(pc) << ": " << poptStrerror( rc ) 
				<< ' ' << poptBadOption( pc, 0 ) << std::endl ;
			return 1 ; 
	}

	if( !output_file ) throw "no output file" ;
	if( nthreads <= 0 ) throw "invalid thread number" ;

	// no-op if output exists and overwriting wasn't asked for
    if( !clobber && 0 == access( output_file, F_OK ) ) return 0 ;

	Config conf = get_default_config( config_file ) ;

	Holder< StreamBundle > comp = new Compose() ;
	{
		Holder< StreamBundle > ins = new ConcatStream() ;
		if( poptPeekArg( pc ) ) 
		{
			while( const char* arg = poptGetArg( pc ) )
				ins->add_stream( new UniversalReader(
							arg, 0, solexa_scale, fastq_origin ) ) ;
		}
		else
		{
			ins->add_stream( new UniversalReader(
						"-", 0, solexa_scale, fastq_origin ) ) ;
		}
		comp->add_stream( ins ) ;
	}
	poptFreeContext( pc ) ;

	for( size_t i = 0 ; i != more_opts.size() ; ++i )
    {
		if( (more_opts[i].first & 0xf) == opt_housekeep )
		{
            comp->add_stream( new Housekeeper( conf ) ) ;
		}
		else if( (more_opts[i].first & 0xf)  == opt_genome )
        {
            if( !conf.has_aligner() ) throw "no aligner configuration---cannot start." ;

			if( more_opts[i].first > 0x1F ) {
				vector< StreamHolder > v ;
				for( int k = 0 ; k != (more_opts[i].first >> 4) ; ++k )
					v.push_back( new Mapper( conf.aligner(), more_opts[i].second ) ) ;
				comp->add_stream( new ConcurrentStream( v.begin(), v.end() ) ) ;
			}
			else
				comp->add_stream( new Mapper( conf.aligner(), more_opts[i].second ) ) ;
        }
        else
        {
			if( more_opts[i].first > 0x1F ) {
				vector< StreamHolder > v ;
				for( int k = 0 ; k != (more_opts[i].first >> 4) ; ++k )
					v.push_back( new Indexer( conf, more_opts[i].second ) ) ;
				comp->add_stream( new ConcurrentStream( v.begin(), v.end() ) ) ;
			}
			else
				comp->add_stream( new Indexer( conf, more_opts[i].second ) ) ;
        }
    }

	{
		StreamHolder outs = new ChunkedWriter( output_file, 25 ) ; // prefer speed over compression

		output::Header ohdr ;
		*ohdr.mutable_config() = conf ;
		ohdr.set_version( PACKAGE_VERSION ) ;
		for( const char **arg = argv ; arg != argv+argc ; ++arg ) *ohdr.add_command_line() = *arg ;
		if( const char *jobid = getenv( "SGE_JOB_ID" ) ) ohdr.set_sge_job_id( atoi( jobid ) ) ;
		if( const char *taskid = getenv( "SGE_TASK_ID" ) ) ohdr.set_sge_task_id( atoi( taskid ) ) ;
		comp->put_header( ohdr ) ; 

		outs->put_header( comp->fetch_header() ) ;
		while( !exit_with && comp->get_state() == Stream::have_output && outs->get_state() == Stream::need_input )
			outs->put_result( comp->fetch_result() ) ;

		output::Footer ofoot = comp->fetch_footer() ;
		exit_with |= ofoot.exit_code() ;
		ofoot.set_exit_code( exit_with ) ;
		outs->put_footer( ofoot ) ;
	}
	return exit_with ;
}

