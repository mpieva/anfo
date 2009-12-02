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

#include "align.h"
#include "anfo_common.h"
#include "compress_stream.h"
#include "conffile.h"
#include "index.h"
#include "stream.h"
#include "queue.h"
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

//! \page anfo_standalone Standalone ANFO executable
//! What's not obvious from the command line help:  ANFO can run in
//! multiple threads, if you have an SMP machine, that is definitely
//! what you want, but even a single CPU machine can benefit from
//! multithreading while the index is being fetched to memory.
//! The -p and -x options can request any number of worker threads for
//! aligning and indexing, and you get two IO threads in addition.  That
//! means not requesting multithreading still gives you four threads in
//! total.
//!
//! \todo We want an E-value...
//! \todo We want more than just the best match.  Think about a sensible
//!       way to configure this.

struct AlignmentWorkload
{
	int pmax ;
	std::deque< alignment_type > ol ;
	std::auto_ptr< output::Result > r ;
	std::auto_ptr< QSequence > ps ;
} ;

struct CommonData
{
	Queue<output::Result*, 8> output_queue ;
	Queue<AlignmentWorkload*, 8> intermed_queue ;
	Queue<output::Result*, 8> input_queue ;
	streams::ChunkedWriter output_stream ;
	Mapper mapper ;

	CommonData( const Config& conf, const char* fn )
		: output_stream( fn, 25 ), mapper( conf ) {}
} ;

void* run_output_thread( void* p )
{
	CommonData *q = (CommonData*)p ;
	while( output::Result *r = q->output_queue.dequeue() )
	{
		q->output_stream.put_result( *r ) ;
		delete r ;
	}
	return 0 ;
}

void* run_indexer_thread( void* cd_ )
{
	CommonData *cd = (CommonData*)cd_ ;
	while( output::Result *r = cd->input_queue.dequeue() )
	{
		std::auto_ptr< AlignmentWorkload > w ( new AlignmentWorkload ) ;
		w->r.reset( r ) ;
		w->ps.reset( new QSequence() ) ;
		w->pmax = cd->mapper.index_sequence( *w->r, *w->ps, w->ol ) ;
		if( w->pmax!= INT_MAX ) cd->intermed_queue.enqueue( w.release() ) ;
		else                    cd->output_queue.enqueue( w->r.release() ) ;
	}
	cd->input_queue.enqueue(0) ;
	return 0 ;
}

void* run_worker_thread( void* cd_ )
{
	CommonData *cd = (CommonData*)cd_ ;
	while( AlignmentWorkload *w = cd->intermed_queue.dequeue() )
	{
		cd->mapper.process_sequence( *w->ps, w->pmax, w->ol, *w->r ) ;
		cd->output_queue.enqueue( w->r.release() ) ;
		delete w ;
	}
	cd->intermed_queue.enqueue(0) ;
	return 0 ;
}

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
	enum option_tags { opt_none, opt_version } ;

	const char* config_file = 0 ;
	const char* output_file = 0 ; 
	int nthreads = 1 ;
	int nxthreads = 1 ;
	int solexa_scale = 0 ;
	int clobber = 0 ;
	int fastq_origin = 33 ;
	int task_id = 0 ;
	if( const char *t = getenv( "SGE_TASK_ID" ) ) task_id = atoi( t ) -1 ; 
	if( const char *t = getenv( "NSLOTS" ) ) { nthreads = atoi( t ) ; nxthreads = (3+nthreads) / 4 ; }


	struct poptOption options[] = {
		{ "version",     'V', POPT_ARG_NONE,   0,            opt_version, "Print version number and exit", 0 },
		{ "config",      'c', POPT_ARG_STRING, &config_file, opt_none,    "Read config from FILE", "FILE" },
		{ "threads",     'p', POPT_ARG_INT,    &nthreads,    opt_none,    "Run in N parallel worker threads", "N" },
		{ "ixthreads",   'x', POPT_ARG_INT,    &nxthreads,   opt_none,    "Run in N parallel indexer threads", "N" },
		{ "output",      'o', POPT_ARG_STRING, &output_file, opt_none,    "Write output to FILE", "FILE" },
		{ "clobber",     'C', POPT_ARG_NONE,   &clobber,     opt_none,    "Overwrite existing output file", 0 },
		{ "nommap",       0 , POPT_ARG_NONE,   &Metagenome::nommap, opt_none,     "Don't use mmap(), read() indexes instead", 0 },
		{ "quiet",       'q', POPT_ARG_VAL,    &console.loglevel, Console::error, "Suppress most output", 0 },
		{ "verbose",     'v', POPT_ARG_VAL,    &console.loglevel, Console::info,  "Produce more output", 0 },
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
	if( nthreads <= 0 ) throw "invalid thread number" ;

	Config conf = get_default_config( config_file ) ;
	std::string ofile = expand( output_file, task_id ) ;

	// no-op if output exists and overwriting wasn't asked for
    if( !clobber && 0 == access( ofile.c_str(), F_OK ) ) return 0 ;

	CommonData common_data( conf, (ofile+".#new#").c_str() ) ;

	deque<string> files ;
	while( const char* arg = poptGetArg( pc ) ) files.push_back( expand( arg, task_id ) ) ;

	poptFreeContext( pc ) ;
	if( files.empty() ) files.push_back( "-" ) ; 

	output::Header ohdr ;
	*ohdr.mutable_config() = conf ;
	ohdr.set_version( PACKAGE_VERSION ) ;

	for( const char **arg = argv ; arg != argv+argc ; ++arg ) *ohdr.add_command_line() = *arg ;
	common_data.output_stream.put_header( ohdr ) ;

	// Running in multiple threads.  The main thread will read the
	// input and enqueue it, then signal end of input by adding a null
	// pointer.  It will then wait for the worker threads, signal end of
	// output, wait for the output process.

	pthread_t output_thread ;
	pthread_create( &output_thread, 0, run_output_thread, &common_data ) ;

	pthread_t worker_thread[ nthreads ] ;
	for( int i = 0 ; i != nthreads ; ++i )
		pthread_create( worker_thread+i, 0, run_worker_thread, &common_data ) ;

	pthread_t indexer_thread[ nxthreads ] ;
	for( int i = 0 ; i != nxthreads ; ++i )
		pthread_create( indexer_thread+i, 0, run_indexer_thread, &common_data ) ;

	for( size_t total_count = 0 ; !exit_with && !files.empty() ; files.pop_front() )
	{
		Holder< streams::Stream > inp( 
				streams::make_input_stream( files.front().c_str(), solexa_scale, fastq_origin ) ) ;

		for( ; !exit_with && inp->get_state() == streams::Stream::have_output ; ++total_count )
			common_data.input_queue.enqueue( new output::Result( inp->fetch_result() ) ) ;
	}
	
	// no more input, wait for indexer(s)
	common_data.input_queue.enqueue( 0 ) ;
	for( int i = 0 ; i != nxthreads ; ++i )
		pthread_join( indexer_thread[i], 0 ) ;

	// done with indexes... wait for the worker(s)
	common_data.intermed_queue.enqueue( 0 ) ;
	for( int i = 0 ; i != nthreads ; ++i )
		pthread_join( worker_thread[i], 0 ) ;

	// end output and wait for it to finish
	common_data.output_queue.enqueue( 0 ) ;
	pthread_join( output_thread, 0 ) ;

	clog << endl ;
	output::Footer ofoot ;
	ofoot.set_exit_code( exit_with ) ;
	common_data.output_stream.put_footer( ofoot ) ;
	if( !exit_with ) std::rename( (ofile+".#new#").c_str(), ofile.c_str() ) ;
	return 0 ;
}

