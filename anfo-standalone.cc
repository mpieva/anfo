#include "config.h"

#include "align.h"
#include "anfo_common.h"
#include "compress_stream.h"
#include "conffile.h"
#include "index.h"
#include "outputfile.h"
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
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>

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
//! \todo Test this: the canonical test case is homo sapiens, chr 21.
//! \todo Memory management and pointer/reference conventions are
//!       somewhat wonky in here.  Deserves a thourough audit.
//! \todo Investigate whether decompressing (LZF?) a complete index to
//!       memory works better than the mmap used now.

struct AlignmentWorkload
{
	int pmax ;
	std::deque< alignment_type > ol ;
	std::auto_ptr< output::Result > r ;
	std::auto_ptr< QSequence > ps ;
} ;

struct CommonData
{
	Queue<output::Result*, 4> output_queue ;
	Queue<AlignmentWorkload*, 4> intermed_queue ;
	Queue<QSequence*, 4> input_queue ;
	google::protobuf::io::CodedOutputStream output_stream ;
	Mapper mapper ;

	CommonData( const Config& conf, google::protobuf::io::ZeroCopyOutputStream *zos )
		: output_stream( zos ), mapper( conf ) {}
} ;

void* run_output_thread( void* p )
{
	CommonData *q = (CommonData*)p ;
	while( output::Result *r = q->output_queue.dequeue() )
	{
		write_delimited_message( q->output_stream, 2, *r ) ;
		delete r ;
	}
	return 0 ;
}

void* run_indexer_thread( void* cd_ )
{
	CommonData *cd = (CommonData*)cd_ ;
	while( QSequence *ps = cd->input_queue.dequeue() )
	{
		std::auto_ptr< AlignmentWorkload > w ( new AlignmentWorkload ) ;
		w->ps.reset( ps ) ;
		w->r.reset( new output::Result ) ;
		w->pmax = cd->mapper.index_sequence( *w->ps, *w->r, w->ol ) ;
		if( w->pmax!= INT_MAX ) 
		{
			cd->intermed_queue.enqueue( w.release() ) ;
		}
		else
		{
			cd->output_queue.enqueue( w->r.release() ) ;
		}
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

int main_( int argc, const char * argv[] )
{
	GOOGLE_PROTOBUF_VERIFY_VERSION ;
	enum option_tags { opt_none, opt_version, opt_quiet } ;

	const char* config_file = 0 ;
	const char* output_file = 0 ; 
	int nthreads = 1 ;
	int nxthreads = 1 ;
	int solexa_scale = 0 ;
	int log_params = 0 ;
	int fastq_origin = 33 ;

	struct poptOption options[] = {
		{ "version",     'V', POPT_ARG_NONE,   0,            opt_version, "Print version number and exit", 0 },
		{ "config",      'c', POPT_ARG_STRING, &config_file, opt_none,    "Read config from FILE", "FILE" },
		{ "threads",     'p', POPT_ARG_INT,    &nthreads,    opt_none,    "Run in N parallel worker threads", "N" },
		{ "ixthreads",   'x', POPT_ARG_INT,    &nxthreads,   opt_none,    "Run in N parallel indexer threads", "N" },
		{ "output",      'o', POPT_ARG_STRING, &output_file, opt_none,    "Write output to FILE", "FILE" },
		{ "quiet",       'q', POPT_ARG_NONE,   0,            opt_quiet,   "Don't show progress reports", 0 },
		{ "dump-params",  0 , POPT_ARG_NONE,   &log_params,  opt_none,    "Print out alignment paramters", 0 },
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
	if( nthreads <= 0 ) throw "invalid thread number" ;

	std::string output_file_ = output_file ;
	output_file_ += ".#new#" ;

	google::protobuf::io::FileOutputStream fos( strcmp( output_file, "-" ) ?
			throw_errno_if_minus1( open( output_file_.c_str(), O_WRONLY | O_CREAT, 0666 ),
				                   "opening", output_file_.c_str() ) : 1 ) ;
	if( strcmp( output_file, "-" ) ) fos.SetCloseOnDelete( true ) ;
	std::auto_ptr< google::protobuf::io::ZeroCopyOutputStream > zos( compress_fast( &fos ) ) ;

	Config conf = get_default_config( config_file ) ;
	CommonData common_data( conf, zos.get() ) ;

	deque<string> files ;
	while( const char* arg = poptGetArg( pc ) ) files.push_back( arg ) ;
	poptFreeContext( pc ) ;
	if( files.empty() ) files.push_back( "-" ) ; 

	int64_t total_size = 0, total_done = 0 ;
	for( size_t i = 0 ; i != files.size() && total_size != -1 ; ++i )
	{
		struct stat s ;
		bool good = files[i] != "-" && !stat( files[i].c_str(), &s ) && S_ISREG( s.st_mode ) ;
		if( good ) total_size += s.st_size ; else total_size = -1 ;
	}

	output::Header ohdr ;
	common_data.output_stream.WriteRaw( "ANFO", 4 ) ; // signature

	*ohdr.mutable_config() = conf ;
	ohdr.set_version( PACKAGE_VERSION ) ;

	for( const char **arg = argv ; arg != argv+argc ; ++arg ) *ohdr.add_command_line() = *arg ;
	write_delimited_message( common_data.output_stream, 1, ohdr ) ;

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
		int inp_fd = files.front().empty() || files.front() == "-" ? 0 :
			throw_errno_if_minus1( open( files.front().c_str(), O_RDONLY ),
					"opening ", files.front().c_str() ) ;

		FileInputStream raw_inp( inp_fd ) ;
		std::auto_ptr<ZeroCopyInputStream> inp( decompress( &raw_inp ) ) ;

		for(;; ++total_count )
		{
			std::auto_ptr<QSequence> ps( new QSequence ) ;
			if( exit_with || !read_fastq( inp.get(), *ps, solexa_scale, fastq_origin ) ) break ;
			
			stringstream progress ;
			progress << ps->get_name() << " (#" << total_count ;
			if( total_size != -1 ) progress << ", " << (total_done + raw_inp.ByteCount()) * 100 / total_size << '%' ;
			progress << ')' ;

			clog << '\r' << progress.str() << "\33[K" << flush ;
			set_proc_title( progress.str().c_str() ) ;

			common_data.input_queue.enqueue( ps.release() ) ;
		}
		if( total_size != -1 ) total_done += raw_inp.ByteCount() ;
		if( inp_fd ) close( inp_fd ) ;
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
	write_delimited_message( common_data.output_stream, 3, ofoot ) ;

	if( !exit_with && strcmp( output_file, "-" ) )
		throw_errno_if_minus1(
				rename( output_file_.c_str(), output_file ),
				"renaming output" ) ;
	return 0 ;
}

