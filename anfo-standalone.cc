#include "config.h"

#include "align.h"
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
#include <csignal>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <string>

#include <fnmatch.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>

//! \brief Configures type of alignment to be done.
//! This is a hack and only intended as a stopgap.
// typedef flat_alignment alignment_type ;
typedef simple_adna alignment_type ;

using namespace config ;
using namespace std ;
using namespace google::protobuf::io ;

volatile int exit_with = 0 ;
extern "C" RETSIGTYPE sig_handler( int sig ) { exit_with = sig + 128 ; }
	
//! \page anfo_executable Standalone ANFO executable
//! This is work in progress; it may morph into an ANFO executable to be
//! run directly from the command line.  Right now it reads a FASTA or
//! FASTQ file and maps it against whatever index is configured.
//!
//! What's not obvious from the command line help:  ANFO can run in
//! multiple threads, if you have an SMP machine, that is definitely
//! what you want, but even a single CPU machine can benefit from
//! multithreading as long as not the whole index is cached in memory.
//! The exception is that some hostile environments (Sun Grid Engine)
//! tend to frown upon multithreaded programs with nontrivial memory
//! consumption.  Therefore, the -p (--threads) option can request any
//! number of worker threads, and you get two IO threads in addition.
//! Setting -p 1 gives you one worker and two IOs, but as a special
//! case, -p 0 turns off multithreading altogether.
//!
//! \todo We want an E-value...
//! \todo We want more than just the best match.  Think about a sensible
//!       way to configure this.
//! \todo Test this: the canonical test case is homo sapiens, chr 21.
//! \todo Memory management and pointer/reference conventions are
//!       somewhat wonky in here.  Deserves a thourough audit.

Policy select_policy( const Config &c, const QSequence &ps )
{
	Policy p ;
	for( int i = 0 ; i != c.policy_size() ; ++i )
	{
		const Policy &pi = c.policy(i) ;
		if( ( !pi.has_minlength() || pi.minlength() <= ps.length() ) &&
			( !pi.has_maxlength() || pi.maxlength() >= ps.length() ) &&
			( !pi.has_name_pattern() || 0 == fnmatch( pi.name_pattern().c_str(), ps.get_name().c_str(), 0 ) ) )
			p.MergeFrom( pi ) ;
	}
	return p ;
}

/*! \brief tries to set proc title for ps
 *
 * This is borderline malpractice: since Linux is missing the
 * appropriate API, we directly overwrite argv[0].  This \e will trash
 * argv, so don't try to access that after setting a title, and this \e
 * may trash other things, should my assumption that argv is terminated
 * by "\0\0" turn out wrong.
 *
 * Ignoring the above, it is quite helpful, though...
 * \param title new program title to be displayed
 */

void set_proc_title( const char *title ) 
{
	extern char* __progname_full ;
	extern char* __progname ;
	static char* pe = 0 ;
	static char* pa = 0 ;

	if( !pe ) {
		char* p = __progname ;
		pa = __progname_full ;
		for( pe = pa ; pe[0] || pe[1] ; ++pe ) ;
		while( *p && pa != pe ) *pa++ = *p++ ;
		if( pa != pe ) *pa++ = ':' ;
		if( pa != pe ) *pa++ = ' ' ;
	}

	char* pf = pa ;
	if( pf != pe ) {
		while( *title && pf != pe ) *pf++ = *title++ ;
		while( pf != pe ) *pf++ = 0 ;
	}
}

typedef map< string, CompactGenome > Genomes ; 
typedef map< string, FixedIndex > Indices ;

struct AlignmentWorkload
{
	int pmax ;
	std::deque< alignment_type > *ol ;
	output::Result *r ;
	QSequence *ps ;
} ;

struct CommonData
{
	Queue<output::Result*, 4> output_queue ;
	Queue<AlignmentWorkload*, 4> intermed_queue ;
	Queue<QSequence*, 4> input_queue ;
	google::protobuf::io::CodedOutputStream *output_stream ;
	Config mi ;
	Genomes genomes ;
	Indices indices ;
} ;

void* run_output_thread( void* p )
{
	CommonData *q = (CommonData*)p ;
	while( output::Result *r = q->output_queue.dequeue() )
	{
		write_delimited_message( *q->output_stream, 2, *r ) ;
		delete r ;
	}
	return 0 ;
}

struct reference_overlaps {
	DnaP x, y ;
	reference_overlaps( DnaP u, DnaP v ) : x(u), y(v) {}
	bool operator()( const alignment_type& a ) {
		return a.reference >= x && a.reference <= y ; }
} ;

int index_sequence( CommonData *cd, QSequence &ps, output::Result *r, std::deque< alignment_type >& ol )
{
	r->set_seqid( ps.get_name() ) ;
	if( !ps.get_descr().empty() ) r->set_description( ps.get_descr() ) ;
	switch( ps.get_validity() ) 
	{
		case QSequence::bases_with_quality:
			r->set_quality( ps.qualities() ) ;

		case QSequence::bases_only:
			r->set_sequence( ps.as_string() ) ;
			break ;

		case QSequence::bases_with_qualities:
			r->set_sequence( ps.as_string() ) ;

		case QSequence::qualities_only:
			std::string *qs[] = { r->mutable_four_quality()->mutable_quality_a(),
				                  r->mutable_four_quality()->mutable_quality_c(),
				                  r->mutable_four_quality()->mutable_quality_t(),
				                  r->mutable_four_quality()->mutable_quality_g() } ;
			ps.four_qualities( qs ) ;
	}

	// XXX: set trim points?

	Policy p = select_policy( cd->mi, ps ) ;

	int num_raw = 0, num_comb = 0, num_clumps = 0 ;
	for( int i = 0 ; i != p.use_compact_index_size() ; ++i )
	{
		const CompactIndexSpec &cis = p.use_compact_index(i) ;
		const FixedIndex &ix = cd->indices[ cis.name() ] ;
		const CompactGenome &g = cd->genomes[ ix.ci_.genome_name() ] ;
		assert( ix ) ; assert( g.get_base() ) ;

		vector<Seed> seeds ;
		num_raw += ix.lookup( ps, seeds, cis.has_cutoff() ? cis.cutoff() : numeric_limits<uint32_t>::max() ) ;
		num_comb += seeds.size() ;
		select_seeds( seeds, p.max_diag_skew(), p.max_gap(), p.min_seed_len(), g.get_contig_map() ) ;
		num_clumps += seeds.size() ;

		setup_alignments( g, ps, seeds.begin(), seeds.end(), ol ) ;
	}
	r->set_num_raw_seeds( num_raw ) ;
	r->set_num_grown_seeds( num_comb ) ;
	r->set_num_clumps( num_clumps ) ;

	if( !p.has_max_penalty_per_nuc() )
	{
		r->set_reason( output::no_policy ) ;
		return INT_MAX ;
	}
	else if( ol.empty() ) 
	{
		r->set_reason( output::no_seeds ) ;
		return INT_MAX ;
	}
	else if( p.has_repeat_threshold() && ol.size() >= p.repeat_threshold() )
	{
		r->set_reason( output::too_many_seeds ) ;
		return INT_MAX ;
	}
	else return p.max_penalty_per_nuc() ;
}

void process_sequence( CommonData *cd, QSequence *ps, double max_penalty_per_nuc, std::deque< alignment_type > &ol, output::Result *r )
{
	uint32_t o, c, t, max_penalty = (uint32_t)( max_penalty_per_nuc * ps->length() ) ;
	alignment_type::ClosedSet cl ;
	alignment_type best = find_cheapest( ol, cl, max_penalty, &o, &c, &t ) ;
	r->set_open_nodes_after_alignment( o ) ;
	r->set_closed_nodes_after_alignment( c ) ;
	r->set_tracked_closed_nodes_after_alignment( t ) ;
	if( !best )
	{
		r->set_reason( output::bad_alignment ) ;
	}
	else
	{
		int penalty = best.penalty ;

		deque< pair< alignment_type, const alignment_type* > > ol_ ;
		reset( best ) ;
		greedy( best ) ;
		(enter_bt<alignment_type>( ol_ ))( best ) ;
		DnaP minpos, maxpos ;
		std::vector<uint8_t> t = find_cheapest( ol_, minpos, maxpos ) ;

		output::Hit *h = r->mutable_best_to_genome() ;

		for( Genomes::const_iterator g = cd->genomes.begin(), ge = cd->genomes.end() ; g != ge ; ++g )
		{
			uint32_t start_pos ;
			int32_t len = maxpos - minpos - 1 ;
			if( const Sequence *sequ = g->second.translate_back( minpos+1, start_pos ) )
			{
				h->set_genome_file( g->first ) ;
				if( g->second.g_.has_name() ) 
					h->set_genome_name( g->second.g_.name() ) ;

				h->set_sequence( sequ->name() ) ;
				if( sequ->has_taxid() ) h->set_taxid( sequ->taxid() ) ;
				else if( g->second.g_.has_taxid() ) h->set_taxid( g->second.g_.taxid() ) ;

				h->set_start_pos( start_pos ) ;
				h->set_aln_length( minpos.is_reversed() ? -len : len ) ;
				break ;
			}
		}

		h->mutable_cigar()->assign( t.begin(), t.end() ) ;
		h->set_score( penalty ) ;
		// XXX: h->set_evalue

		//! \todo Find second best hit and similar stuff.
		//! We want the distance to the next best hit; also,
		//! unless already found, we want the best hit to some
		//! selected genome(s).

		// XXX set diff_to_next_species, diff_to_next_order
		// XXX find another best hit (genome only)

		// get rid of overlaps of that first alignment, then look
		// for the next one
		// XXX this is cumbersome... need a better PQueue impl...
		// XXX make distance configurable
		ol.erase( 
				std::remove_if( ol.begin(), ol.end(), reference_overlaps( minpos, maxpos ) ),
				ol.end() ) ;
		make_heap( ol.begin(), ol.end() ) ;
		alignment_type second_best = find_cheapest( ol, cl, max_penalty ) ;
		if( second_best ) r->set_diff_to_next( second_best.penalty - penalty ) ;
	}
}

void* run_indexer_thread( void* cd_ )
{
	CommonData *cd = (CommonData*)cd_ ;
	while( QSequence *ps = cd->input_queue.dequeue() )
	{
		output::Result *r = new output::Result ;
		std::deque< alignment_type > *ol = new std::deque< alignment_type >() ;
		int pmax = index_sequence( cd, *ps, r, *ol ) ;
		if( pmax!= INT_MAX ) 
		{
			AlignmentWorkload *w = new AlignmentWorkload ;
			w->pmax = pmax ;
			w->ps = ps ;
			w->ol = ol ;
			w->r = r ;
			cd->intermed_queue.enqueue( w ) ;
		}
		else
		{
			cd->output_queue.enqueue( r ) ;
			delete ps ;
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
		process_sequence( cd, w->ps, w->pmax, *w->ol, w->r ) ;
		cd->output_queue.enqueue( w->r ) ;
		delete w->ol ;
		delete w->ps ;
		delete w ;
	}
	cd->intermed_queue.enqueue(0) ;
	return 0 ;
}

void expand_placeholcer( string &s, int x )
{
	stringstream ss ; ss << x ;
	for( size_t p ; string::npos != (p = s.find( "$$" )) ; )
		s.replace( p, 2, ss.str() ) ;
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
	int stride = 1 ;
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

	CommonData common_data ;
	if( config_file ) common_data.mi = parse_text_config( config_file ) ;
	else if( !access( "anfo.cfg", F_OK ) ) common_data.mi = parse_text_config( "anfo.cfg" ) ;
	else if( !access( ".anfo.cfg", F_OK ) ) common_data.mi = parse_text_config( ".anfo.cfg" ) ;
	else {
		std::string f = getenv( "HOME" ) + std::string( ".anfo.cfg" ) ;
		if( !access( f.c_str(), F_OK ) ) common_data.mi = parse_text_config( f.c_str() ) ;
		else throw "no config file found" ;
	}

	if( common_data.mi.has_aligner() ) simple_adna::configure( common_data.mi.aligner(), log_params ? &cout : 0 ) ;
	if( !common_data.mi.policy_size() ) throw "no policies---nothing to do." ;
	if( !output_file ) throw "no output file" ;
	if( nthreads < 0 ) throw "invalid thread number" ;


	unsigned slicenum = 0, total_slices = 1 ;
	if( const char* tid = getenv("SGE_TASK_ID") ) slicenum = atoi(tid) -1 ;
	if( const char* tid = getenv("SGE_TASK_LAST") ) if( stride == 1 ) stride = atoi(tid) ;

	for( int i = 0 ; i != common_data.mi.policy_size() ; ++i )
	{
		for( int j = 0 ; j != common_data.mi.policy(i).use_compact_index_size() ; ++j )
		{
			CompactIndexSpec &ixs = *common_data.mi.mutable_policy(i)->mutable_use_compact_index(j) ;
			if( ixs.has_number_of_slices() ) 
			{
				if( ixs.number_of_slices() != total_slices )
				{
					if( total_slices == 1 ) total_slices = ixs.number_of_slices() ;
					else throw "multiple differently sliced indices won't work" ;
				}
				expand_placeholcer( *ixs.mutable_name(), slicenum % total_slices ) ;
			}

			FixedIndex &ix = common_data.indices[ ixs.name() ] ;
			if( !ix ) {
				FixedIndex( ixs.name(), common_data.mi, MADV_WILLNEED ).swap( ix ) ;
				const string& genome_name = ix.ci_.genome_name() ; 
				CompactGenome &g = common_data.genomes[ genome_name ] ;
				if( !g.get_base() )
				{
					CompactGenome( genome_name, common_data.mi, MADV_WILLNEED ).swap( g ) ;
				}
			}
		}
	}
	slicenum /= total_slices ;
	stride /= total_slices ;

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

	std::string output_file_ = output_file ;
	output_file_ += ".#new#" ;

	google::protobuf::io::FileOutputStream fos( strcmp( output_file, "-" ) ?
			throw_errno_if_minus1( open( output_file_.c_str(), O_WRONLY | O_CREAT, 0666 ),
				                   "opening", output_file_.c_str() ) : 1 ) ;
	if( strcmp( output_file, "-" ) ) fos.SetCloseOnDelete( true ) ;
	std::auto_ptr< google::protobuf::io::ZeroCopyOutputStream > zos( compress_fast( &fos ) ) ;

	output::Header ohdr ;
	google::protobuf::io::CodedOutputStream cos( zos.get() ) ;
	cos.WriteRaw( "ANFO", 4 ) ; // signature
	*ohdr.mutable_config() = common_data.mi ;
	ohdr.set_version( PACKAGE_VERSION ) ;
	if( stride > 1 ) 
	{
		ohdr.set_sge_slicing_stride( stride ) ;
		ohdr.add_sge_slicing_index( slicenum ) ;
	}
	for( const char **arg = argv ; arg != argv+argc ; ++arg ) *ohdr.add_command_line() = *arg ;
	if( const char *jobid = getenv( "SGE_JOB_ID" ) ) ohdr.set_sge_job_id( atoi( jobid ) ) ;
	if( const char *taskid = getenv( "SGE_TASK_ID" ) ) ohdr.set_sge_task_id( atoi( taskid ) ) ;
	write_delimited_message( cos, 1, ohdr ) ;

	signal( SIGUSR1, sig_handler ) ;
	signal( SIGUSR2, sig_handler ) ;
	signal( SIGTERM, sig_handler ) ;
	signal( SIGINT, sig_handler ) ;

	// Running in multiple threads.  The main thread will read the
	// input and enqueue it, then signal end of input by adding a null
	// pointer.  It will then wait for the worker threads, signal end of
	// output, wait for the output process.

	common_data.output_stream = &cos ;
	pthread_t output_thread ;
	if( nthreads ) pthread_create( &output_thread, 0, run_output_thread, &common_data ) ;

	pthread_t worker_thread[ nthreads ] ;
	for( int i = 0 ; i != nthreads ; ++i )
		pthread_create( worker_thread+i, 0, run_worker_thread, &common_data ) ;

	pthread_t indexer_thread[ nxthreads ] ;
	if( nthreads )
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

			if( total_count % stride == slicenum ) {
				if( nthreads ) common_data.input_queue.enqueue( ps.release() ) ;
				else 
				{
					output::Result r ;
					std::deque< alignment_type > ol ;
					int pmax = index_sequence( &common_data, *ps, &r, ol ) ;
					if( pmax != INT_MAX ) process_sequence( &common_data, ps.get(), pmax, ol, &r ) ;
					write_delimited_message( *common_data.output_stream, 2, r ) ;
				}
			}
		}
		if( total_size != -1 ) total_done += raw_inp.ByteCount() ;
		if( inp_fd ) close( inp_fd ) ;
	}
	if( nthreads )
	{
		// no more input, wait for indexer(s)
		common_data.input_queue.enqueue( 0 ) ;
		if( nthreads )
			for( int i = 0 ; i != nxthreads ; ++i )
				pthread_join( indexer_thread[i], 0 ) ;

		// done with indexes... wait for the worker(s)
		common_data.intermed_queue.enqueue( 0 ) ;
		for( int i = 0 ; i != nthreads ; ++i )
			pthread_join( worker_thread[i], 0 ) ;

		// end output and wait for it to finish
		common_data.output_queue.enqueue( 0 ) ;
		pthread_join( output_thread, 0 ) ;
	}

	clog << endl ;
	output::Footer ofoot ;
	ofoot.set_exit_code( exit_with ) ;
	write_delimited_message( cos, 3, ofoot ) ;

	if( !exit_with && strcmp( output_file, "-" ) )
		throw_errno_if_minus1(
				rename( output_file_.c_str(), output_file ),
				"renaming output" ) ;
	return 0 ;
}

