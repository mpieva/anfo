#include "align.h"
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
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <string>

#include <unistd.h>
#include <sys/stat.h>

using namespace config ;
using namespace std ;

//! \page anfo_executable Standalone ANFO executable
//! This is work in progress; it may morph into an ANFO executable to be
//! run directly from the command line.  Right now it reads a FASTA or
//! FASTQ file and maps it against whatever index is configured.
//!
//! \todo We want an E-value...
//! \todo Commandline is missing, all this is very inflexible.
//! \todo We want more than just the best match.  Think about a sensible
//!       way to configure this.
//! \todo Test this: the canonical test case is homo sapiens, chr 21.

Policy select_policy( const Config &c, const QSequence &ps )
{
	Policy p ;
	for( int i = 0 ; i != c.policy_size() ; ++i )
	{
		const Policy &pi = c.policy(i) ;
		if( ( !pi.has_minlength() || pi.minlength() <= ps.length() ) &&
			( !pi.has_maxlength() || pi.maxlength() >= ps.length() ) )
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
 * \param title new rpogram title to be displayed
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

struct CommonData
{
	Queue<output::Result*, 10> output_queue ;
	Queue<QSequence*, 10> input_queue ;
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
		write_delimited_message( *q->output_stream, *r ) ;
		delete r ;
	}
	return 0 ;
}

struct reference_overlaps {
	DnaP x, y ;
	reference_overlaps( DnaP u, DnaP v ) : x(u), y(v) {}
	bool operator()( const flat_alignment& a ) {
		return a.reference >= x && a.reference <= y ; }
} ;

void* run_worker_thread( void* cd_ )
{
	CommonData *cd = (CommonData*)cd_ ;
	while( QSequence *ps = cd->input_queue.dequeue() )
	{
		output::Result *r = new output::Result ;
		r->set_seqid( ps->get_name() ) ;
		if( !ps->get_descr().empty() ) r->set_description( ps->get_descr() ) ;
		r->set_sequence( ps->as_string() ) ;
		// XXX: set trim points?

		Policy p = select_policy( cd->mi, *ps ) ;

		deque<flat_alignment> ol ;
		int num_raw = 0, num_comb = 0, num_clumps = 0 ;
		for( int i = 0 ; i != p.use_compact_index_size() ; ++i )
		{
			const CompactIndexSpec &cis = p.use_compact_index(i) ;
			const FixedIndex &ix = cd->indices[ cis.name() ] ;
			const CompactGenome &g = cd->genomes[ ix.ci_.genome_name() ] ;
			assert( ix ) ; assert( g.get_base() ) ;

			vector<Seed> seeds ;
			num_raw += ix.lookup( *ps, seeds, cis.has_cutoff() ? cis.cutoff() : numeric_limits<uint32_t>::max() ) ;
			num_comb += seeds.size() ;
			select_seeds( seeds, p.max_diag_skew(), p.max_gap(), p.min_seed_len(), g.get_contig_map() ) ;
			num_clumps += seeds.size() ;

			setup_alignments( g, *ps, seeds.begin(), seeds.end(), ol ) ;
		}
		r->set_num_raw_seeds( num_raw ) ;
		r->set_num_grown_seeds( num_comb ) ;
		r->set_num_clumps( num_clumps ) ;

		if( !p.has_max_penalty_per_nuc() )
		{
			r->set_reason( output::no_policy ) ;
		}
		else if( ol.empty() ) 
		{
			r->set_reason( output::no_seeds ) ;
		}
		else if( p.has_repeat_threshold() && ol.size() >= p.repeat_threshold() )
		{
			r->set_reason( output::too_many_seeds ) ;
		}
		else
		{
			flat_alignment best = find_cheapest( ol, p.max_penalty_per_nuc() * ps->length() / 1000 ) ;
			if( !best )
			{
				r->set_reason( output::bad_alignment ) ;
			}
			else
			{
				int penalty = best.penalty ;

				deque< pair< flat_alignment, const flat_alignment* > > ol_ ;
				reset( best ) ;
				greedy( best ) ;
				(enter_bt<flat_alignment>( ol_ ))( best ) ;
				Trace t = find_cheapest( ol_ ) ;

				output::Hit *h = r->mutable_best_to_genome() ;

				for( Genomes::const_iterator g = cd->genomes.begin(), ge = cd->genomes.end() ; g != ge ; ++g )
				{
					uint32_t start_pos ;
					int32_t len = t.maxpos - t.minpos - 1 ;
					if( const Sequence *sequ = g->second.translate_back( t.minpos+1, start_pos ) )
					{
						h->set_genome( g->first ) ;
						h->set_sequence( sequ->name() ) ;
						if( sequ->has_taxid() ) h->set_taxid( sequ->taxid() ) ;
						else if( g->second.g_.has_taxid() ) h->set_taxid( g->second.g_.taxid() ) ;

						if( t.minpos.is_reversed() )
						{
							h->set_start_pos( start_pos - len + 1 ) ;
							h->set_aln_length( -len ) ;
						}
						else
						{
							h->set_start_pos( start_pos ) ;
							h->set_aln_length( len ) ;
						}
						break ;
					}
				}

				for( Trace_::const_iterator i = t.trace.begin(), e = t.trace.end() ; i != e ; ++i )
				{
					h->mutable_ref()->push_back( from_ambicode( i->first ) ) ;
					h->mutable_qry()->push_back( from_ambicode( i->second ) ) ;
					h->mutable_con()->push_back( i->first == i->second ? '*' : ' ' ) ;
				}
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
				// XXX this is cumbersome... need a better PQueue
				// impl...
				// XXX make distance configurable
				ol.erase( 
						std::remove_if( ol.begin(), ol.end(), reference_overlaps( t.minpos, t.maxpos ) ),
						ol.end() ) ;
				make_heap( ol.begin(), ol.end() ) ;
				flat_alignment second_best = find_cheapest( ol, penalty + 10 ) ;
				if( second_best ) r->set_diff_to_next( second_best.penalty - penalty ) ;
			}
		}
		cd->output_queue.enqueue( r ) ;
		delete ps ;
	}
	cd->input_queue.enqueue(0) ;
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

int main_( int argc, const char * argv[] )
{
	GOOGLE_PROTOBUF_VERIFY_VERSION ;
	enum option_tags { opt_none, opt_version, opt_quiet } ;

	const char* config_file = 0 ;
	const char* output_file = 0 ; 
	int nthreads = 1 ;

	struct poptOption options[] = {
		{ "version",     'V', POPT_ARG_NONE,   0,            opt_version, "Print version number and exit", 0 },
		{ "config",      'c', POPT_ARG_STRING, &config_file, opt_none,    "Read config from FILE", "FILE" },
		{ "threads",     'p', POPT_ARG_INT,    &nthreads,    opt_none,    "Run in N parallel threads", "N" },
		{ "output",      'o', POPT_ARG_STRING, &output_file, opt_none,    "Write output to FILE", "FILE" },
		{ "quiet",       'q', POPT_ARG_NONE,   0,            opt_quiet,   "Don't show progress reports", 0 },
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
			std::cout << poptGetInvocationName(pc) << ", revision " << VERSION << std::endl ;
			return 0 ;

		default:
			std::cerr << argv[0] << ": " << poptStrerror( rc ) 
				<< ' ' << poptBadOption( pc, 0 ) << std::endl ;
			return 1 ; 
	}

	if( !output_file ) throw "no output file" ;
	if( nthreads < 1 ) throw "invalid thread number" ;

	CommonData common_data ;
	if( config_file ) common_data.mi = parse_text_config( config_file ) ;
	else if( !access( "anfo.cfg", F_OK ) ) common_data.mi = parse_text_config( "anfo.cfg" ) ;
	else if( !access( ".anfo.cfg", F_OK ) ) common_data.mi = parse_text_config( ".anfo.cfg" ) ;
	else {
		std::string f = getenv( "HOME" ) + std::string( ".anfo.cfg" ) ;
		if( !access( f.c_str(), F_OK ) ) common_data.mi = parse_text_config( f.c_str() ) ;
		else throw "no config file found" ;
	}

	if( !common_data.mi.policy_size() ) throw "no policies---nothing to do." ;

	output::Header ohd ;
	for( int i = 0 ; i != common_data.mi.policy_size() ; ++i )
	{
		for( int j = 0 ; j != common_data.mi.policy(i).use_compact_index_size() ; ++j )
		{
			const CompactIndexSpec &ixs = common_data.mi.policy(i).use_compact_index(j) ;
			FixedIndex &ix = common_data.indices[ ixs.name() ] ;
			if( !ix ) {
				FixedIndex( ixs.name(), common_data.mi ).swap( ix ) ;

				CompactGenome &g = common_data.genomes[ ix.ci_.genome_name() ] ;
				if( !g.get_base() )
				{
					*ohd.add_genome() = ix.ci_.genome_name() ;
					CompactGenome( ix.ci_.genome_name(), common_data.mi ).swap( g ) ;
				}
			}
		}
	}

	deque<string> files ;
	while( const char* arg = poptGetArg( pc ) ) files.push_back( arg ) ;
	if( files.empty() ) files.push_back( "-" ) ; 

	int64_t total_size = 0, total_done = 0 ;
	for( size_t i = 0 ; i != files.size() && total_size != -1 ; ++i )
	{
		struct stat s ;
		bool good = files[i] != "-" && !stat( files[i].c_str(), &s ) && S_ISREG( s.st_mode ) ;
		if( good ) total_size += s.st_size ; else total_size = -1 ;
	}

	ofstream output_file_stream ;
	ostream& output_stream = strcmp( output_file, "-" ) ?
		(output_file_stream.open( output_file ), output_file_stream) : cout ;

	google::protobuf::io::OstreamOutputStream oos( &output_stream ) ;
	google::protobuf::io::CodedOutputStream cos( &oos ) ;
	write_delimited_message( cos, ohd ) ;

	// Running in multiple threads.  The main thread will read the
	// input and enqueue it, then signal end of input by adding a null
	// pointer.  It will then wait for the worker threads, signal end of
	// output, wait for the output process.

	common_data.output_stream = &cos ;
	pthread_t output_thread ;
	pthread_create( &output_thread, 0, run_output_thread, &common_data ) ;

	pthread_t worker_thread[ nthreads ] ;
	for( int i = 0 ; i != nthreads ; ++i )
		pthread_create( worker_thread+i, 0, run_worker_thread, &common_data ) ;

	for( int total_count = 1 ; !files.empty() ; files.pop_front() )
	{
		ifstream input_file_stream ;
		istream *inp = &cin ;
		if( !files.front().empty() && files.front() != "-" )
		{
			input_file_stream.open( files.front().c_str() ) ;
			inp = &input_file_stream ;
		}

		for(;;)
		{
			QSequence *ps = new QSequence ;
			if( !read_fastq( *inp, *ps ) ) {
				delete ps ;
				break ;
			}
			stringstream progress ;
			progress << ps->get_name() << " (#" << total_count ;
			if( total_size != -1 )
				progress << ", " << (total_done + input_file_stream.tellg()) * 100 / total_size << '%' ;
			progress << ')' ;

			clog << '\r' << progress.str() << "\33[K" << flush ;
			set_proc_title( progress.str().c_str() ) ;

			common_data.input_queue.enqueue( ps ) ;
			++total_count ;
		}
		if( total_size != -1 ) total_done += input_file_stream.tellg() ;
	}
	common_data.input_queue.enqueue( 0 ) ;

	// done with input... wait for the workers
	for( int i = 0 ; i != nthreads ; ++i )
		pthread_join( worker_thread[i], 0 ) ;

	// end output
	common_data.output_queue.enqueue( 0 ) ;
	pthread_join( output_thread, 0 ) ;

	clog << endl ;
	return 0 ;
}

