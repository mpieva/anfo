#include "align.h"
#include "conffile.h"
#include "index.h"
#include "outputfile.h"
#include "util.h"

#include "output.pb.h"

#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/text_format.h>

#include <popt.h>

#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <string>

#include <unistd.h>

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
//! \todo Make this run in N threads (on the order of four).
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

int main_( int argc, const char * argv[] )
{
	enum option_tags { opt_none, opt_version, opt_quiet } ;

	const char* config_file = 0 ;
	const char* output_file = 0 ; 

	struct poptOption options[] = {
		{ "version",     'V', POPT_ARG_NONE,   0,            opt_version, "Print version number and exit", 0 },
		{ "config",      'c', POPT_ARG_STRING, &config_file, opt_none,    "Read config from FILE", "FILE" },
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

	Config mi ;
	if( config_file ) mi = parse_text_config( config_file ) ;
	else if( !access( "anfo.cfg", F_OK ) ) mi = parse_text_config( "anfo.cfg" ) ;
	else if( !access( ".anfo.cfg", F_OK ) ) mi = parse_text_config( ".anfo.cfg" ) ;
	else {
		std::string f = getenv( "HOME" ) + std::string( ".anfo.cfg" ) ;
		if( !access( f.c_str(), F_OK ) ) mi = parse_text_config( f.c_str() ) ;
		else throw "no config file found" ;
	}

	typedef map< string, CompactGenome > Genomes ; 
	typedef map< string, FixedIndex > Indices ;
	Genomes genomes ;
	Indices indices ;
	if( !mi.policy_size() ) throw "no policies---nothing to do." ;

	output::Header ohd ;
	for( int i = 0 ; i != mi.policy_size() ; ++i )
	{
		for( int j = 0 ; j != mi.policy(i).use_compact_index_size() ; ++j )
		{
			const CompactIndexSpec &ixs = mi.policy(i).use_compact_index(j) ;
			FixedIndex &ix = indices[ ixs.name() ] ;
			if( !ix ) {
				FixedIndex( ixs.name(), mi ).swap( ix ) ;

				CompactGenome &g = genomes[ ix.ci_.genome_name() ] ;
				if( !g.get_base() )
				{
					*ohd.add_genome() = ix.ci_.genome_name() ;
					CompactGenome( ix.ci_.genome_name(), mi ).swap( g ) ;
				}
			}
		}
	}

	ofstream output_file_stream ;
	ostream& output_stream = strcmp( output_file, "-" ) ?
		(output_file_stream.open( output_file ), output_file_stream) : cout ;

	google::protobuf::io::OstreamOutputStream oos( &output_stream ) ;
	google::protobuf::io::CodedOutputStream cos( &oos ) ;
	write_delimited_message( cos, ohd ) ;

	deque<string> files ;
	while( const char* arg = poptGetArg( pc ) ) files.push_back( arg ) ;
	if( files.empty() ) files.push_back( "-" ) ;

	for( int total_count = 1 ; !files.empty() ; files.pop_front() )
	{
		ifstream input_file_stream ;
		istream& inp = files.front().empty() || files.front() == "-" 
			? cin : (input_file_stream.open( files.front().c_str() ), input_file_stream) ;

		for( QSequence ps ; read_fastq( inp, ps ) ; ++total_count )
		{
			output::Result r ;
			stringstream progress ;
			progress << ps.get_name() << " (#" << total_count << ")" ;
			clog << '\r' << progress.str() << "\33[K" << flush ;
			set_proc_title( progress.str().c_str() ) ;

			r.set_seqid( ps.get_name() ) ;
			if( !ps.get_descr().empty() ) r.set_description( ps.get_descr() ) ;
			r.set_sequence( ps.as_string() ) ;
			// XXX: set trim points?

			Policy p = select_policy( mi, ps ) ;

			deque<flat_alignment> ol ;
			int num_raw = 0, num_comb = 0, num_clumps = 0 ;
			for( int i = 0 ; i != p.use_compact_index_size() ; ++i )
			{
				const CompactIndexSpec &cis = p.use_compact_index(i) ;
				const FixedIndex &ix = indices[ cis.name() ] ;
				const CompactGenome &g = genomes[ ix.ci_.genome_name() ] ;
				assert( ix ) ; assert( g.get_base() ) ;

				vector<Seed> seeds ;
				num_raw += ix.lookup( ps, seeds, cis.has_cutoff() ? cis.cutoff() : numeric_limits<uint32_t>::max() ) ;
				num_comb += seeds.size() ;
				select_seeds( seeds, p.max_diag_skew(), p.max_gap(), p.min_seed_len(), g.get_contig_map() ) ;
				num_clumps += seeds.size() ;

				setup_alignments( g, ps, seeds.begin(), seeds.end(), ol ) ;
			}
			r.set_num_raw_seeds( num_raw ) ;
			r.set_num_grown_seeds( num_comb ) ;
			r.set_num_clumps( num_clumps ) ;

			if( !p.has_max_penalty_per_nuc() )
			{
				r.set_reason( output::no_policy ) ;
			}
			else if( ol.empty() ) 
			{
				r.set_reason( output::no_seeds ) ;
			}
			else if( p.has_repeat_threshold() && ol.size() >= p.repeat_threshold() )
			{
				r.set_reason( output::too_many_seeds ) ;
			}
			else
			{
				flat_alignment best = find_cheapest( ol, p.max_penalty_per_nuc() * ps.length() / 1000 ) ;
				if( !best )
				{
					r.set_reason( output::bad_alignment ) ;
				}
				else
				{
					int penalty = best.penalty ;

					deque< pair< flat_alignment, const flat_alignment* > > ol_ ;
					reset( best ) ;
					greedy( best ) ;
					(enter_bt<flat_alignment>( ol_ ))( best ) ;
					Trace t = find_cheapest( ol_ ) ;

					output::Hit *h = r.mutable_best_hit() ;

					for( Genomes::const_iterator g = genomes.begin(), ge = genomes.end() ; g != ge ; ++g )
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
				}
			}

			write_delimited_message( cos, r ) ;
		}
	}
	clog << endl ;
	return 0 ;
}

