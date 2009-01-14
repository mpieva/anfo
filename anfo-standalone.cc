#include "align.h"
#include "conffile.h"
#include "index.h"
#include "outputfile.h"
#include "util.h"

#include "output.pb.h"

#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/text_format.h>

#include <fstream>

#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>

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

config::Policy select_policy( const config::Config &c, const QSequence &ps )
{
	config::Policy p ;
	for( int i = 0 ; i != c.policy_size() ; ++i )
	{
		const config::Policy &pi = c.policy(i) ;
		if( ( !pi.has_minlength() || pi.minlength() <= ps.length() ) &&
			( !pi.has_maxlength() || pi.maxlength() >= ps.length() ) )
			p.MergeFrom( pi ) ;
	}
	return p ;
}

int main_( int argc, const char * argv[] )
{
	config::Config mi = parse_text_config( "config.txt" ) ;

	typedef std::map< std::string, CompactGenome > Genomes ; 
	typedef std::map< std::string, FixedIndex > Indices ;
	Genomes genomes ;
	Indices indices ;
	if( !mi.policy_size() ) throw "no policies---nothing to do." ;

	output::Header ohd ;
	for( int i = 0 ; i != mi.policy_size() ; ++i )
	{
		for( int j = 0 ; j != mi.policy(i).use_compact_index_size() ; ++j )
		{
			const config::CompactIndexSpec &ixs = mi.policy(i).use_compact_index(j) ;
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

	std::ofstream config_out( "output.tab" ) ;
	google::protobuf::io::OstreamOutputStream oos( &config_out ) ;
	google::protobuf::io::CodedOutputStream cos( &oos ) ;
	write_delimited_message( cos, ohd ) ;

	std::ifstream inp( "input.fa" ) ;
	QSequence ps ;
	while( read_fastq( inp, ps ) )
	{
		output::Result r ;
		std::clog << '\r' << ps.get_name() << "\33[K" << std::flush ;

		r.set_seqid( ps.get_name() ) ;
		if( !ps.get_descr().empty() ) r.set_description( ps.get_descr() ) ;
		r.set_sequence( ps.as_string() ) ;
		// XXX: set trim points?

		const config::Policy& p = select_policy( mi, ps ) ;

		deque<flat_alignment> ol ;
		int num_raw = 0, num_comb = 0, num_clumps = 0 ;
		for( int i = 0 ; i != p.use_compact_index_size() ; ++i )
		{
			const config::CompactIndexSpec &cis = p.use_compact_index(i) ;
			const FixedIndex &ix = indices[ cis.name() ] ;
			const CompactGenome &g = genomes[ ix.ci_.genome_name() ] ;
			assert( ix ) ; assert( g.get_base() ) ;

			vector<Seed> seeds ;
			num_raw += ix.lookup( ps, seeds, cis.has_cutoff() ? cis.cutoff() : std::numeric_limits<uint32_t>::max() ) ;
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
					if( const config::Sequence *sequ = g->second.translate_back( t.minpos+1, start_pos ) )
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
	std::clog << std::endl ;
	return 0 ;
}

