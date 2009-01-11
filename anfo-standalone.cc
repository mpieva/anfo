#include "align.h"
#include "conffile.h"
#include "index.h"
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

metaindex::Policy select_policy( const metaindex::Config &c, const QSequence &ps )
{
	metaindex::Policy p ;
	for( int i = 0 ; i != c.policy_size() ; ++i )
	{
		const metaindex::Policy &pi = c.policy(i) ;
		if( ( !pi.has_minlength() || pi.minlength() <= ps.length() ) &&
			( !pi.has_maxlength() || pi.maxlength() >= ps.length() ) )
			p.MergeFrom( pi ) ;
	}
	return p ;
}

void write_separator( google::protobuf::io::ZeroCopyOutputStream& s )
{
	void *p ;
	int sz ;
	if( !s.Next( &p, &sz ) || sz < 4 ) throw "write error" ;
	((char*)p)[0] = '\n' ;
	((char*)p)[1] = 0x1e ; // RS 
	((char*)p)[2] = '\n' ;
	((char*)p)[3] = '\n' ;
	s.BackUp( sz - 4 ) ;
}

int main_( int argc, const char * argv[] )
{
	metaindex::Config mi ;
	merge_text_config( "config.txt", mi ) ;

	typedef std::map< std::string, CompactGenome > Genomes ; 
	typedef std::map< std::string, std::map< int, FixedIndex > > Indices ;
	Genomes genomes ;
	Indices indices ;
	if( !mi.policy_size() ) throw "no policies---nothing to do." ;

	output::Header ohd ;
	for( int i = 0 ; i != mi.policy_size() ; ++i )
	{
		for( int j = 0 ; j != mi.policy(i).use_compact_index_size() ; ++j )
		{
			const metaindex::CompactIndexSpec &ixs = mi.policy(i).use_compact_index(j) ;
			CompactGenome &g = genomes[ ixs.genome_name() ] ;
			if( !g.get_base() )
			{
				*ohd.add_genome() = ixs.genome_name() ;
				CompactGenome( find_genome( mi, ixs.genome_name() ) ).swap( g ) ;
			}

			FixedIndex &ix = indices[ ixs.genome_name() ][ ixs.wordsize() ] ;
			if( !ix ) {
				const metaindex::CompactIndex& cix = find_compact_index( mi, ixs.genome_name(), ixs.wordsize() ) ;
				FixedIndex( cix.filename().c_str(), cix.wordsize() ).swap( ix ) ;
			}
		}
	}

	int config_out = throw_errno_if_minus1( creat( "output.tab", 0644 ), "writing output" ) ;
	google::protobuf::io::FileOutputStream fos( config_out ) ;
	fos.SetCloseOnDelete( true ) ;

	google::protobuf::TextFormat::Print( ohd, &fos ) ; write_separator( fos ) ; // XXX

	std::ifstream inp( "input.fa" ) ;
	QSequence ps ;
	while( read_fastq( inp, ps ) )
	{
		output::Result r ;
		r.set_seqid( ps.get_name() ) ;
		if( !ps.get_descr().empty() ) r.set_description( ps.get_descr() ) ;
		r.set_sequence( ps.as_string() ) ;
		// XXX: set trim points?

		const metaindex::Policy& p = select_policy( mi, ps ) ;

		deque<flat_alignment> ol ;
		int total_seeds = 0 ;
		for( int i = 0 ; i != p.use_compact_index_size() ; ++i )
		{
			const metaindex::CompactIndexSpec &cis = p.use_compact_index(i) ;

			vector<Seed> seeds ;
			int num_raw = indices[ cis.genome_name() ][ cis.wordsize() ].lookup( ps, seeds,
					cis.has_cutoff() ? cis.cutoff() : std::numeric_limits<uint32_t>::max() ) ;
			int num_comb = seeds.size() ;
			select_seeds( seeds, p.max_diag_skew(), p.max_gap(), p.min_seed_len(),
					genomes[ cis.genome_name() ].get_contig_map() ) ;
			int num_clumps = seeds.size() ;

			cout << ps.get_name() << ": got " << num_raw << " seeds, combined into "
				 << num_comb << " larger ones, clumped into " << num_clumps
				 << " clumps." << endl ;

			setup_alignments( genomes[ cis.genome_name() ], ps, seeds.begin(), seeds.end(), ol ) ;
			total_seeds += num_raw ;
		}
		
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
					if( g->second.translate_back( t.minpos+1, *h->mutable_sequence(), start_pos ) )
					{
						h->set_genome( g->first ) ;
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
				// XXX: h->set_taxid

				//! \todo Find second best hit and similar stuff.
				//! We want the distance to the next best hit; also,
				//! unless already found, we want the best hit to some
				//! selected genome(s).

				// XXX set diff_to_next_species, diff_to_next_order
				// XXX find another best hit (genome only)
			}
		}

		google::protobuf::TextFormat::Print( r, &fos ) ; write_separator( fos ) ; // XXX
	}
	return 0 ;
}

