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

metaindex::Policy select_policy( const metaindex::Config &c, const PreparedSequence &ps )
{
	metaindex::Policy p ;
	for( int i = 0 ; i != c.policy_size() ; ++i )
	{
		const metaindex::Policy &pi = c.policy(i) ;
		if( ( !pi.has_minlength() || pi.minlength() >= ps.length() ) &&
			( !pi.has_maxlength() || pi.maxlength() <= ps.length() ) )
			p.MergeFrom( pi ) ;
	}
	return p ;
}

int main_( int argc, const char * argv[] )
{
	metaindex::Config mi ;
	merge_text_config( "index.txt", mi ) ;

	std::map< std::string, CompactGenome > genomes ;
	std::map< std::string, std::map< int, FixedIndex > > indices ;
	if( !mi.policy_size() ) throw "no policies---nothing to do." ;

	for( int i = 0 ; i != mi.policy_size() ; ++i )
	{
		for( int j = 0 ; j != mi.policy(i).use_compact_index_size() ; ++j )
		{
			const metaindex::CompactIndexSpec &ixs = mi.policy(i).use_compact_index(j) ;
			CompactGenome &g = genomes[ ixs.genome_name() ] ;
			if( !g.get_base() ) {
				const metaindex::Genome &gdef = find_genome( mi, ixs.genome_name() ) ;
				CompactGenome( gdef.filename().c_str() ).swap( g ) ;
			}
			FixedIndex &ix = indices[ ixs.genome_name() ][ ixs.wordsize() ] ;
			if( !ix ) {
				const metaindex::CompactIndex& cix = find_compact_index( mi, ixs.genome_name(), ixs.wordsize() ) ;
				FixedIndex( cix.filename().c_str(), cix.wordsize() ).swap( ix ) ;
			}
		}
	}



	// cout << "Found index " << cix.filename() << " with wordsize " << cix.wordsize() << " and " ;
	// if( cix.has_cutoff() ) cout << "cutoff " << cix.cutoff() << '.' ;
	                  // else cout << "no cutoff." ;
	// cout << "\nFound genome " << gdef.filename() << " of length " << gdef.total_size() << endl ;

	// CompactGenome g( gdef.filename().c_str() ) ;

	int config_out = throw_errno_if_minus1( creat( "output.tab", 0644 ), "writing output" ) ;
	google::protobuf::io::FileOutputStream fos( config_out ) ;
	fos.SetCloseOnDelete( true ) ;

	std::ifstream inp( "input.fa" ) ;
	while( inp )
	{
		std::string name, seq, line ;
		while( inp && inp.peek() != '>' ) inp.ignore( INT_MAX, '\n' ) ;
		inp.get() ;
		inp >> name ;
		inp.ignore( INT_MAX, '\n' ) ;

		while( inp && inp.peek() != '>' ) 
		{
			getline( inp, line ) ;
			seq.append( line ) ;
		}

		PreparedSequence ps( seq.c_str() ) ;
		const metaindex::Policy&      p    = select_policy( mi, ps ) ;
		// const metaindex::CompactInde& cix  = find_compact_index( mi, mi.compact_index(0) ;
		// const metaindex::Genome&      gdef = find_genome( mi, cix.genome_name() ) ;

		deque<flat_alignment> ol ;
		for( int i = 0 ; i != p.use_compact_index_size() ; ++i )
		{
			const metaindex::CompactIndexSpec &cis = p.use_compact_index(i) ;

			vector<Seed> seeds ;
			int num_raw = indices[ cis.genome_name() ][ cis.wordsize() ].lookup( seq, seeds,
					cis.has_cutoff() ? cis.cutoff() : std::numeric_limits<uint32_t>::max() ) ;
			int num_comb = seeds.size() ;
			select_seeds( seeds, p.max_diag_skew(), p.max_gap(), p.min_seed_len() ) ;
			int num_clumps = seeds.size() ;

			cout << name << ": got " << num_raw << " seeds, combined into "
				 << num_comb << " larger ones, clumped into " << num_clumps
				 << " clumps." << endl ;

			setup_alignments( genomes[ cis.genome_name() ], ps, seeds.begin(), seeds.end(), ol ) ;
		}
		
		if( ol.empty() ) continue ; // XXX
		if( p.has_repeat_threshold() && ol.size() >= p.repeat_threshold() ) continue ; // XXX

		flat_alignment best = find_cheapest( ol, p.max_penalty_per_nuc() * ps.length() / 1000 ) ;
		if( !best ) continue ; // XXX
		int penalty = best.penalty ;

		// cout << "Done near " << best.reference - g.get_base() << " costing " << best.penalty << ':' << endl ;

		deque< pair< flat_alignment, const flat_alignment* > > ol_ ;
		reset( best ) ;
		greedy( best ) ;
		(enter_bt<flat_alignment>( ol_ ))( best ) ;
		Trace t = find_cheapest( ol_ ) ;

		output::Result r ;
		r.set_seqid( name ) ;
		output::Hit *h = r.mutable_best_hit() ;
		// XXX: h->set_refseq( gdef.name() ) ;
		// XXX: h->set_offset( best.reference  + best.ref_offs - g.get_base() ) ;
		h->set_score( penalty ) ;
		for( Trace::const_iterator i = t.begin(), e = t.end() ; i != e ; ++i )
		{
			h->mutable_ref()->push_back( from_ambicode( i->first ) ) ;
			h->mutable_qry()->push_back( from_ambicode( i->second ) ) ;
		}

		google::protobuf::TextFormat::Print( r, &fos ) ;

		// cout <<  << endl ;
	}
	return 0 ;
}

