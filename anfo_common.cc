#include "anfo_common.h"

#include "trim.h"
#include "outputfile.h"

#include "config.pb.h"
#include "output.pb.h"

//! \todo We want an E-value...
//! \todo We want more than just the best match.  Think about a sensible
//!       way to configure this.
//! \todo Test this: the canonical test case is homo sapiens, chr 21.

using namespace config ;
using namespace output ;
using namespace std; 

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

Mapper::Mapper( const config::Config &config ) : mi(config)
{
	if( mi.has_aligner() ) simple_adna::configure( mi.aligner(), 0 ) ;
	if( !mi.policy_size() ) throw "no policies---nothing to do." ;

	for( int i = 0 ; i != mi.policy_size() ; ++i )
	{
		for( int j = 0 ; j != mi.policy(i).use_compact_index_size() ; ++j )
		{
			CompactIndexSpec &ixs = *mi.mutable_policy(i)->mutable_use_compact_index(j) ;
			FixedIndex &ix = indices[ ixs.name() ] ;
			if( !ix ) {
				FixedIndex( ixs.name(), mi, MADV_WILLNEED ).swap( ix ) ;
				const string& genome_name = ix.ci_.genome_name() ; 
				CompactGenome &g = genomes[ genome_name ] ;
				if( !g.get_base() ) {
					CompactGenome( genome_name, mi, MADV_WILLNEED ).swap( g ) ;
				}
			}
		}
	}
}

static const int maxd = 32 ;
static const int minscore = 4 ;

int Mapper::index_sequence( QSequence &ps, output::Result &r, std::deque< alignment_type >& ol )
{
	r.set_seqid( ps.get_name() ) ;
	if( !ps.get_descr().empty() ) r.set_description( ps.get_descr() ) ;
	switch( ps.get_validity() ) 
	{
		case QSequence::bases_with_quality:
			r.set_quality( ps.qualities() ) ;

		case QSequence::bases_only:
			r.set_sequence( ps.as_string() ) ;
			break ;

		case QSequence::bases_with_qualities:
			r.set_sequence( ps.as_string() ) ;

		case QSequence::qualities_only:
			std::string *qs[] = { r.mutable_four_quality()->mutable_quality_a(),
				                  r.mutable_four_quality()->mutable_quality_c(),
				                  r.mutable_four_quality()->mutable_quality_t(),
				                  r.mutable_four_quality()->mutable_quality_g() } ;
			ps.four_qualities( qs ) ;
	}

	// trim adapters, set trim points
	// How does this work?  We create an overlap alignment, then
	// calculate an alignment score from the number of differences.  The
	// formula equals a match score of 1 and a mismatch score of -3.
	// Gaps score (mat-mis)/2 to allow for the simple algorithm.  If the
	// score is good enough, we trim.

	if( mi.trim_right_size() || mi.trim_left_size() ) ps.drop_trailing_ns() ;

	for( int i = 0 ; i != mi.trim_right_size() ; ++i )
	{
		int ymax, xmax = mi.trim_right(i).size() ;
		int diff = align( ps.rbegin(), ps.rend(), mi.trim_right(i).rbegin(), mi.trim_right(i).rend(),
				          maxd, overlap, 0, &ymax ) ;
		int score = xmax + ymax - 8 * diff ;
		if( diff < maxd && score >= minscore && ymax > 0 )
		{
			r.set_trim_right( ps.length() - ymax ) ;
			ps.trim_right( ps.length() - ymax ) ;
		}
	}

	for( int i = 0 ; i != mi.trim_left_size() ; ++i )
	{
		int ymax, xmax = mi.trim_left(i).size() ;
		int diff = align( ps.begin(), ps.end(), mi.trim_left(i).begin(), mi.trim_left(i).end(),
				          maxd, overlap, 0, &ymax ) ;
		int score = xmax + ymax - 8 * diff ;
		if( diff < maxd && score >= minscore && ymax > 0 )
		{
			r.set_trim_left( r.trim_left() + ymax ) ;
			ps.trim_left( ymax ) ;
		}
	}

	Policy p = select_policy( mi, ps ) ;

	int num_raw = 0, num_comb = 0, num_clumps = 0, num_useless = 0 ;
	for( int i = 0 ; i != p.use_compact_index_size() ; ++i )
	{
		const CompactIndexSpec &cis = p.use_compact_index(i) ;
		const FixedIndex &ix = indices[ cis.name() ] ;
		const CompactGenome &g = genomes[ ix.ci_.genome_name() ] ;
		assert( ix ) ; assert( g.get_base() ) ;

		vector<Seed> seeds ;
		num_raw += ix.lookupS( 
				ps, seeds, cis.allow_near_perfect(), &num_useless,
				cis.has_cutoff() ? cis.cutoff() : numeric_limits<uint32_t>::max() ) ;
		num_comb += seeds.size() ;
		select_seeds( seeds, p.max_diag_skew(), p.max_gap(), p.min_seed_len(), g.get_contig_map() ) ;
		num_clumps += seeds.size() ;

		setup_alignments( g, ps, seeds.begin(), seeds.end(), ol ) ;
	}
	r.set_num_raw_seeds( num_raw ) ;
	r.set_num_useless( num_useless ) ;
	r.set_num_grown_seeds( num_comb ) ;
	r.set_num_clumps( num_clumps ) ;

	if( !p.has_max_penalty_per_nuc() )
	{
		r.set_reason( output::no_policy ) ;
		return INT_MAX ;
	}
	else if( ol.empty() ) 
	{
		r.set_reason( num_useless ? output::repeats_only : output::no_seeds ) ;
		return INT_MAX ;
	}
	else if( p.has_repeat_threshold() && ol.size() >= p.repeat_threshold() )
	{
		r.set_reason( output::too_many_seeds ) ;
		return INT_MAX ;
	}
	else return p.max_penalty_per_nuc() ;
}

bool Mapper::translate_to_genome_coords( DnaP pos, uint32_t &xpos, const config::Sequence** s_out, const config::Genome** g_out, std::string* g_file )
{
	for( Genomes::const_iterator g = genomes.begin(), ge = genomes.end() ; g != ge ; ++g )
	{
		if( const Sequence *sequ = g->second.translate_back( pos, xpos ) )
		{
			if( s_out ) *s_out = sequ ;
			if( g_file ) *g_file = g->first ;
			if( g_out) *g_out = &g->second.g_ ;
			return true ;
		}
	}
	return false ;
}

void Mapper::process_sequence( const QSequence &ps, double max_penalty_per_nuc, std::deque< alignment_type > &ol, output::Result &r )
{
	uint32_t o, c, t, max_penalty = (uint32_t)( max_penalty_per_nuc * ps.length() ) ;
	alignment_type::ClosedSet cl ;
	alignment_type best = find_cheapest( ol, cl, max_penalty, &o, &c, &t ) ;
	r.set_open_nodes_after_alignment( o ) ;
	r.set_closed_nodes_after_alignment( c ) ;
	r.set_tracked_closed_nodes_after_alignment( t ) ;
	if( !best )
	{
		r.set_reason( output::bad_alignment ) ;
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
		int32_t len = maxpos - minpos - 1 ;

		output::Hit *h = r.mutable_best_to_genome() ;

		uint32_t start_pos ;
		std::string genome_file ;
		const Sequence *sequ ;
		const Genome *genome ;

		if( translate_to_genome_coords( minpos+1, start_pos, &sequ, &genome, &genome_file ) )
		{
			h->set_genome_file( genome_file ) ;
			if( genome->has_name() ) h->set_genome_name( genome->name() ) ;
			h->set_sequence( sequ->name() ) ;
			if( sequ->has_taxid() ) h->set_taxid( sequ->taxid() ) ;
			else if( genome->has_taxid() ) h->set_taxid( genome->taxid() ) ;
		}

		h->set_start_pos( minpos.is_reversed() ? start_pos-len+1 : start_pos ) ;
		h->set_aln_length( minpos.is_reversed() ? -len : len ) ;
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
		if( second_best ) r.set_diff_to_next( second_best.penalty - penalty ) ;
	}
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

