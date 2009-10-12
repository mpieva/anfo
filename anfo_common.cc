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

#include "config.h"
#include "anfo_common.h"

#include "stream.h"
#include "trim.h"

#include "config.pb.h"
#include "output.pb.h"

//! \todo We want an E-value...
//! \todo We want more than just the best match.  Think about a sensible
//!       way to configure this.
//! \todo Test this: the canonical test case is homo sapiens, chr 21.

using namespace config ;
using namespace output ;
using namespace std; 

Policy select_policy( const Config &c, const Read &r )
{
	Policy p ;
	for( int i = 0 ; i != c.policy_size() ; ++i )
	{
		const Policy &pi = c.policy(i) ;
		if( ( !pi.has_minlength() || pi.minlength() <= r.sequence().length() ) &&
			( !pi.has_maxlength() || pi.maxlength() >= r.sequence().length() ) &&
			( !pi.has_name_pattern() || 0 == fnmatch( pi.name_pattern().c_str(), r.seqid().c_str(), 0 ) ) )
			p.MergeFrom( pi ) ;
	}
	return p ;
}

Mapper::Mapper( const config::Config &config ) : mi(config)
{
	if( !mi.has_aligner() ) throw "no aligner configuration---cannot start." ;
	if( !mi.policy_size() ) throw "no policies---nothing to do." ;

    for( int i = mi.genome_path_size() ; i != 0 ; --i )
        Metagenome::add_path( mi.genome_path(i-1) ) ;

	simple_adna::configure( mi.aligner() ) ;
	for( int i = 0 ; i != mi.policy_size() ; ++i )
	{
		for( int j = 0 ; j != mi.policy(i).use_compact_index_size() ; ++j )
		{
			CompactIndexSpec &ixs = *mi.mutable_policy(i)->mutable_use_compact_index(j) ;
			FixedIndex &ix = indices[ ixs.name() ] ;
			if( !ix ) {
				FixedIndex( ixs.name(), mi, MADV_WILLNEED ).swap( ix ) ;
				const string& genome_name = ix.ci_.genome_name() ; 
				Metagenome::find_genome( genome_name, Metagenome::persistent ) ;
			}
		}
	}
}

static const int maxd = 32 ;
static const int minscore = 4 ;

int Mapper::index_sequence( output::Result &r, QSequence &qs, std::deque< alignment_type >& ol )
{
	// trim adapters, set trim points
	// How does this work?  We create an overlap alignment, then
	// calculate an alignment score from the number of differences.  The
	// formula equals a match score of 1 and a mismatch score of -3.
	// Gaps score (mat-mis)/2 to allow for the simple algorithm.  If the
	// score is good enough, we trim.

	if( mi.trim_right_size() || mi.trim_left_size() ) 
	{
		unsigned n = r.read().sequence().length() ; 
		while( n && ( r.read().sequence()[n-1] == 'n' || r.read().sequence()[n-1] == 'N') ) --n ;
		r.mutable_read()->set_trim_right( n ) ;
	}

	for( int i = 0 ; i != mi.trim_right_size() ; ++i )
	{
		int ymax, xmax = mi.trim_right(i).size() ;
		// XXX: already got a trim point
		int diff = align(
				r.read().sequence().rbegin(), r.read().sequence().rend(),
				mi.trim_right(i).rbegin(), mi.trim_right(i).rend(),
				maxd, overlap, 0, &ymax ) ;
		int score = xmax + ymax - 8 * diff ;
		if( diff < maxd && score >= minscore && ymax > 0 )
			r.mutable_read()->set_trim_right( r.read().sequence().length() - ymax ) ;
	}

	for( int i = 0 ; i != mi.trim_left_size() ; ++i )
	{
		int ymax, xmax = mi.trim_left(i).size() ;
		// XXX already got trim points
		int diff = align(
				r.read().sequence().begin(), r.read().sequence().end(),
				mi.trim_left(i).begin(), mi.trim_left(i).end(),
				maxd, overlap, 0, &ymax ) ;
		int score = xmax + ymax - 8 * diff ;
		if( diff < maxd && score >= minscore && ymax > 0 )
			r.mutable_read()->set_trim_left( r.read().trim_left() + ymax ) ;
	}

	Policy p = select_policy( mi, r.read() ) ;
	QSequence( r.read() ).swap( qs ) ;

	int num_raw = 0, num_comb = 0, num_clumps = 0, num_useless = 0 ;
	for( int i = 0 ; i != p.use_compact_index_size() ; ++i )
	{
		const CompactIndexSpec &cis = p.use_compact_index(i) ;
		const FixedIndex &ix = indices[ cis.name() ] ;
		const CompactGenome &g = Metagenome::find_genome( ix.ci_.genome_name(), Metagenome::persistent ) ;
		assert( ix ) ; assert( g.get_base() ) ;

		vector<Seed> seeds ;
		num_raw += ix.lookupS( 
				r.read().sequence(), seeds, cis.allow_near_perfect(), &num_useless,
				cis.has_cutoff() ? cis.cutoff() : numeric_limits<uint32_t>::max() ) ;
		num_comb += seeds.size() ;
		select_seeds( seeds, p.max_diag_skew(), p.max_gap(), p.min_seed_len(), g.get_contig_map() ) ;
		num_clumps += seeds.size() ;

		setup_alignments( g, qs, seeds.begin(), seeds.end(), ol ) ;
	}
	AlnStats *as = r.mutable_aln_stats() ;
	as->set_num_raw_seeds( num_raw ) ;
	as->set_num_useless( num_useless ) ;
	as->set_num_grown_seeds( num_comb ) ;
	as->set_num_clumps( num_clumps ) ;

	if( !p.has_max_penalty_per_nuc() )
	{
		as->set_reason( output::no_policy ) ;
		return INT_MAX ;
	}
	else if( ol.empty() ) 
	{
		as->set_reason( num_useless ? output::repeats_only : output::no_seeds ) ;
		return INT_MAX ;
	}
	else if( p.has_repeat_threshold() && ol.size() >= p.repeat_threshold() )
	{
		as->set_reason( output::too_many_seeds ) ;
		return INT_MAX ;
	}
	else return p.max_penalty_per_nuc() ;
}

void Mapper::process_sequence( const QSequence &ps, double max_penalty_per_nuc, std::deque< alignment_type > &ol, output::Result &r )
{
	uint32_t o, c, tt, max_penalty = (uint32_t)( max_penalty_per_nuc * ps.length() ) ;
	alignment_type::ClosedSet cl ;
	alignment_type best = find_cheapest( ol, cl, max_penalty, &o, &c, &tt ) ;

	AlnStats *as = r.mutable_aln_stats() ;
	as->set_open_nodes_after_alignment( o ) ;
	as->set_closed_nodes_after_alignment( c ) ;
	as->set_tracked_closed_nodes_after_alignment( tt ) ;
	if( !best )
	{
		as->set_reason( output::bad_alignment ) ;
	}
	else
	{
		int penalty = best.penalty ;

		deque< pair< alignment_type, const alignment_type* > > ol_ ;
		reset( best ) ;
		greedy( best ) ;
		(enter_bt<alignment_type>( ol_ ))( best ) ;
		DnaP minpos, maxpos ;
		std::vector<unsigned> t = find_cheapest( ol_, minpos, maxpos ) ;
		int32_t len = maxpos - minpos - 1 ;

		output::Hit *h = r.add_hit() ;

		uint32_t start_pos ;
		const Sequence *sequ ;
		const Genome *genome ;

		if( !Metagenome::translate_to_genome_coords( minpos+1, start_pos, &sequ, &genome ) )
			throw "Not supposed to happen:  invalid alignment coordinates" ;
		
		if( genome->has_name() ) h->set_genome_name( genome->name() ) ;
		h->set_sequence( sequ->name() ) ;
		if( sequ->has_taxid() ) h->set_taxid( sequ->taxid() ) ;
		else if( genome->has_taxid() ) h->set_taxid( genome->taxid() ) ;

		h->set_start_pos( minpos.is_reversed() ? start_pos-len+1 : start_pos ) ;
		h->set_aln_length( minpos.is_reversed() ? -len : len ) ;
		h->set_score( penalty ) ;
		std::copy( t.begin(), t.end(), RepeatedFieldBackInserter( h->mutable_cigar() ) ) ;

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

		// search long enough to make sensible mapping quality possible
		uint32_t max_penalty_2 = 254 + penalty ;
		alignment_type second_best = find_cheapest( ol, cl, max_penalty_2, &o, &c, &tt ) ;
		as->set_open_nodes_after_alignment( o ) ;
		as->set_closed_nodes_after_alignment( c ) ;
		as->set_tracked_closed_nodes_after_alignment( tt ) ;
		if( second_best ) h->set_diff_to_next( second_best.penalty - penalty ) ;
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

