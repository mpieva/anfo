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

string expand( const string& s, int x )
{
    if( s.size() <= 1 ) return s ;

    char u = 'a' + (x%26), v = 'a' + (x/26) ;
    string r ;
    size_t i = 0 ;
    for( ; i+1 < s.size() ; ++i )
    {
        if( s[i] == '%' && s[i+1] == '%' ) {
            r.push_back( v ) ;
            r.push_back( u ) ;
            ++i ;
        }
        else r.push_back( s[i] ) ;
    }
    if( i < s.size() ) r.push_back( s[i] ) ;
    return r ;
}

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

Mapper::Mapper( const config::Config &config, const string& index_name ) :
	mi_(config), index_name_(index_name) 
{
	if( !mi_.has_aligner() ) throw "no aligner configuration---cannot start." ;
	if( !mi_.policy_size() ) throw "no policies---nothing to do." ;

	simple_adna::pb = adna_parblock( mi_.aligner() ) ;

	FixedIndex( index_name ).swap( index_ ) ;
	genome_ = Metagenome::find_genome( index_.ci_.genome_name() ) ;
}

void Mapper::index_sequence( output::Result &r ) // XXX , QSequence &qs ) // , std::deque< alignment_type >& ol )
{
	static const int minscore = 4 ;

	// trim adapters, set trim points
	// How does this work?  We create an overlap alignment, then
	// calculate an alignment score from the number of differences.  The
	// formula equals a match score of 1 and a mismatch score of -3.
	// Gaps score (mat-mis)/2 to allow for the simple algorithm.  If the
	// score is good enough, we trim.

	const output::Read &rd = r.read() ;
	const string &seq = rd.sequence() ;

	if( mi_.trim_right_size() || mi_.trim_left_size() ) 
	{
		unsigned n = rd.has_trim_right() ? rd.trim_right() : seq.length() ; 
		while( n && ( seq[n-1] == 'n' || seq[n-1] == 'N') ) --n ;
		r.mutable_read()->set_trim_right( n ) ;
	}

	for( int i = 0 ; i != mi_.trim_right_size() ; ++i )
	{
		int ymax, score = overlap_align(
				seq.rbegin() + ( rd.has_trim_right() ? seq.length() - rd.trim_right() : 0 ), seq.rend(),
				mi_.trim_right(i).rbegin(), mi_.trim_right(i).rend(), &ymax ) ;
		if( score >= minscore && ymax > 0 )
			r.mutable_read()->set_trim_right(
					(rd.has_trim_right() ? rd.trim_right() : seq.length()) - ymax ) ;
	}

	for( int i = 0 ; i != mi_.trim_left_size() ; ++i )
	{
		int ymax, score = overlap_align(
				seq.begin() + rd.trim_left(), seq.end(),
				mi_.trim_left(i).begin(), mi_.trim_left(i).end(), &ymax ) ;
		if( score >= minscore && ymax > 0 )
			r.mutable_read()->set_trim_left( rd.trim_left() + ymax ) ;
	}

	Policy p = select_policy( mi_, r.read() ) ;
	// XXX QSequence( rd ).swap( qs ) ;

	int num_raw = 0, num_comb = 0, num_clumps = 0, num_useless = 0 ;
	for( int i = 0 ; i != p.use_compact_index_size() ; ++i )
	{
		const CompactIndexSpec &cis = p.use_compact_index(i) ;
		if( cis.name() == index_name_ ) {
			FixedIndex::LookupParams params ;
			params.cutoff = cis.has_cutoff() ? cis.cutoff() : numeric_limits<uint32_t>::max() ;
			params.allow_mismatches = cis.allow_near_perfect() ;
			params.wordsize = index_.ci_.wordsize() ;
			params.stride = index_.ci_.stride() ;

			// const FixedIndex &ix = indices[ cis.name() ] ;
			// console.output( Console::warning, "opening genome " + ix.ci_.genome_name() ) ;
			// GenomeHolder g = Metagenome::find_genome( ix.ci_.genome_name() ) ;
			// XXX console.output( Console::warning, "done" ) ;
			// if( !ix || !g->get_base() )
			// throw "couldn't find requested index " + cis.name()
			// + " or genome " + ix.ci_.genome_name() ;

			Seeds *ss = 0 ;
			for( int i = 0 ; !ss && i != r.seeds_size() ; ++i )
				if( r.seeds(i).genome_name() == index_.ci_.genome_name() )
					ss = r.mutable_seeds(i) ;
			if( !ss ) ss = r.add_seeds() ;

			ss->set_max_penalty_per_nuc( max( p.max_penalty_per_nuc(), ss->max_penalty_per_nuc() ) ) ;
			if( p.has_repeat_threshold() ) ss->set_repeat_threshold( p.repeat_threshold() ) ;

			// XXX soooo ugly (and slow)
			PreSeeds seeds ;
			num_raw += index_.lookupS( 
					seq.substr( rd.trim_left(), rd.has_trim_right() ? rd.trim_right() : std::string::npos ),
					seeds, params, &num_useless ) ;
					
			// num_comb += deep_count( seeds ) ;
			num_clumps += cis.allow_near_perfect() 
				? combine_seeds( seeds, p.min_seed_len(), ss )
				: select_seeds( seeds, p.max_diag_skew(), p.max_gap(), p.min_seed_len(), genome_->get_contig_map(), ss ) ;

			// XXX num_clumps += seeds.size() ;

			// g->add_ref() ; // XXX dirty hack
			// setup_alignments( *g, qs, seeds.begin(), seeds.end(), ol ) ;
		}
	}
	AlnStats *as = r.add_aln_stats() ;
	as->set_num_raw_seeds( num_raw ) ;
	as->set_num_useless( num_useless ) ;
	// as->set_num_grown_seeds( num_comb ) ;
	as->set_num_clumps( num_clumps ) ;
	if( p.has_tag() ) as->set_tag( p.tag() ) ;
}

void Mapper::process_sequence( /*const QSequence &ps, double max_penalty_per_nuc, std::deque< alignment_type > &ol,*/ output::Result &r )
{
	AlnStats *as = r.aln_stats_size() ? r.mutable_aln_stats( r.aln_stats_size()-1 ) : r.add_aln_stats() ;

	// not seeded means no policy (or logic bug, but let's not go there...)
	if( !r.seeds_size() ) {
		as->set_reason( output::no_policy ) ;
		return ;
	}

	const Seeds &ss = r.seeds( r.seeds_size()-1 ) ;

	// invariant violated, we won't deal with that
	if( ss.ref_positions_size() != ss.query_positions_size() )
		throw "invalid seeds: coordinates must come in pairs" ;

	if( !ss.ref_positions_size() ) {
		as->set_reason( as->num_useless() ? output::repeats_only : output::no_seeds ) ;
		return ;
	}

	if( ss.has_repeat_threshold() && (unsigned)ss.ref_positions_size() >= ss.repeat_threshold() ) {
		as->set_reason( output::too_many_seeds ) ;
		return ;
	}

	GenomeHolder g = Metagenome::find_genome( ss.genome_name() ) ;

	QSequence qs( r.read() ) ;
	uint32_t o, c, tt, max_penalty = (uint32_t)( ss.max_penalty_per_nuc() * qs.length() ) ;

	std::deque< alignment_type > ol ;
	for( int i = 0 ; i != ss.ref_positions_size() ; ++i )
		ol.push_back( alignment_type( *g, qs, ss.ref_positions(i), ss.query_positions(i) ) ) ;

	alignment_type::ClosedSet cl ;
	alignment_type best = find_cheapest( ol, cl, max_penalty, &o, &c, &tt ) ;
	r.mutable_seeds()->RemoveLast() ;

	as->set_open_nodes_after_alignment( o ) ;
	as->set_closed_nodes_after_alignment( c ) ;
	as->set_tracked_closed_nodes_after_alignment( tt ) ;

	if( !best ) {
		as->set_reason( output::bad_alignment ) ;
		return ;
	}

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

