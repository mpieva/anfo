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

namespace streams {

//! \brief selects a policy appropriate for a given read
//! This simply checks whether a read is within a given length bin and
//! matches a name pattern, if given.  All policies found this way are
//! merged.  Checking for too many Ns was considered, but isn't actually
//! worthwhile since Ns don't contribute to seeds anyway.
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

string effective_sequence( const Read& rd )
{
	return rd.sequence().substr(
			rd.trim_left(),
			rd.has_trim_right() ? rd.trim_right() : std::string::npos
			) ;
}

void trim_cigar_left( google::protobuf::RepeatedField<uint32_t> &cig, unsigned len )
{
	for( int p = 0 ; p != cig.size() ; ++p )
	{
		switch( cigar_op( cig.Get( p ) ) )
		{
			case Hit::Match:
			case Hit::Mismatch:
			case Hit::Insert:
			case Hit::SoftClip:
				if( cigar_len( cig.Get( p ) ) > len )
				{
					cig.Set( p, mk_cigar(
								cigar_op( cig.Get( p ) ),
								cigar_len( cig.Get( p ) ) - len ) ) ;
					// remove beginning and get out
					if( !p ) return ;
					int q = 0 ;
					while( p != cig.size() ) cig.Set( q++, cig.Get( p++ ) ) ;
					while( q != cig.size() ) cig.RemoveLast() ;
					return ;
				}
				else len -= cigar_len( cig.Get( p ) ) ;
				break ;

			case Hit::Delete:
			case Hit::Skip:
			case Hit::Pad:
			case Hit::HardClip:
				break ;
		}
	}
	// nothing left if we get here
	cig.Clear() ;
}

void trim_cigar_right( google::protobuf::RepeatedField<uint32_t> &cig, unsigned len )
{
	for( int p = cig.size() ; p ; --p )
	{
		switch( cigar_op( cig.Get( p-1 ) ) )
		{
			case Hit::Match:
			case Hit::Mismatch:
			case Hit::Insert:
			case Hit::SoftClip:
				if( cigar_len( cig.Get( p-1 ) ) > len )
				{
					cig.Set( p-1, mk_cigar(
								cigar_op( cig.Get( p-1 ) ),
								cigar_len( cig.Get( p-1 ) ) - len ) ) ;
					// remove beginning and get out
					while( p != cig.size() ) cig.RemoveLast() ;
					return ;
				}
				else len -= cigar_len( cig.Get( p-1 ) ) ;
				break ;

			case Hit::Delete:
			case Hit::Skip:
			case Hit::Pad:
			case Hit::HardClip:
				break ;
		}
	}
	// nothing left if we get here
	cig.Clear() ;
}

void Housekeeper::put_result( const Result& rs )
{
	Stream::put_result( rs ) ;

	// trim adapters, set trim points
	// How does this work?  We create an overlap alignment, then
	// calculate an alignment score from the number of differences.  The
	// formula equals a match score of 1 and a mismatch score of -3.
	// Gaps score (mat-mis)/2 to allow for the simple algorithm.  If the
	// score is good enough, we trim.

	const output::Read &rd = res_.read() ;
	const string &seq = rd.sequence() ;

	unsigned r0 = rd.has_trim_right() ? rd.trim_right() : seq.length() ; 
	unsigned l0 = rd.trim_left() ;
	unsigned r = r0, l = l0 ;

	// first step: always trim trailing Ns (often caused by 454
	// failed quality trimming)
	while( r && ( seq[r-1] == 'n' || seq[r-1] == 'N') ) --r ;

	for( size_t i = 0 ; i != trim_right_.size() ; ++i )
	{
		int ymax, score = overlap_align(
				seq.rbegin() + ( seq.length() - r ), seq.rend(),
				trim_right_[i].rbegin(), trim_right_[i].rend(), &ymax ) ;
		if( score >= minscore_ && ymax > 0 ) r -= ymax ;
	}

	for( size_t i = 0 ; i != trim_left_.size() ; ++i )
	{
		int ymax, score = overlap_align(
				seq.begin() + l, seq.end(),
				trim_left_[i].begin(), trim_left_[i].end(), &ymax ) ;
		if( score >= minscore_ && ymax > 0 ) l += ymax ;
	}

	if( l > l0 )
	{
		for( int i = 0 ; i != res_.hit_size() ; ++i )
			trim_cigar_left( *res_.mutable_hit(i)->mutable_cigar(), l-l0 ) ;
		res_.mutable_read()->set_trim_right( l ) ;
	}
	if( r < r0 ) 
	{
		for( int i = 0 ; i != res_.hit_size() ; ++i )
			trim_cigar_right( *res_.mutable_hit(i)->mutable_cigar(), r-r0 ) ;
		res_.mutable_read()->set_trim_right( r ) ;
	}
}

Indexer::Indexer( const config::Config &config, const string& index_name ) :
	conf_( config ), index_name_( index_name ), index_( index_name )
{
	if( !conf_.policy_size() ) throw "no policies---nothing to do." ;
}

void Indexer::put_header( const Header& h )
{
	Stream::put_header( h ) ;
	hdr_.mutable_config()->mutable_policy()->MergeFrom( conf_.policy() ) ;
}

void Indexer::put_result( const Result& r )
{
	Stream::put_result( r ) ;
	Policy p = select_policy( conf_, res_.read()  ) ;
	AlnStats *as = res_.add_aln_stats() ;
	as->set_tag( index_.metadata().genome_name() ) ;

	int num_raw = 0, num_clumps = 0, num_useless = 0 ;
	for( int i = 0 ; i != p.use_compact_index_size() ; ++i )
	{
		const CompactIndexSpec &cis = p.use_compact_index(i) ;
		if( cis.name() == index_name_ || cis.name() + ".idx" == index_name_ ) {
            FixedIndex::LookupParams params ;
			params.cutoff = cis.has_cutoff() ? cis.cutoff() : numeric_limits<uint32_t>::max() ;
			params.allow_mismatches = cis.allow_near_perfect() ;
			params.wordsize = cis.has_wordsize() ? cis.wordsize() : index_.metadata().wordsize() ;
			params.stride = cis.has_stride() ? cis.stride() : index_.metadata().stride() ;

			assert( params.wordsize <= index_.metadata().wordsize() ) ;
			assert( index_.metadata().stride() % params.stride == 0 ) ;

			Seeds *ss = 0 ;
			for( int i = 0 ; !ss && i != res_.seeds_size() ; ++i )
				if( res_.seeds(i).genome_name() == index_.metadata().genome_name() )
					ss = res_.mutable_seeds(i) ;
			if( !ss ) 
            {
                ss = res_.add_seeds() ;
                ss->set_genome_name( index_.metadata().genome_name() ) ;
            }

			PreSeeds seeds ;
			num_raw += index_.lookupS(
					effective_sequence( res_.read() ), seeds, params, &num_useless ) ;
					
			num_clumps += cis.allow_near_perfect() 
				? combine_seeds( seeds, p.min_seed_len(), ss )
				: select_seeds( seeds, p.max_diag_skew(), p.max_gap(), p.min_seed_len(), index_.gaps(), ss ) ;

            if( p.has_repeat_threshold() && (unsigned)ss->ref_positions_size() >= p.repeat_threshold() )
            {
                ss->clear_ref_positions() ;
                ss->clear_query_positions() ;
                as->set_reason( output::too_many_seeds ) ;
            }
		}
	}
	as->set_num_raw_seeds( num_raw ) ;
	as->set_num_useless( num_useless ) ;
	as->set_num_clumps( num_clumps ) ;
}

Mapper::Mapper( const config::Aligner &config, const string& genome_name ) :
	conf_(config), genome_( Metagenome::find_genome( genome_name ) )
{
	simple_adna::pb = adna_parblock( conf_ ) ;
}

void Mapper::put_header( const Header& h )
{
	Stream::put_header( h ) ;
	hdr_.mutable_config()->mutable_aligner()->MergeFrom( conf_ ) ;
}

static inline bool icompare( const string& a, const string& b )
{
    if( a.size() != b.size() ) return false ;
    for( size_t i = 0 ; i != a.size() ; ++i )
        if( tolower( a[i] ) != tolower( b[i] ) ) return false ;
    return true ;
}

void Mapper::put_result( const Result& r )
{
	Stream::put_result( r ) ;

	AlnStats *as = 0 ;
    for( int i = 0 ; i != res_.aln_stats_size() ; ++i )
	{
		if( icompare( res_.aln_stats(i).tag(), genome_->name() )
				|| icompare( res_.aln_stats(i).tag(), genome_->name() + ".dna"  ) )
		{
			as = res_.mutable_aln_stats(i) ;
		}
	}
	
	if( !as ) {
		as = res_.add_aln_stats() ;
		as->set_tag( genome_->name() ) ;
	}

	Seeds ss ;
    for( int i = 0 ; i != res_.seeds_size() ; ++i )
	{
		if( icompare( res_.seeds(i).genome_name(), genome_->name() )
				|| icompare( res_.seeds(i).genome_name(), genome_->name() + ".dna"  ) )
		{
			ss.Swap( res_.mutable_seeds(i) ) ;
			res_.mutable_seeds()->SwapElements( i, res_.seeds_size()-1 ) ;
			res_.mutable_seeds()->RemoveLast() ;
			break ;
		}
	}

    // is it already clear that we cannot do anything?
    if( as->has_reason() ) return ;

	// not seeded means no policy (or logic bug, but let's not go there...)
	if( ss.genome_name().empty() ) {
		as->set_reason( output::no_policy ) ;
		return ;
	}

	// invariant violated, we won't deal with that
	if( ss.ref_positions_size() != ss.query_positions_size() )
		throw "invalid seeds: coordinates must come in pairs" ;

	if( !ss.ref_positions_size() ) {
		as->set_reason( as->num_useless() ? output::repeats_only : output::no_seeds ) ;
		return ;
	}

	QSequence qs( res_.read() ) ;
	uint32_t o, c, tt, max_penalty = (uint32_t)( conf_.max_penalty_per_nuc() * qs.length() ) ;

	std::deque< alignment_type > ol ;
	for( int i = 0 ; i != ss.ref_positions_size() ; ++i )
		ol.push_back( alignment_type( *genome_, qs, ss.ref_positions(i), ss.query_positions(i) ) ) ;

	alignment_type::ClosedSet cl ;
	alignment_type best = find_cheapest( ol, cl, max_penalty, &o, &c, &tt ) ;
	res_.mutable_seeds()->RemoveLast() ;

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

	output::Hit *h = res_.add_hit() ;

	uint32_t start_pos ;
    const config::Sequence *sequ = genome_->translate_back( minpos+1, start_pos ) ;
	if( !sequ ) throw "Not supposed to happen:  invalid alignment coordinates" ;

	if( genome_->g_.has_name() ) h->set_genome_name( genome_->g_.name() ) ;
	h->set_sequence( sequ->name() ) ;
	if( sequ->has_taxid() ) h->set_taxid( sequ->taxid() ) ;
	else if( genome_->g_.has_taxid() ) h->set_taxid( genome_->g_.taxid() ) ;

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

	// get rid of overlaps of that first alignment, then look
	// for the next one
	// XXX this is cumbersome... need a better PQueue impl... or a
	// better algorithm
	ol.erase( 
			std::remove_if( ol.begin(), ol.end(), reference_overlaps( minpos, maxpos ) ),
			ol.end() ) ;
	make_heap( ol.begin(), ol.end() ) ;

	// search long enough to make sensible mapping quality possible
	uint32_t max_penalty_2 = conf_.max_mapq() + penalty ;
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

} ; // namespace 
