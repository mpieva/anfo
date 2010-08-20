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
            rd.has_trim_right() ? rd.trim_right() - rd.trim_left(): std::string::npos
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
	conf_( config ), index_name_( index_name ), index_( MetaIndex::add_ref( index_name ) )
{
	if( !conf_.policy_size() ) throw "no policies---nothing to do." ;
}

Indexer::~Indexer()
{
	MetaIndex::free_ref( index_ ) ;
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
				ss->set_max_mismatches( p.max_mismatches_in_seed() ) ;
				ss->set_min_seed_length( p.min_seed_len() ) ;
            }

			PreSeeds seeds ;
			num_raw += index_.lookupS(
					effective_sequence( res_.read() ), seeds, params, &num_useless ) ;
					
			num_clumps += combine_seeds( seeds, p.min_seed_len(), ss ) ;

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
	parblock_ = adna_parblock( conf_ ) ;
}

void Mapper::put_header( const Header& h )
{
	Stream::put_header( h ) ;
	hdr_.mutable_config()->mutable_aligner()->MergeFrom( conf_ ) ;
}

// Check for a stretch of at least minsize matches with up to mm
// mismatches.  The pointers should be clean: both point to the start of
// a decently matching region of length seedlen.
static bool check_seed_quality( DnaP ref, const QSequence::Base *qry, int seedlen, int mm, int minsize )
{
	int lengths[3] = { 0,0,0 } ;
	for( ; seedlen ; --seedlen, ++ref, ++qry )
	{
		if( *ref == qry->ambicode )
		{
			++lengths[0] ;
			++lengths[1] ;
			++lengths[2] ;
		}
		else
		{
			lengths[2] = lengths[1] + 1 ;
			lengths[1] = lengths[0] + 1 ;
			lengths[0] = 0 ;
		}
		if( lengths[mm] == minsize ) return true ;
	}
	return false ;
}

void Mapper::put_result( const Result& r )
{
	Stream::put_result( r ) ;

	AlnStats *as = 0 ;
    for( int i = 0 ; i != res_.aln_stats_size() ; ++i )
	{
		if( res_.aln_stats(i).tag() == genome_->name() )
			as = res_.mutable_aln_stats(i) ;
	}
	
	if( !as ) {
		as = res_.add_aln_stats() ;
		as->set_tag( genome_->name() ) ;
	}

	Seeds ss ;
    for( int i = 0 ; i != res_.seeds_size() ; ++i )
	{
		if( res_.seeds(i).genome_name() == genome_->name() )
		{
			ss.Swap( res_.mutable_seeds(i) ) ;
            if( i != res_.seeds_size()-1 ) res_.mutable_seeds()->SwapElements( i, res_.seeds_size()-1 ) ;
            if( res_.seeds_size() > 1 ) res_.mutable_seeds()->RemoveLast() ;
            else res_.clear_seeds() ;
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
	if( ss.ref_positions_size() != ss.query_positions_size() ||
			ss.ref_positions_size() != ss.seed_sizes_size() )
		throw "invalid seeds: coordinates must come in triples" ;

	if( !ss.ref_positions_size() ) {
		as->set_reason( as->num_useless() ? output::repeats_only : output::no_seeds ) ;
		return ;
	}

	QSequence qs( res_.read() ) ;
	Logdom max_penalty = Logdom::from_phred( conf_.max_penalty_per_nuc() * qs.length() ) ;

	// do actual alignments:
	// 1 initialize each one, make sure the mismatch count is low enough
	//   [need to pass mismatch limit through from Indexer]
	// 2 iterate, increasing the limit if necessary:
	// 2a evaluate alignments, keep the two best scores and associated seeds
	// 2b remove seeds that already produced an alignment
	// 3 redo winning alignment and backtrace it

	std::deque< SeededAlignment > seedlist ;
	for( int i = 0 ; i != ss.ref_positions_size() ; ++i )
	{
		// Reconstruct pointers to seed region.  We want the region to
		// be oriented forwards on the query.  That means:
		// - for a forward seed, the ref position is correct, the query
		//   is one too big because of the virtual gap at 0
		// - for a reverse seed we have coordinates for the *end*
		//   position of the seed region, and we need to reverse the
		//   pointer on the reference
		DnaP reference = genome_->get_base() + ss.ref_positions(i) ;
		int32_t qoffs = ss.query_positions(i) ;
		uint32_t size = ss.seed_sizes(i) ;
			
		if( check_seed_quality( reference, qs.start() + qoffs, size,
					ss.max_mismatches(), ss.min_seed_length() ) )
			seedlist.push_back( SeededAlignment( parblock_, reference, qs, qoffs, size ) ) ;
	}

	// iteration: we track the best score along with its seed and the
	// second best score
	// XXX: if we find the first alignment close to the limit, we must
	// increase(!) the limit to best_score*maxq and do another iteration
	// (actually the whole limit business is shaky right now)
	
	SeededAlignment best_seed ;
	ExtendBothEnds best_ext ;

	Logdom best_score = Logdom::null(),
		   runnerup_score = Logdom::null(),
	       limit = Logdom::from_phred( 60 ) ;

	while( !seedlist.empty() )
	{
		// std::cerr << "Starting " << seedlist.size() << " alignments pass at limit " << limit.to_phred() << std::endl ;
		std::deque< SeededAlignment >::iterator
			cur_aln( seedlist.begin() ), end_aln( seedlist.end() ), out_aln( seedlist.begin() ) ;
		while( cur_aln != end_aln )
		{
			ExtendBothEnds extension( parblock_, qs, *cur_aln, limit ) ;
			Logdom score = extension.score_ ;

			// found an alignment?  finish with this seed
			if( score.is_finite() )
			{
				// new best score?
				if( score > best_score ) {
					runnerup_score = best_score ;
					best_score = score ;
					best_seed = *cur_aln ;
					best_ext.swap( extension ) ;
					limit = max( limit, max( best_score * Logdom::from_phred( conf_.max_mapq() ), runnerup_score ) ) ;
				}
				// new second best score?
				else if( score > runnerup_score ) {
					runnerup_score = score ;
					limit = max( limit, runnerup_score ) ;
				}
				// this seed is exhausted, we never need it again
				++cur_aln ;
				// std::cerr << "Got " << score.to_phred() << ", new limit is " << limit.to_phred() 
					// << ", top two are " << best_score.to_phred() << " and " << runnerup_score.to_phred() << std::endl ;
			}
			else // no alignment: move to next seed
			{
				// std::cerr << "Nothing yet" << std::endl ;
				if( out_aln != cur_aln ) *out_aln = *cur_aln ;
				out_aln++, cur_aln++ ;
			}
		}
		// two hits -> we're done (regardless of score, we got
		// everything)
		if( runnerup_score.is_finite() ) break ;

		// what's the current absolute limit?  if we got a hit, it's
		// worse by max mapq, else it's the max penalty
		Logdom abs_max = best_score.is_finite()
			? best_score * Logdom::from_phred( conf_.max_mapq() ) : max_penalty ;

		// already over? -> we're done
		if( limit <= abs_max ) break ;

		// get rid of exhausted seeds, increase limit
		seedlist.erase( out_aln, end_aln ) ;

		// new limit: twice as high, but not unnecessarily high
		limit = max( limit*limit, abs_max ) ;
	}

	if( !best_score.is_finite() )
	{
		as->set_reason( output::bad_alignment ) ;
		return ;
	}

	DnaP minpos, maxpos ;
	std::vector<unsigned> t = best_ext.backtrace( best_seed, minpos, maxpos ) ;
	int32_t len = best_seed.qoffs_ > 0 ? maxpos - minpos : minpos - maxpos ;
	assert( maxpos > minpos && !maxpos.is_reversed() && !minpos.is_reversed() ) ;

	output::Hit *h = res_.add_hit() ;

	uint32_t start_pos ;
	const config::Sequence *sequ = genome_->translate_back( minpos, start_pos ) ;
	if( !sequ ) throw "Not supposed to happen:  invalid alignment coordinates" ;

	if( genome_->g_.has_name() ) h->set_genome_name( genome_->g_.name() ) ;
	h->set_sequence( sequ->name() ) ;
	if( sequ->has_taxid() ) h->set_taxid( sequ->taxid() ) ;
	else if( genome_->g_.has_taxid() ) h->set_taxid( genome_->g_.taxid() ) ;

	h->set_start_pos( start_pos ) ;
	h->set_aln_length( len ) ;
	h->set_score( best_score.to_phred() ) ;
	std::copy( t.begin(), t.end(), RepeatedFieldBackInserter( h->mutable_cigar() ) ) ;

	//! \todo calculate a real E-value
	// XXX h->set_evalue

	//! \todo Find second best hit and similar stuff.
	//! We want the distance to the next best hit; also,
	//! unless already found, we want the best hit to some
	//! selected genome(s).

	// XXX set diff_to_next_species, diff_to_next_order

	if( runnerup_score.is_finite() ) h->set_diff_to_next( (runnerup_score/best_score).to_phred() ) ;
}

//! \page finding_alns How to find everything we need
//! We look for best hits globally and specifically on one genome.  They
//! will be discovered in the order of decreasing score.  After setup,
//! we can operate in a loop of cleaning out the stuff we don't need
//! anymore and finding more alignments.
//!
//! Cleanup: XXX this is outdated
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
