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

#define NDEBUG

#include "align.h"

#include "config.pb.h"

#include <ostream>
#include <iomanip>

namespace {
	void print_subst_mat( std::ostream& s, const subst_mat mat )
	{
		s << "   r\\q " ;
		for( Ambicode qry = 1 ; qry != 16 ; ++qry )
		{
			s << std::setw(2) << from_ambicode( qry ) << "  " ;
		}
		for( Ambicode ref = 1 ; ref != 16 ; ++ref )
		{
			s << "\n   " << from_ambicode( ref ) << "   " ;
			for( Ambicode qry = 1 ; qry != 16 ; ++qry )
			{
				s << std::setw(2) << mat[ref][qry].to_phred() << "  " ;
			}
		}
		s << '\n' ;
	}
}

//! \brief sets up aDNA alignments according to parameter set
//! Penalties are logs of probabilities.  We calculate a probability in
//! a straightforward manner, then take the log, then multiply by (say)
//! 10 and round to an integer (effectively shifting the base).
//! Afterwards we could shift the score up so the minimum score is zero,
//! but it will be zero to within accuracy anyway.  Here we
//! precalculate the penalties and log-transform them, however, the
//! substitution matrices are left in the linear scale, because they
//! need to be added later.
//!
//! For ambiguity codes, we sum over all the possible codes, divide by
//! the number of possibilities, and later (when aligning) go to the log
//! domain.  Note that matches end up costing nothing (almost nothing,
//! but accuracy isn't good enough to pick up on that), matching to an N
//! *does* cost something, and matching a T to a C is dirt cheap while
//! we believe we are in single stranded DNA.
//!
//! For gaps, we already have probabilities.  Taking the log is all it
//! takes to make scores from that.
//!
//! For the overhang, we want something that's easy to parameterize, and
//! the single parameter will be the average overhang length.  We assume
//! there's an even chance that we have a 5' overhang (as opposed to a
//! 3' overhang), and if we have a 5' overhang, it's length follows a
//! geometric distribution with parameter p.
//! 
//! The mean overhang length ends up being (1-p)/(2p), which is just
//! half the mean of the geometric distribution.  Solving for p derives
//! the parameter from the configuration, but the code operates on \c
//! p_comp, which is just 1-p.  Entering an overhang involves a penalty
//! (since the probability is only one half), but not entering it
//! involves the same penalty.  The difference in penalties is zero,
//! which also means a perfect alignment stays at a zero score.
//!
//! The penalty for extending an overhang by one nucleotide is just the
//! probability of it being longer.
//!
//! \param conf set of configuration parameters
//! \param out if not NULL, receives a protocol of the calculated
//!            internal parameters
adna_parblock::adna_parblock( const config::Aligner& conf )
{
	gap_ext_penalty = Logdom::from_float( conf.gap_extension_rate() ) ;
	gap_open_penalty = conf.has_gap_open_rate() ? Logdom::from_float( conf.gap_open_rate() ) : gap_ext_penalty ;

	if( conf.has_mean_overhang_length() )
	{
		double p_comp = 1 - 1 / (2*conf.mean_overhang_length() + 1) ;
		overhang_enter_penalty = Logdom::one() ;
		overhang_ext_penalty = Logdom::from_float( p_comp ) ;
	}
	else
	{
		overhang_enter_penalty = Logdom::null() ;
		overhang_ext_penalty = Logdom::one() ;
	}

	double tv = conf.rate_of_transversions() ;
	double ts = conf.has_rate_of_transitions() ? conf.rate_of_transitions() : tv ;
	double dds = conf.rate_of_ds_deamination() ;
	double dss = conf.rate_of_ss_deamination() ;
	double r = 1-ts-tv-tv ;

	double prim_mat_ds[4][4] = { // to A C T G
		{      r,    tv,     tv,    ts }, // from A
		{     tv, r-dds, ts+dds,    tv }, // from C
		{     tv,    ts,      r,    tv }, // from T
		{ ts+dds,    tv,     tv, r-dds }  // from G
	} ;
	double prim_mat_ss[4][4] = { // to A C T G
		{      r,    tv,     tv,    ts }, // from A
		{     tv, r-dss, ts+dss,    tv }, // from C
		{     tv,    ts,      r,    tv }, // from T
		{ ts+dss,    tv,     tv, r-dss }  // from G
	} ;

	// filling the matrices
	// The entry for any ambiguity code is the sum of the entries in the
	// prim. matrices over all possible pairs of encoded nucleotides,
	// divided by the number of pairs, then converted to log domain.
	// (Yeah, it's inefficient.  It also runs only once per invocation,
	// so sue me.)
	for( Ambicode ref = 1 ; ref != 16 ; ++ref )
	{
		for( Ambicode qry = 1 ; qry != 16 ; ++qry )
		{
			double p_ds = 0, p_ss = 0 ;
			int npairs = 0 ;
			for( int ref_n = 0 ; ref_n != 4 ; ++ref_n )
			{
				for( int qry_n = 0 ; qry_n != 4 ; ++qry_n )
				{
					if( (ref & (1<<ref_n)) && (qry & (1<<qry_n)) )
					{
						++npairs ;
						p_ds += prim_mat_ds[ref_n][qry_n] ;
						p_ss += prim_mat_ss[ref_n][qry_n] ;
					}
				}
			}

			ds_mat[ref][qry] = Logdom::from_float( p_ds / npairs ) ;
			ss_mat[ref][qry] = Logdom::from_float( p_ss / npairs ) ;
		}
	}
}

std::ostream& operator << ( std::ostream& s, const adna_parblock& p )
{
	s << "\n\e[1mALIGNMENT PARAMETERS\e[0m:\n\n"
		<< "  \e[4mDouble Stranded Substitution Matrix\e[0m:\n\n" ;
	print_subst_mat( s, p.ds_mat ) ;
	s << "\n  \e[4mSingle Stranded Substitution Matrix\e[0m:" ;
	if( p.ss_mat[1][1] != p.ds_mat[1][1] ) {
		s << "\n\n" ; print_subst_mat( s, p.ss_mat ) ;
	} else s << " N/A\n" ;
	s << "\n  \e[4mGap Open Penalty\e[0m: " << p.gap_open_penalty.to_phred()
		<< "\n  \e[4mGap Extension Penalty\e[0m: " << p.gap_ext_penalty.to_phred()
		<< "\n\n  \e[4mOverhang Penalty\e[0m: " ;
	if( !p.overhang_enter_penalty.is_finite() ) s << "N/A" ; else s << p.overhang_enter_penalty.to_phred() ;
	s << "\n  \e[4mOverhang Extension Penalty\e[0m: " ;
	if( !p.overhang_ext_penalty.is_finite() ) s << "N/A" ; else s << p.overhang_ext_penalty.to_phred() ;
	return s << '\n' << std::endl ;
}


SeededAlignment::SeededAlignment(
		const adna_parblock& pb,
		DnaP reference,
		const QSequence& query,
		int qoffs,
		int size )
	: reference_( reference ), score_( Logdom::one() ), qoffs_( qoffs ), size_( size )
{
	// greedy initialization: run over the seed, accumulating a
	// score.  then extend greedily as long as there are
	// matches, store the resulting score.
	for( size_t i = 0 ; i <= size_ ; ++i )
		score_ *= pb.subst_penalty( 0, reference_[i], query[qoffs_+i] ) ;

	for( ; reference_[size_] && query[qoffs_+size_].ambicode &&
			reference_[size_] == query[qoffs_+size_].ambicode ; ++size_ )
		score_ *= pb.subst_penalty( 0, reference_[size_], query[qoffs_+size_] ) ;

	for( ; reference_[-1] && query[qoffs_-1].ambicode &&
			reference_[-1] == query[qoffs_-1].ambicode ;
			++size_, --reference_, --qoffs_ )
		score_ *= pb.subst_penalty( 0, reference_[-1], query[qoffs_-1] ) ;
}

// multiple calls to put must happen in increasing order of x
inline void ExtendAlignment::put( int s, int os, int x, int xo, int y, int yo, Logdom z )
{
	if( z >= limit_ )
	{
		if( mins_[ y+yo ] == maxs_[ y+yo ] ) mins_[y+yo] = maxs_[y+yo] = x ;
		for( ; maxs_[y+yo] != x ; ++maxs_[y+yo] )
			for( int s = 0 ; s != adna_parblock::num_states ; ++s )
				cells_[ width_*(y+yo) + maxs_[y+yo] ][ s ].score = Logdom::null() ;

		assert( y+yo < width_+1 ) ;
		assert( x+xo < width_ ) ;
		assert( width_*(y+yo) + x+xo < (1+width_) * width_ ) ;

		Cell& c = cells_[ width_*(y+yo) + x+xo ][ s ] ;
		if( c.score < z )
		{
			c.score = z ;
			c.from_state = os ;
			c.from_x_offset = xo ;
			c.from_y_offset = yo ;
		}
	}
}

static inline int query_length( const QSequence::Base *query )
{
	int r = 0 ;
	while( query->ambicode ) ++query, ++r ;
	return r ;
}

// Alignment proper: the intial greedy matching must have been done,
// here we extend one side of this into a full alignment, as long as it
// doesn't score more than a prescribed limit.
ExtendAlignment::ExtendAlignment( const adna_parblock& pb, DnaP reference, const QSequence::Base *query, Logdom limit ) :
	width_( query_length( query )+1 ), cells_( (1+width_)*width_ ), mins_( 1+width_ ), maxs_( 1+width_ ), 
	limit_( limit ), result_( Logdom::null() )
{
	if( limit <= Logdom::one() ) {
		mins_[0] = maxs_[0] = 0 ;
		put( 0, 0, 0, 0, 0, 0, Logdom::one() ) ;
		for( size_t y = 0 ; y != width_ && mins_[y] != maxs_[y] ; ++y )
		{
			assert( y <= width_ ) ;
			assert( mins_[y] >= 0 ) ;
			assert( maxs_[y] <= width_ ) ;
			assert( mins_[y] <= maxs_[y] ) ;

			mins_[y+1] = maxs_[y+1] = 0 ;
			// std::cerr << y << std::endl ;
			// expand the current row for each state in turn... of course,
			// each state is a special case.
			for( uint32_t x = mins_[y] ; x != width_-1 && x != maxs_[y] ; ++x )
			{
				for( size_t s = 0 ; s != adna_parblock::num_states ; ++s )
				{
					assert( width_*y + x < (1+width_) * width_ ) ;
					assert( width_*y + x < cells_.size() ) ;

					Logdom score = cells_[ width_*y + x ][ s ].score ;
					if( score.is_finite() ) extend( pb, score, s, x, y, reference+x, query+y ) ;
				}
			}
		}
	}
}

// what to do?  
// If in matching state, we know there's no immediate match, so we can...
// - mismatch
// - detect deamination and change to SS state (while matching)
// - open ref gap
// - open query gap
//
// If a gap is open, we can...
// - extend it
// - close it
//
// If we hit a gap symbol, we must...
// - start over at second half in initial state

void ExtendAlignment::extend( const adna_parblock &pb_, Logdom score, int s, int x, int y, DnaP ref, const QSequence::Base *qry )
{
	// Note the penalties: The appropriate substitution penalty is
	// applied whenever we (mis-)match two codes, the gap open penalties
	// are applied when opening/extending a gap, the
	// overhang_enter_penalty is applied when changing to SS mode and
	// the overhang_ext_penalty is applied whenever moving along the
	// query while single stranded, even when a gap is open!  This gives
	// correct scores for a geometric distribution of overhang lengths.
	if( !*ref || !qry->ambicode ) 
	{
		// We hit a gap in either the reference or the query.  Whatever
		// is left of the query (if any) must be penalized.  To do this,
		// we virtually extend the reference with Ns and align to those.
		// This is a white lie in that it will overestimate the real
		// penalty, but that's okay, because such an alignment isn't all
		// that interesting in reality anyway.  Afterwards we're
		// finished and adjust result and limit accordingly.
		for( int yy = 0 ; qry[ yy ].ambicode ; ++yy )
		{
			score *= pb_.subst_penalty( s, *ref, qry[ yy ] ) ;
			if( s & adna_parblock::mask_ss ) score *= pb_.overhang_ext_penalty ;
		}
		if( score > result_ ) {
			result_ = limit_ = score ;
			max_s_ = s ;
			max_x_ = x ;
			max_y_ = y ;
		}
	}
	else if( (s & adna_parblock::mask_gaps) == 0 )
	{
		// no gaps open --> mismatch, open either gap, enter SS
		put( s, s, x, 1, y, 1, score * pb_.subst_penalty( s, *ref, *qry ) 
				* ( s & adna_parblock::mask_ss ? pb_.overhang_ext_penalty : Logdom::one() ) ) ;
		if( *ref != qry->ambicode ) {	// only on a mismatch try anything fancy
			put( s | adna_parblock::mask_gap_qry, s, x, 1, y, 0, score * pb_.gap_open_penalty ) ;
			put( s | adna_parblock::mask_gap_ref, s, x, 0, y, 1, score * pb_.gap_open_penalty ) ;
			if( pb_.overhang_enter_penalty.is_finite() && (s & adna_parblock::mask_ss) == 0 )
			{
				// To enter single stranded we require that the penalty for
				// doing so is immediately recovered by the better match.
				// This is easily the case for the observed deamination
				// rates in aDNA.
				Logdom p0 = pb_.subst_penalty( s, *ref, *qry ) ;
				Logdom p4 = pb_.subst_penalty( s | adna_parblock::mask_ss, *ref, *qry ) 
					* pb_.overhang_enter_penalty * pb_.overhang_ext_penalty ;
				if( p4 > p0 ) put( s | adna_parblock::mask_ss, s, x, 1, y, 1, score * p4 ) ;
			}
		}
	}
	else if( (s & adna_parblock::mask_gaps) == adna_parblock::mask_gap_ref )
	{
		put( s, s, x, 0, y, 1, score * pb_.gap_ext_penalty *
				( s & adna_parblock::mask_ss ? pb_.overhang_ext_penalty : Logdom::one() ) ) ;
		put( s & ~adna_parblock::mask_gap_ref, s, x, 1, y, 1, score * pb_.subst_penalty( s, *ref, *qry ) *
				( s & adna_parblock::mask_ss ? pb_.overhang_ext_penalty : Logdom::one() ) ) ;
	}
	else
	{
		put( s, s, x, 1, y, 0, score * pb_.gap_ext_penalty ) ;
		put( s & ~adna_parblock::mask_gap_qry, s, x, 1, y, 1, score * pb_.subst_penalty( s, *ref, *qry ) *
				( s & adna_parblock::mask_ss ? pb_.overhang_ext_penalty : Logdom::one() ) ) ;
	}
}

void ExtendAlignment::backtrace( std::vector<unsigned>& out ) const
{
	for( size_t x = max_x_, y = max_y_, s = max_s_ ; x || y ; )
	{
		const Cell& c = cells_[ width_*y + x ][ s ] ;
		if( !c.from_x_offset && !c.from_y_offset ) throw "stuck in backtracing" ;
		else if( !c.from_x_offset ) streams::push_i( out, c.from_y_offset ) ;
		else if( !c.from_y_offset ) streams::push_d( out, c.from_x_offset ) ;
		else if( c.from_y_offset == c.from_x_offset ) streams::push_m( out, c.from_x_offset ) ;
		else throw "inconsistency in backtracing" ;

		x -= c.from_x_offset ;
		y -= c.from_y_offset ;
		s = c.from_state ;
	}
}

// Extension of both sides.  We first run a forward extension at half
// the limit.  If this succeeds, we do the backwards extension limited
// to whatever is left.  If forward extension fails, we do backwards
// extension to half the limit, then add forward using up what's left.
//
// Only one alignment is produced, but we make sure it is the cheapest
// one.  The score may exceed the limit, if we happen to finish right
// when stepping over the limit.  If really nothing is found, we return
// Logdom::null().

ExtendBothEnds::ExtendBothEnds(
		const adna_parblock& pb,
		const QSequence& query,
		const SeededAlignment& seed,
		Logdom limit ) :
	forwards_(
			pb,
			seed.reference_ + seed.size_,
			query.start() + seed.qoffs_ + seed.size_,
			( limit / seed.score_ ).sqrt() ),
	backwards_(
			pb,
			seed.reference_.reverse() + 1,
			query.start() - seed.qoffs_ + 1,
			limit / forwards_.get_result() / seed.score_ ),
	score_( forwards_.get_result() * backwards_.get_result() * seed.score_ )
{
	if( score_.is_finite() ) return ;

	ExtendAlignment backwards2(
			pb,
			seed.reference_.reverse() + 1,
			query.start() - seed.qoffs_ + 1,
			(limit / seed.score_).sqrt() ) ;
	ExtendAlignment forwards2(
			pb,
			seed.reference_ + seed.size_,
			query.start() + seed.qoffs_ + seed.size_,
			limit / backwards2.get_result() / seed.score_ ) ;
	Logdom score2 = forwards2.get_result() * backwards2.get_result() * seed.score_ ;
	if( score2.is_finite() )
	{
		score_ = score2 ;
		forwards_.swap( forwards2 ) ;
		backwards_.swap( backwards2 ) ;
	}
}

std::vector<unsigned> ExtendBothEnds::backtrace( const SeededAlignment& seed, DnaP &minpos, DnaP &maxpos ) const 
{
	minpos = seed.reference_ - backwards_.max_x() ;
	maxpos = seed.reference_ + seed.size_ + forwards_.max_x() ;

	std::vector<unsigned> trace ;
	backwards_.backtrace( trace ) ;

	trace.push_back( 0 ) ;
	streams::push_m( trace, seed.size_ ) ;
	trace.push_back( 0 ) ;

	std::vector<unsigned> rtrace ;
	forwards_.backtrace( rtrace ) ;
	std::copy( rtrace.rbegin(), rtrace.rend(), back_inserter( trace ) ) ;
	return trace ;
}

