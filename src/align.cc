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

