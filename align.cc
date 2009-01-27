#include "align.h"

#include <cmath>
// #include <iostream>

subst_mat simple_adna::ds_mat ;
subst_mat simple_adna::ss_mat ;
uint32_t simple_adna::overhang_ext_penalty ;
uint32_t simple_adna::overhang_enter_penalty ;
uint32_t simple_adna::gap_open_penalty ;
uint32_t simple_adna::gap_ext_penalty ;

uint32_t to_log_dom( double p ) { return (uint32_t)( -10 * std::log( p ) + 0.5 ) ; }

//! \brief sets up aDNA alignments according to parameter set
//! Penalties are logs of probabilities.  We calculate a probability in
//! a straightforward manner, then take the log, then multiply by (say)
//! 10 and round to an integer (effectively shifting the base).
//! Afterwards we could shift the score up so the minimum score is zero,
//! but it will be zero to within accuracy anyway.
//!
//! For ambiguity codes, we sum over all the possible codes, divide by
//! the number of possibilities, the go to the log domain.  Note that
//! matches end up costing nothing, matching to an N *does* cost
//! something, and matching a T to a C might be dirt cheap.
//!
//! For gaps, we already have probabilities.  Taking the log is all it
//! takes.
//!
//! For the overhang, we want something that's easy to parameterize.  We
//! assume a geometric distribution of overhang length with parameter p,
//! but half the overhangs end up pointing in the wrong direction, which
//! maps them to zero length.  The mean overhang length ends up being
//! (1-p)/(2p), and we take that as parameter.
void simple_adna::configure( const config::Aligner& conf ) 
{
	gap_ext_penalty = to_log_dom( conf.gap_extension_rate() ) ;
	gap_open_penalty = conf.has_gap_open_rate() ? to_log_dom( conf.gap_open_rate() ) : gap_ext_penalty ;

	if( conf.has_mean_overhang_length() )
	{
		double p_comp = 1 - 1 / (2*conf.mean_overhang_length() + 1) ;
		overhang_enter_penalty = to_log_dom( p_comp / 2 ) ;
		overhang_ext_penalty = to_log_dom( p_comp ) ;
		// std::cerr << p_comp << ' ' << overhang_ext_penalty << ' ' << overhang_enter_penalty << '\n' << '\n' ;
	}
	else
	{
		overhang_enter_penalty = ~0U ;
		overhang_ext_penalty = ~0U ;
	}

	// std::cerr << gap_ext_penalty << ' ' << gap_open_penalty << '\n' 

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

#if 0
	for( int i = 0 ; i != 4 ; ++i )
	{
		for( int j = 0 ; j != 4 ; ++j )
		{
			std::cerr << to_log_dom(prim_mat_ds[i][j]) << " (" << to_log_dom(prim_mat_ss[i][j]) << ") " ;
		}
		std::cerr << '\n' ;
	}
	std::cerr << '\n' ;
#endif

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

			ds_mat[ref][qry] = to_log_dom( p_ds / npairs ) ;
			ss_mat[ref][qry] = to_log_dom( p_ss / npairs ) ;

			// std::cerr << ds_mat[ref][qry] << " (" << ss_mat[ref][qry] << ")\t" ;
		}
		// std::cerr << '\n' ;
	}
	// std::cerr << '\n' ;
	// exit(0) ;
}

