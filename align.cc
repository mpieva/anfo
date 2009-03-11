#include "align.h"

#include <ostream>
#include <iomanip>

subst_mat simple_adna::ds_mat ;
subst_mat simple_adna::ss_mat ;
uint32_t simple_adna::overhang_ext_penalty ;
uint32_t simple_adna::overhang_enter_penalty ;
uint32_t simple_adna::gap_open_penalty ;
uint32_t simple_adna::gap_ext_penalty ;

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
				s << std::setw(2) << to_log_dom( mat[ref][qry] ) << "  " ;
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
//! but accuracy isn't good enough to pick that up), matching to an N
//! *does* cost something, and matching a T to a C might be dirt cheap.
//!
//! For gaps, we already have probabilities.  Taking the log is all it
//! takes to make scores from that.
//!
//! For the overhang, we want something that's easy to parameterize.  We
//! assume a geometric distribution of overhang length with parameter p,
//! but half the overhangs end up pointing in the wrong direction, which
//! maps them to zero length.  The mean overhang length ends up being
//! (1-p)/(2p), and we take that as parameter.  Entering an overhang
//! involves a penalty (since the probability is less than one), but so
//! does not entering it.  We simply take the difference of the two,
//! leaving the penalty for a perfect alignment at zero.
//!
//! \param conf set of configuration parameters
//! \param out if not NULL, receives a protocol of the calculated
//!            internal parameters
void simple_adna::configure( const config::Aligner& conf, std::ostream *out ) 
{
	gap_ext_penalty = to_log_dom( conf.gap_extension_rate() ) ;
	gap_open_penalty = conf.has_gap_open_rate() ? to_log_dom( conf.gap_open_rate() ) : gap_ext_penalty ;

	if( conf.has_mean_overhang_length() )
	{
		double p_comp = 1 - 1 / (2*conf.mean_overhang_length() + 1) ;
		overhang_enter_penalty = to_log_dom( p_comp / 2 ) - to_log_dom( 1 - p_comp / 2 ) ;
		overhang_ext_penalty = to_log_dom( p_comp ) ;
	}
	else
	{
		overhang_enter_penalty = ~0U ;
		overhang_ext_penalty = ~0U ;
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

			ds_mat[ref][qry] = p_ds / npairs ;
			ss_mat[ref][qry] = p_ss / npairs ;
		}
	}
	
	if( out ) 
	{
		*out << "\n\e[1mALIGNMENT PARAMETERS\e[0m:\n\n"
			 << "  \e[4mDouble Stranded Substitution Matrix\e[0m:\n\n" ;
		print_subst_mat( *out, ds_mat ) ;
		*out << "\n  \e[4mSingle Stranded Substitution Matrix\e[0m:" ;
		if( conf.has_rate_of_ss_deamination() ) {
			*out << "\n\n" ; print_subst_mat( *out, ss_mat ) ;
		} else *out << " N/A\n" ;
		*out << "\n  \e[4mGap Open Penalty\e[0m: " << gap_open_penalty
			 << "\n  \e[4mGap Extension Penalty\e[0m: " << gap_ext_penalty
			 << "\n\n  \e[4mOverhang Penalty\e[0m: " ;
		if( overhang_enter_penalty == ~0U ) *out << "N/A" ; else *out << overhang_enter_penalty ;
		*out << "\n  \e[4mOverhang Extension Penalty\e[0m: " ;
		if( overhang_ext_penalty == ~0U ) *out << "N/A" ; else *out << overhang_ext_penalty ;
		*out << '\n' << std::endl ;
	}
}

