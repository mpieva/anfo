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

#ifndef INCLUDED_ALIGN_FWD_H
#define INCLUDED_ALIGN_FWD_H

#include "logdom.h"
#include "sequence.h"

//! \page adna_alignment Operations on Ancient DNA Alignments
//! This type of alignment models aDNA as two additional states, one
//! with higher deamination rate, with an asymmetric substitution
//! matrix, and with affine gap costs.  Additional parameters are all
//! static.
//!
//! \todo We actually don't need the parameter set to be static, and it
//!       may even pay to have it configurable so differently
//!       parameterized alignments an be mixed in a single run.
//!
//! @{

//! \brief substitution matrix
//! We prepare a full matrix of 16 ambiguity codes vs. 16 ambiguity
//! codes, this avoid the need for expensive additions in the log-domain
//! should the need to align ambiguity codes arise.  First index is
//! "from" (reference code), second index is "to" (query code).
typedef Logdom subst_mat[16][16] ;

namespace config { class Aligner ; } ;

//! \brief parameters for simple_adna
//! One such parblock is a static variable for the aligner proper,
//! various support tools shunt additional structures around.
struct adna_parblock
{
	enum {
		mask_ss      = 1,
		mask_gap_ref = 2,
		mask_gap_qry = 4,

		mask_gaps = mask_gap_ref | mask_gap_ref,
		num_states = 6
	} ;

	adna_parblock() {} 
	adna_parblock( const config::Aligner& conf ) ;

	//! \brief DS substitution matrix, forward direction
	subst_mat ds_mat ;

	//! \brief SS substitution matrix, forward direction
	//! Deamination shows up as C->T as it is best understood this way.
	//! To process reverse-complemented deamination, we have to do
	//! rev-complemented lookups while actually moving in the forward
	//! (5'->3') direction.  \see simple_adna::subst_penalty()
	subst_mat ss_mat ;

	//! \brief Penalty for extending an overhang.
	//! Having a constant penalty for the overhang length models its
	//! length distribution as geometric.
	Logdom overhang_ext_penalty ;

	//! \brief penalty for entering SS state
	//! This is essentially the probability of having an overhang at
	//! all.
	Logdom overhang_enter_penalty ;

	//! \brief gap open penalty
	Logdom gap_open_penalty ;

	//! \brief gap extension penalty
	Logdom gap_ext_penalty ;


	Logdom subst_penalty( int s, Ambicode r, const QSequence::Base &qry ) const
	{
		Logdom prob = Logdom::null() ;
		for( uint8_t p = 0 ; p != 4 ; ++p )
			prob += ( s & mask_ss ? ss_mat[r][1<<p] : ds_mat[r][1<<p] )
				* Logdom::from_phred( qry.qscores[p] ) ;
		return prob ;
	}
} ;

std::ostream& operator << ( std::ostream&, const adna_parblock& ) ;

#endif
