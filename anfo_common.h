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

#ifndef INCLUDED_ANFO_COMMON_H
#define INCLUDED_ANFO_COMMON_H

// Commonalities between different command line utilities, factored out.

#include "align.h"
#include "index.h"

#include "config.pb.h"
#include "output.pb.h"

#include <map>
#include <string>

#include <fnmatch.h>

//! \brief Configures type of alignment to be done.
//! This is a hack and only intended as a stopgap.
// typedef flat_alignment alignment_type ;
typedef simple_adna alignment_type ;

struct reference_overlaps {
	DnaP x, y ;
	reference_overlaps( DnaP u, DnaP v ) : x(u), y(v) {}
	bool operator()( const alignment_type& a ) {
		return a.reference >= x && a.reference <= y ; }
} ;

//! \brief configured mapper
//! This is what you would think of as an instance of ANFO running.
//! Wrapping it in a class makes dragging along all the configuration
//! data somewhat easier.
class Mapper
{
	private:
		config::Config mi ;
		// Genomes genomes ;
		Indices indices ;

	public:
		Mapper( const config::Config &config ) ;

		int index_sequence( output::Result &r, QSequence &qs, std::deque< alignment_type >& ol ) ;
		void process_sequence( const QSequence &ps, double max_penalty_per_nuc, std::deque< alignment_type > &ol, output::Result &r ) ;
} ;

#endif
