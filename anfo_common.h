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
