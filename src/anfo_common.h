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
#include "stream.h"

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

namespace streams {

//! \brief takes care of stuff needed before aligning
//! This is for "housekeeping" that needs to be done before indexing
//! and alignment, most notably adapter/artefact trimming.

class Housekeeper : public Stream
{
	private:
		const int minscore_ ;
		vector< string > trim_left_ ;
		vector< string > trim_right_ ;

	public:
		Housekeeper( const config::Config &config ) :
			minscore_( config.min_trim_score() ),
			trim_left_( config.trim_left().begin(), config.trim_left().end() ),
			trim_right_( config.trim_right().begin(), config.trim_right().end() )
		{}

		Housekeeper( vector<string> left, vector<string> right, int minscore )
			: minscore_( minscore_ ), trim_left_( left ), trim_right_( right ) {}

		virtual void put_result( const Result& ) ;
} ;

//! \brief configured Indexer
//! Knows about one index, will look at each record, check the policy
//! and if appropriate, will do an index lookup.  Constructed seeds are
//! stored and passed on.  The policy is taken from a config file and
//! stored in the header.
class Indexer : public Stream
{
	private:
		config::Config conf_ ;
		std::string index_name_ ;
		const FixedIndex& index_ ;

	public:
		Indexer( const config::Config &config, const std::string& index_name ) ;
		virtual ~Indexer() ;

		virtual void put_header( const Header& ) ;
		virtual void put_result( const Result& ) ;
} ;

//! \brief configured mapper
//! Knows about one genome and an aligner configuration.  Looks at each
//! record, takes the seed collection for its genome and aligns.  Seeds
//! are removed, the result is added and passed on.
class Mapper : public Stream
{
	private:
		config::Aligner conf_ ;
		GenomeHolder genome_ ;

	public:
		Mapper( const config::Aligner &config, const std::string& genome_name ) ;

		virtual void put_header( const Header& ) ;
		virtual void put_result( const Result& ) ;
} ;

} ; // namespace streams

std::string expand( const std::string&, int ) ;

#endif
