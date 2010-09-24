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

//! \page index_test Testing the indexing phase
//! This binary just runs the indexer on whatever input it is given,
//! printing the resulting seeds in various detail.  If the input looks
//! as if it contains a coordinate (simulated data does), this program
//! will check whether a seed in the correct region is found at all.

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include "align.h"
#include "anfo_common.h"
#include "compress_stream.h"
#include "conffile.h"
#include "index.h"
#include "stream.h"
#include "util.h"

#include <popt.h>

#include <limits>
#include <string>
#include <iostream>
#include <iomanip>

using namespace config ;
using namespace std ;
using namespace google::protobuf::io ;

struct Subject
{
	string chrom, seq ;
	int start, len ;
	bool strand ;
} ;

istream& operator >> ( istream& s, Subject& t )
{
	string strand ;
	s >> t.chrom >> strand >> t.start >> t.len >> t.seq ;
	t.strand = !strand.empty() && strand[0] == '-' ;
	return s ; 
}

WRAPPED_MAIN
{
	GOOGLE_PROTOBUF_VERIFY_VERSION ;
	enum option_tags { opt_none, opt_version } ;

	FixedIndex::LookupParams params ;
	params.cutoff = numeric_limits<uint32_t>::max() ;
	params.allow_mismatches = 0 ;
	int lower_limit = 9 ;

	struct poptOption options[] = {
		{ "version",     'V', POPT_ARG_NONE,   0,            opt_version, "Print version number and exit", 0 },
		{ "lower", 0, POPT_ARG_INT, &lower_limit, 0, "Set lower limit for word length", 0 },
		{ "mismatches", 0, POPT_ARG_INT, &params.allow_mismatches, 0, "Set number of allowed mismatches in seed", 0 },
		POPT_AUTOHELP POPT_TABLEEND
	} ;

	console.loglevel = Console::info ;

	poptContext pc = poptGetContext( "anfo", argc, argv, options, 0 ) ;
	poptSetOtherOptionHelp( pc, "[OPTION...] [sequence-file...]" ) ;
	poptReadDefaultConfig( pc, 0 ) ;

	if( argc <= 1 ) { poptPrintHelp( pc, stderr, 0 ) ; return 1 ; }
	for( int rc = poptGetNextOpt( pc ) ; rc > 0 ; rc = poptGetNextOpt(pc) ) switch( rc )
	{
		case opt_version:
			std::cout << poptGetInvocationName(pc) << ", revision " << PACKAGE_VERSION << std::endl ;
			return 0 ;

		default:
			std::clog << poptGetInvocationName(pc) << ": " << poptStrerror( rc ) 
				<< ' ' << poptBadOption( pc, 0 ) << std::endl ;
			return 1 ; 
	}

	FixedIndex index( "hg18_chr1_12.idx" ) ;
	GenomeHolder genome = Metagenome::find_genome( index.metadata().genome_name() ) ;

	vector<Subject> subjects ;
	while( const char* arg = poptGetArg( pc ) ) 
	{
		ifstream ifs( arg ) ;
		Subject subj ;
		while( ifs >> subj ) subjects.push_back( subj ) ;
	}
		
	for( params.wordsize = lower_limit ; params.wordsize <= 12 ; ++params.wordsize )
	{
		for( params.stride = 4 ; params.stride <= 8 ; params.stride *= 2 )
		{
			for( int min_seed_len = 24 ; min_seed_len <= 48 ; min_seed_len+=4 )
			{
				cout << params.wordsize << '\t' << params.stride << '\t' << min_seed_len << '\t' << flush ;

				unsigned total_seeds = 0, total_seqs = 0, total_hits = 0 ;

				for( vector<Subject>::const_iterator subj = subjects.begin() ; subj != subjects.end() ; ++subj )
				{
					PreSeeds seeds ;
					int num_useless ;
					index.lookupS( subj->seq, seeds, params, &num_useless ) ;

					output::Seeds ss ;
					total_seeds += combine_seeds( seeds, min_seed_len, &ss ) ;
					// XXX total_seeds += params.allow_mismatches 
						// XXX ? combine_seeds( seeds, min_seed_len, &ss ) 
						// XXX : select_seeds( seeds, 2 /*p.max_diag_skew()*/, 4 /*p.max_gap()*/, min_seed_len, index.gaps(), &ss ) ;
					++total_seqs ;

					unsigned left = genome->find_pos( subj->chrom, subj->start ) - genome->get_base() ;
					for( int i = 0 ; i != ss.ref_positions_size() ; ++i )
					{
						unsigned ref = ss.ref_positions(i) ;
						int qry = ss.query_positions(i) ;
						if( subj->strand == (qry < 0) && ref >= left && ref < left+subj->len )
						{
							++total_hits ;
							break ;
						}
					}
				}
				cout << total_seqs << '\t' << total_seeds << '\t' << total_hits << endl ;
			}
		}
	}

	// cout << endl ;

	poptFreeContext( pc ) ;
	return 0 ;
}
