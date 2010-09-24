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

/*! \page shotgun sequencing simulation
 *
 * Reads a DNA file, breaks out sequences of a certain length
 * distribution, applies a given level of mutations and sequencing
 * error, and spits out the result as FastQ.
 */

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include "index.h"
#include "util.h"
#include "config.pb.h"

#include <google/protobuf/text_format.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>

#include <popt.h>

#include <iostream>
#include <iomanip>
#include <sstream>

#include <fcntl.h>
#include <unistd.h>

int norm( double m, double s )
{
	return m + s * sqrt(2) * erfc( 2*drand48() - 1 ) ;
}

//! \page simulated shotgun sequencing
//!
//! Simulates shotgun sequencing with normally distributed fragment
//! length, constant mutation and error rate, and affine gap
//! probabilities.
//!
//! \todo Simulate somewhat variable quality and error rate.

WRAPPED_MAIN
{
	GOOGLE_PROTOBUF_VERIFY_VERSION ;
	enum option_tags { opt_none, opt_version } ;
	const char* genome_name = 0 ;
	int verbose = 0 ;
	int amount = 100000 ;
	int mean_length = 50 ;				// Neanderthal?
	double length_dev = 15 ;
	int quality = 40 ;					// rather typical
	double mutation_rate = 0.014 ;		// approx. human/chimp divergence
	double gap_rate = 0.001 ; 			// pulled out of my ass
	double gap_ext_rate = 0.01 ; 		// likewise

	config::Config mi ;
	struct poptOption options[] = {
		{ "version",     'V', POPT_ARG_NONE,   0,            opt_version, "Print version number and exit", 0 },
		{ "verbose",     'v', POPT_ARG_NONE,   &verbose,     opt_none,    "Make more noise while working", 0 },
		{ "genome",      'g', POPT_ARG_STRING, &genome_name, opt_none,    "Read genome GENOME", "GENOME" },
		{ "amount",      'n', POPT_ARG_INT,    &amount,      opt_none,    "Generate N sequences", "N" },
		{ "quality",     'q', POPT_ARG_INT,    &quality,     opt_none,    "Set sequencing quality to Q", "Q" },
		{ "mutations",   'm', POPT_ARG_DOUBLE, &mutation_rate, opt_none,  "Set mutation (SNP) rate to R", "R" },
		{ "gap-opens",   'O', POPT_ARG_DOUBLE, &gap_rate,    opt_none,    "Set gap open rate to R", "R" },
		{ "gap-extensions",'E',POPT_ARG_DOUBLE,&gap_ext_rate,opt_none,    "Set gap extension rate to R", "R" },
		{ "length",      'l', POPT_ARG_INT,    &mean_length, opt_none,    "Set desired length to L", "L" },
		{ "sigma",       's', POPT_ARG_DOUBLE, &length_dev,  opt_none,    "Set standard deviation of length to S","S"},
		POPT_AUTOHELP POPT_TABLEEND
	} ;

	poptContext pc = poptGetContext( "random_sample", argc, argv, options, 0 ) ;
	poptSetOtherOptionHelp( pc, "[OPTION...]" ) ;
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

	if( !genome_name ) throw "missing --genome option" ;
	if( poptGetArg( pc ) ) throw "unexpected non-option argument" ;

	GenomeHolder genome = Metagenome::find_genome( genome_name ) ;
	for( int i = 0 ; i != amount ; )
	{
		int length = norm( mean_length, length_dev ) ;
		if( length > 0 ) {
			uint32_t start = (genome->raw_size() - length) * drand48() ;
			DnaP ps = genome->get_base() + start, p = ps, pe = p + length ;
			while( p != pe && *p ) ++p ;
			if( p == pe ) 
			{
				uint32_t offset ;
				const config::Sequence *seq = genome->translate_back( ps, offset ) ;

				if( rand() % 1 ) 
				{
					ps = pe.reverse()+1 ;
					pe = ps + length ;
					length = -length ;
				}

				std::cout << '@' << seq->name() << ':' << offset 
					<< "+-"[ length < 0 ] << length << '\n' ;
				std::string qs ;

				for( p = ps ; p != pe ; )
				{
					double pp = drand48() ;
					if( pp < gap_rate ) // make a deletion
					{
						do
						{
							++p ;
						} while( p != pe && drand48() < gap_ext_rate ) ;
					}
					else {
						if( pp < 2*gap_rate ) // make an insertion
						{
							do {
								std::cout << (1 << (rand()&3)) ;
								qs.push_back( 33+quality ) ;
							} while( drand48() < gap_ext_rate ) ;
						}

						// mutation and/or seq. error?
						Ambicode b = drand48() < mutation_rate + pow(10,-quality/10) ? 1 << (rand()&3) : *p ;
						std::cout << from_ambicode(b) ;
						qs.push_back( 33+quality ) ;
						++p ;
					}
				}
				std::cout << "\n+\n" << qs << '\n' ;
				++i ;
			}
		}
	}

	return 0 ;
}


