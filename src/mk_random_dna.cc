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

#include "conffile.h"
#include "config.pb.h"
#include "index.h"
#include "sequence.h"
#include "util.h"

#include <popt.h>

#include <cctype>
#include <fstream>
#include <iostream>

//! \brief decodes a FASTA file
//!
//! This class is a state machine that recognizes a FASTA file.
//! Nucleotides are split into contigs and written out to an \c ostream,
//! names and coordinates are remembered and can later be written to a
//! config file.
//!
//! \todo The state machine is implemented using a pointer to member
//!       functions.  This is straight-forward, but also a bit hard to
//!       read and it may be slower than necessary.  This could be
//!       enhanced by including a loop in each of the involved
//!       functions.

class FastaFake
{
	private:
		unsigned position_ ;
		unsigned source_position_ ;
		uint8_t one_nt_ ;
		bool verbose_ ;

		std::ostream& dna_ ;
		config::Genome cur_genome_ ;
		config::Sequence *cur_sequence_ ;
		config::Contig *cur_contig_ ;
		
		
	public:
		//! \brief sets up a fresh decoder
		//! The decoder is initialized with an output stream, a genome
		//! metadata structure and some parameters.  The output stream
		//! is assumed to be empty, which means that a signature is
		//! written first and all coordinates will be relative to the
		//! current position.  The metadata should be properly
		//! initialized, if it contains sequence definitions, they will
		//! be lost.
		//!
		//! \param dna output stream for compact DNA data
		//! \param genome metadata about genome, will be filled in
		//! \param maxn maximum number of Ns in a row that is not a gap
		//! \param verbose whether to produce progress reports to \c
		//!        stderr
		FastaFake( std::ostream& dna, const config::Genome &genome, bool verbose = false )
			: position_( 24 ), source_position_( 0 ), one_nt_( 0 )
			, verbose_( verbose ), dna_( dna )
			, cur_genome_( genome ), cur_sequence_( 0 ), cur_contig_( 0 )
		{
			uint32_t dummy = 0, signature = CompactGenome::signature ;
			dna_.write( (const char*)&signature, 4 ) ;
			dna_.write( (const char*)&dummy, 4 ) ;
			dna_.write( (const char*)&dummy, 4 ) ;
			cur_genome_.clear_sequence() ;
		}

		void fake( int length )
		{
			cur_sequence_ = cur_genome_.add_sequence() ;
			*cur_sequence_->mutable_name() = cur_genome_.name() ;
			put_nt( '-' ) ;
			cur_contig_ = cur_sequence_->add_contig() ;
			cur_contig_->set_offset( position_ ) ;
			cur_contig_->set_range_start( source_position_ ) ;
			report() ;

			for( int i = 0 ; i != length ; ++i, ++source_position_ )
				put_nt( "ACGT"[ rand() & 3 ] ) ;

			cur_contig_->set_range_end( source_position_ ) ;
			report() ;
			put_nt( '-' ) ;
			cur_genome_.set_total_size( position_ ) ;
			if( position_ & 1 ) put_nt( '-' ) ;
			if( verbose_ ) std::clog << "\33[K" << std::flush ;

			uint32_t ix_start = dna_.tellp() ;
			if( !cur_genome_.SerializeToOstream( &dna_ ) ) 
				throw "couldn't serialize config" ;
			uint32_t ix_len = (uint32_t)dna_.tellp() - ix_start ;
			dna_.seekp( 4 ) ;
			dna_.write( (char*)&ix_start, 4 ) ;
			dna_.write( (char*)&ix_len, 4 ) ; 
		}

		void report() const {
			if( verbose_ )
				std::clog << cur_sequence_->name() << " @ " << source_position_
					      << "\33[K\r" << std::flush ;
		}

		void put_nt( char c )
		{
			if( position_ & 1 ) dna_.put( (to_ambicode(c) << 4) | one_nt_ ) ;
			else                one_nt_ = to_ambicode(c) ;
			++position_ ;
			if( !(position_ & 0xfffff ) ) report() ;
		}
} ;

std::string drop_suffix( const std::string& suf, const std::string& s )
{
	if( s.substr( s.size()-suf.size() ) == suf ) return s.substr( 0, s.size()-suf.size() ) ;
	else return s ;
}

WRAPPED_MAIN
{
	GOOGLE_PROTOBUF_VERIFY_VERSION ;
	enum option_tags { opt_none, opt_version } ;

	const char* output_file = 0 ;
	const char* description = 0 ;
	const char* genome_name = 0 ;
	int verbose = 0 ;
	int genome_size = 0 ;

	struct poptOption options[] = {
		{ "version",     'V', POPT_ARG_NONE,   0,            opt_version, "Print version number and exit", 0 },
		{ "output",      'o', POPT_ARG_STRING, &output_file, opt_none,    "Write DNA output to FILE", "FILE" },
		{ "genome",      'g', POPT_ARG_STRING, &genome_name, opt_none,    "Set genome name to NAME", "NAME" },
		{ "size",        's', POPT_ARG_INT,    &genome_size, opt_none,    "Set genome size to N", "N"},
		{ "description", 'd', POPT_ARG_STRING, &description, opt_none,    "Add TEXT as description to genome", "TEXT" },
		{ "verbose",     'v', POPT_ARG_NONE,   &verbose,     opt_none,    "Make more noise while working", 0 },
		POPT_AUTOHELP POPT_TABLEEND
	} ;

	poptContext pc = poptGetContext( "mk_random_dna", argc, argv, options, 0 ) ;
	poptSetOtherOptionHelp( pc, "[OPTION...]" ) ;
	poptReadDefaultConfig( pc, 0 ) ;

	if( argc <= 1 ) { poptPrintHelp( pc, stderr, 0 ) ; return 1 ; }
	for( int rc = poptGetNextOpt( pc ) ; rc > 0 ; rc = poptGetNextOpt(pc) ) switch( rc )
	{
		case opt_version:
			std::cout << poptGetInvocationName(pc) << ", revision " << PACKAGE_VERSION << std::endl ;
			return 0 ;

		default:
			std::clog << argv[0] << ": " << poptStrerror( rc ) 
				<< ' ' << poptBadOption( pc, 0 ) << std::endl ;
			return 1 ; 
	}

	if( !genome_name ) throw "missing --genome option" ;
	config::Genome g ;

	if( description ) g.set_description( description ) ;
	g.set_name( genome_name ) ;

	std::ofstream output_stream( 
			(output_file ? output_file : (g.name() + ".dna")).c_str(),
			std::ios::trunc ) ;
	FastaFake fd( output_stream, g, verbose ) ;
	fd.fake( genome_size ) ;
	return 0 ;
}


