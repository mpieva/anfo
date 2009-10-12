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

/*!
 * \page fasta_to_dna FASTA to DNA file converter
 *
 * \c fa2dna reads a set of FASTA files with one ore more DNA or RNA
 * sequences, possibly containing ambiguity codes, from files or stdin,
 * writes it in a compact representation into a new file and writes meta
 * information into another file.
 *
 * We chop the genome into contigs.  Every sequence start also starts a
 * new contig, stretches of more than 2 Ns (configurable, but more Ns
 * become nonsensical quite quickly) separate contigs and are not
 * stored.  Ambiguity codes are recognized, and so is 'U' for uracil.
 *
 * A "compact genome" is a file encoding DNA as 4 bits per nucleotide.
 * The bits code for A,C,T,G from LSB to MSB.  Combinations of bits form
 * ambiguity codes (0 is a gap, used as terminator only, 15 is an N),
 * two left rotations therefore reverse-complement a single nucleotide.
 * In a byte, the first nucleotide occupies the 4 LSBs.
 *
 * The "compact genome" starts with the ACSII-signature "DNA1", then the
 * offset to the internal metadata, then the length of the metadata, followed
 * by the raw DNA.  Finally the metadata is attached, in the form of a
 * binary protocol buffer message of type config::Genome.  
 *
 * In the DNA part, each contig is prepended with a single gap, then the
 * contigs are concatenated, then a final single gap is appended.  The
 * DNA is written out with two nucleotides per byte, and the final byte
 * is padded with a gap if necessary.  Note that every contig is always
 * preceded and succeeded by a gap; this is done so we can \c mmap a
 * genome into memory and operate on pointers into it in the knowledge
 * that every contig is terminated by a gap in either direction.
 *
 * \todo A genome may or may not contain information about taxonomy, and
 *       this info may apply to the genome or just some sequences.  Add
 *       support to either set a taxid or look it up in some sort of
 *       data base.
 */

#include "config.h"
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

class FastaDecoder
{
	private:
		void (FastaDecoder::*state_)( char ) ;

		unsigned position_ ;
		unsigned source_position_ ;
		unsigned num_n_  ;
		uint8_t one_nt_ ;
		const unsigned max_n_ ;
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
		FastaDecoder( std::ostream& dna, const config::Genome &genome, unsigned maxn = 2, bool verbose = false )
			: state_( &FastaDecoder::s_idle ), position_( 24 ), num_n_( 0 ), one_nt_( 0 )
			, max_n_( maxn ), verbose_( verbose ), dna_( dna )
			, cur_genome_( genome ), cur_sequence_( 0 ), cur_contig_( 0 )
		{
			uint32_t dummy = 0, signature = CompactGenome::signature ;
			dna_.write( (const char*)&signature, 4 ) ;
			dna_.write( (const char*)&dummy, 4 ) ;
			dna_.write( (const char*)&dummy, 4 ) ;
			cur_genome_.clear_sequence() ;
			cur_genome_.set_maxn( maxn ) ;
		}

		void finish() 
		{
			step( 0 ) ;
			if( verbose_ ) std::clog << "\33[K" << std::flush ;

			uint32_t ix_start = dna_.tellp() ;
			if( !cur_genome_.SerializeToOstream( &dna_ ) ) 
				throw "couldn't serialize config" ;
			uint32_t ix_len = (uint32_t)dna_.tellp() - ix_start ;
			dna_.seekp( 4 ) ;
			dna_.write( (char*)&ix_start, 4 ) ;
			dna_.write( (char*)&ix_len, 4 ) ; 
		}

		void consume( std::istream& s ) 
		{
			int c ;
			while( (c = s.get()) != std::istream::traits_type::eof() ) step( c ) ;
		}

		void step( char c ) { (this->*state_)( c ) ; }

		void report() const {
			if( verbose_ )
				std::clog << cur_sequence_->name() << " @ " << source_position_
					      << "\33[K\r" << std::flush ;
		}

		void begin_new_sequence()
		{
			source_position_ = 0 ;
			cur_sequence_ = cur_genome_.add_sequence() ;
		}
		void finish_sequence() { cur_sequence_ = 0 ; }
		void begin_new_contig() 
		{ 
			put_nt( '-' ) ;
			num_n_ = 0 ;
			cur_contig_ = cur_sequence_->add_contig() ;
			cur_contig_->set_offset( position_ ) ;
			cur_contig_->set_range_start( source_position_ ) ;
			report() ;
		}
		void finish_contig()
		{ 
			cur_contig_->set_range_end( source_position_ ) ;
			report() ;
		}
		void finish_genome()
		{
			put_nt( '-' ) ;
			cur_genome_.set_total_size( position_ ) ;
			if( position_ & 1 ) put_nt( '-' ) ;
		}

		void put_nt( char c )
		{
			if( c == 'n' || c == 'N' ) 
			{
				++num_n_ ;
			}
			else
			{
				for( ; num_n_ ; --num_n_ ) put_nt_( '-' ) ;
				put_nt_( c ) ;
			}
		}

		void put_nt_( char c )
		{
			if( position_ & 1 )
			{
				dna_.put( (to_ambicode(c) << 4) | one_nt_ ) ;
			}
			else
			{
				one_nt_ = to_ambicode(c) ;
			}
			++position_ ;
			if( !(position_ & 0xfffff ) ) report() ;
		}

		void s_idle( char c )
		{
			if( c == '>' ) {
				begin_new_sequence() ;
				state_ = &FastaDecoder::s_name ;
			}
		}

		void s_name( char c )
		{
			if( !c ) finish_genome() ;
			else if( std::isspace(c) || c == '\n' )
			{
				if( c == '\n' ) state_ = &FastaDecoder::s_seq_begin ;
				else state_ = &FastaDecoder::s_info ;
			}
			else cur_sequence_->mutable_name()->push_back( c ) ;
		}

		void s_info( char c )
		{
			if( !c ) finish_genome() ; 
			else if( c == '\n' ) state_ = &FastaDecoder::s_seq_begin ;
		}

		void s_seq_begin( char c )
		{
			if( !c ) finish_genome() ;
			else if( c == ';' ) state_ = &FastaDecoder::s_info ;
			else if( c == '>' ) {
				begin_new_sequence() ;
				state_ = &FastaDecoder::s_name ;
			}
			else if( c == 'n' || c == 'N' ) {
				if( num_n_ == max_n_ ) {
					++source_position_ ;
					num_n_ = 0 ;
					state_ = &FastaDecoder::s_contig_gap ;
				}
				else
				{
					++source_position_ ;
					put_nt( c ) ;
				}
			}
			else if( encodes_nuc( c ) ) {
				begin_new_contig() ;
				++source_position_ ;
				put_nt( c ) ;
				state_ = &FastaDecoder::s_sequence ;
			}
		}

		void s_sequence( char c )
		{
			if( !c ) {
				source_position_ -= num_n_ ;
				finish_contig() ;
				finish_sequence() ;
				finish_genome() ;
			}
			else if( (c == 'N' || c == 'n') && num_n_ == max_n_ ) {
				source_position_ -= num_n_ ;
				finish_contig() ;
				source_position_ += num_n_+1 ;
				num_n_ = 0 ;
				state_ = &FastaDecoder::s_contig_gap ;
			}
			else if( encodes_nuc( c ) ) {
				++source_position_ ;
				put_nt( c ) ;
			}
			else if( c == '>' ) {
				source_position_ -= num_n_ ;
				finish_contig() ;
				finish_sequence() ;
				begin_new_sequence() ;
				state_ = &FastaDecoder::s_name ;
			}
		}

		void s_contig_gap( char c ) 
		{
			if( !c ) {
				finish_sequence() ;
				finish_genome() ;
			}
			else if( c == '>' ) {
				finish_sequence() ;
				begin_new_sequence() ;
				state_ = &FastaDecoder::s_name ;
			}
			else if( encodes_nuc( c ) && c != 'n' && c != 'N' ) {
				begin_new_contig() ;
				++source_position_ ;
				put_nt( c ) ;
				state_ = &FastaDecoder::s_sequence ;
			}
			else if( c == 'n' || c == 'N' ) {
				++source_position_ ;
			}
		}
} ;

std::string drop_suffix( const std::string& suf, const std::string& s )
{
	if( s.substr( s.size()-suf.size() ) == suf ) return s.substr( 0, s.size()-suf.size() ) ;
	else return s ;
}

int main_( int argc, const char * argv[] )
{
	GOOGLE_PROTOBUF_VERIFY_VERSION ;
	enum option_tags { opt_none, opt_version } ;

	const char* output_file = 0 ;
	const char* description = 0 ;
	const char* genome_name = 0 ;
	int max_num_n = 2 ;
	int verbose = 0 ;

	struct poptOption options[] = {
		{ "version",     'V', POPT_ARG_NONE,   0,            opt_version, "Print version number and exit", 0 },
		{ "output",      'o', POPT_ARG_STRING, &output_file, opt_none,    "Write DNA output to FILE", "FILE" },
		{ "maxn",        'm', POPT_ARG_INT,    &max_num_n,   opt_none,    "Treat N consecutive Ns as separator", "N" },
		{ "genome",      'g', POPT_ARG_STRING, &genome_name, opt_none,    "Set genome name to NAME", "NAME" },
		{ "description", 'd', POPT_ARG_STRING, &description, opt_none,    "Add TEXT as description to genome", "TEXT" },
		{ "verbose",     'v', POPT_ARG_NONE,   &verbose,     opt_none,    "Make more noise while working", 0 },
		POPT_AUTOHELP POPT_TABLEEND
	} ;

	poptContext pc = poptGetContext( "fa2dna", argc, argv, options, 0 ) ;
	poptSetOtherOptionHelp( pc, "[OPTION...] [FASTA-file...]" ) ;
	poptReadDefaultConfig( pc, 0 ) ;

	if( argc <= 1 ) { poptPrintHelp( pc, stderr, 0 ) ; return 1 ; }
	for( int rc = poptGetNextOpt( pc ) ; rc > 0 ; rc = poptGetNextOpt(pc) ) switch( rc )
	{
		case opt_version:
			std::cout << poptGetInvocationName(pc) << ", revision " << PACKAGE_VERSION << std::endl ;
			return 0 ;

		default:
			std::cerr << argv[0] << ": " << poptStrerror( rc ) 
				<< ' ' << poptBadOption( pc, 0 ) << std::endl ;
			return 1 ; 
	}

	if( !genome_name ) throw "missing --genome option" ;
	config::Genome g ;

	if( description ) g.set_description( description ) ;
	if( genome_name ) g.set_name( genome_name ) ;
	else if( !poptPeekArg( pc ) ) g.set_name( "genome" ) ;
	else g.set_name( drop_suffix( ".fa", drop_suffix( ".fas", drop_suffix( ".fna", poptPeekArg( pc ))))) ;

	std::ofstream output_stream( 
			(output_file ? output_file : (g.name() + ".dna")).c_str(),
			std::ios::trunc ) ;
	FastaDecoder fd( output_stream, g, max_num_n, verbose ) ;

	if( !poptPeekArg( pc ) ) fd.consume( std::cin ) ;
	while( const char* arg = poptGetArg( pc ) )
	{
		if( strcmp(arg,"-") ) 
		{
			std::ifstream str( arg ) ;
			fd.consume( str ) ;
		}
		else fd.consume( std::cin ) ;
	}

	fd.finish() ;
	return 0 ;
}


