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

#ifndef INCLUDED_OUTPUT_STREAMS_H
#define INCLUDED_OUTPUT_STREAMS_H

#include "index.h"
#include "stream.h"

#include <google/protobuf/io/gzip_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <fstream>
#include <ios>
#include <map>

namespace output { class Hit ; }
namespace config { class Config ; }

namespace streams {

	using namespace google::protobuf::io ;

//! \brief writes in (a modification of) Google's text format
//! This is essentially the human readable version of the native format.
//! Additionally, alignment strings and a CLUSTAL-style 'conservation'
//! line are added, provided they were looked up using the reference
//! genome beforehand, and some fields, notably the CIGAR line, are
//! printed in a more compact format.
class TextWriter : public Stream
{
	private:
		auto_ptr< std::ostream > out_ ;

		void print_msg( const google::protobuf::Message& ) ;

	public:
		TextWriter( const pair< std::ostream*, string > &p ) : out_( p.first ) {}

		virtual void put_header( const Header& ) ;
		virtual void put_result( const Result& ) ;
		virtual void put_footer( const Footer& ) ;
		virtual string type_name() const { return "TextWriter" ; }
} ;

//! \brief writes in SAM format
//! Every hit for every sequence is written.  A sequence without hits is
//! written out anyway.  If that's not desired, the input stream must be
//! filtered appropriately.
//! \todo Writing SAM is probably a bad idea in the long run.  Directly
//!       generating BAM is more sustainable
//! \todo Creation of correctly sorted SAM files is next to impossible.
//!       Maybe creation of BAM files is better anyway.
class SamWriter : public Stream
{
	private:
		enum bad_stuff { goodness = 0, no_hit, multiple_hits, no_seqid, no_seq, bad_cigar, bad_stuff_max } ; 
		static const char *descr[] ;

		std::auto_ptr< std::ostream > out_ ;
		string nm_ ;
		int discarded[bad_stuff_max] ;

		enum bam_flags {
			bam_fpaired        = 1,   // read is paired in sequencing
			bam_fproper_pair   = 2,   // read is mapped in proper pair
			bam_funmap         = 4,   // query seq. is unmapped
			bam_fmunmap        = 8,   // mate is umapped
			bam_freverse      = 16,   // strand of query (0 - fwd, 1 - rev)
			bam_fmreverse     = 32,   // strand of mate
			bam_fread1        = 64,   // read is 1. in pair
			bam_fread2       = 128,   // read is 2. in pair
			bam_fsecondary   = 256,   // alignment is NOT primary
			bam_fqcfail      = 512,   // read fails due low quality
			bam_fdup        = 1024    // read is duplicate
		} ;

		bad_stuff protoHit_2_bam_Hit( const output::Result& ) ;

	public:
		SamWriter( const pair< ostream*, string > &p ) : out_( p.first ), nm_( p.second )
		{
			memset( discarded, 0, sizeof(discarded) ) ;
		}

		virtual void put_header( const Header& h ) 
		{
			Stream::put_header( h ) ;
			*out_ << "@HD\tVN:1.0" ;
			if( h.is_sorted_by_name() ) *out_ << "\tSO:queryname" ;
			*out_ << "\n@PG\tID:ANFO\tVN:" << h.version() << '\n' ;
		}

		virtual void put_result( const Result& res )
		{
			if (bad_stuff r = protoHit_2_bam_Hit( res )) discarded[r]++;
		}

		virtual void put_footer( const Footer& ) ;
} ;


//! \brief writes in aligned FASTA-format
//! Each alignment appears as a pair of sequences, reference first,
//! query last.  Score, coordinates and whether an adapter was trimmed
//! are encoded in the header.  This is considered a legacy format,
//! mostly useful to make substitution graphs from.
class FastaAlnWriter : public Stream
{
	private:
		std::auto_ptr< std::ostream > out_ ;

	public:
		FastaAlnWriter( const pair< ostream*, string > &p ) : out_( p.first ) {}
		virtual void put_header( const Header& ) ;
		virtual void put_result( const Result& ) ;
} ;

//! \brief writes in FASTQ format
//! This should simply reproduce the input to ANFO, so I can throw away
//! the ugly FASTQ files.  Every sequence gets a header, then the
//! sequence (50 bases per line), the header is then *not* repeated and
//! the quality scores follow in the same layout.
class FastqWriter : public Stream
{
	private:
		std::auto_ptr< std::ostream > out_ ;
		bool with_qual_ ;

	public:
		FastqWriter( const pair< ostream*, string > &p, bool q ) : out_( p.first ), with_qual_(q) {}
		virtual void put_result( const Result& ) ;
} ;

class TableWriter : public Stream
{
	private:
		std::auto_ptr< std::ostream > out_ ;

	public:
		TableWriter( const pair< ostream*, string > &p ) : out_( p.first ) {}
		virtual void put_result( const Result& ) ;
} ;

class GenTextAlignment : public Filter
{
	private:
		int context_ ;

	public:
		GenTextAlignment( int context ) : context_( context ) {}
		virtual bool xform( Result& ) ;
} ;

//! \brief writes out coverage depth in WIG format
//! This only works after DuctTaper has been applied, afterwards it
//! extracts the depth of coverage (aka number of observations) per
//! position.
//! Note that the output file will get huge; we're talking about a text
//! based format...
class WigCoverageWriter : public Stream
{
	private:
		std::auto_ptr< std::ostream > out_ ;

	public:
		WigCoverageWriter( const pair< ostream*, string > &p ) : out_( p.first ) {}
		virtual void put_result( const Result& ) ;
} ;

} // namespace

#endif

