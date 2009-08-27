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

//! \brief writes in google's text format
//! This is the human readable version of the native format.
//! Additionally, alignment strings and a CLUSTAL-style 'conservation'
//! line are added, provided the reference genome is available.
class TextWriter : public Stream
{
	private:
		google::protobuf::io::FileOutputStream fos_ ;

		void print_msg( const google::protobuf::Message& ) ;

	public:
		TextWriter( int fd ) : fos_( fd ) {}
		TextWriter( const char* fn ) 
			: fos_( throw_errno_if_minus1( creat( fn, 0666 ), "creating", fn ) )
		{ fos_.SetCloseOnDelete( true ) ; }
		virtual ~TextWriter() {}

		virtual void put_header( const Header& h ) 
		{
			Stream::put_header( h ) ;
			print_msg( h ) ;
		}

		virtual void put_result( const Result& ) ;

		virtual void put_footer( const Footer& f )
		{
			Stream::put_footer( f ) ;
			print_msg( f ) ;
		}
} ;

//! \brief writes in SAM format
class SamWriter : public Stream
{
	private:
		enum bad_stuff { goodness = 0, no_hit, multiple_hits, no_seqid, no_seq, bad_cigar, bad_stuff_max } ; 
		static const char *descr[] ;

		std::auto_ptr< std::filebuf > buf_ ;
		std::ostream out_ ;
		const char *g_ ;
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
		SamWriter( const char* fn, const char* g ) : buf_( new std::filebuf ), out_( buf_.get() ), g_(g)
		{
			buf_->open( fn, std::ios_base::binary | std::ios_base::out | std::ios_base::trunc ) ;
			memset( discarded, 0, sizeof(discarded) ) ;
		}

		SamWriter( std::streambuf *s, const char* g ) : buf_(), out_( s ), g_(g)
		{
			memset( discarded, 0, sizeof(discarded) ) ;
		}

		virtual ~SamWriter() {}

		virtual void put_header( const Header& h ) 
		{
			Stream::put_header( h ) ;
			out_ << "@HD\tVN:1.0" ;
			if( h.has_is_sorted_by_coordinate() && (
						(!g_ && h.is_sorted_by_coordinate().empty())
						|| h.is_sorted_by_coordinate() == g_ ) )
				out_ << "\tSO:coordinate" ;
			else if( h.is_sorted_by_name() ) out_ << "\tSO:queryname" ;
			out_ << "\n@PG\tID:ANFO\tVN:" << h.version() << '\n' ;
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
		std::auto_ptr< std::filebuf > buf_ ;
		std::ostream out_ ;
		const char* g_ ;
		int c_ ;

	public:
		FastaAlnWriter( const char* fn, const char* g, int c ) : buf_( new std::filebuf ), out_( buf_.get() ), g_(g), c_(c)
		{
			buf_->open( fn, std::ios_base::binary | std::ios_base::out | std::ios_base::trunc ) ;
		}

		FastaAlnWriter( std::streambuf *s, const char* g, int c ) : buf_(), out_( s ), g_(g), c_(c) { }
		virtual ~FastaAlnWriter() {}

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
		std::auto_ptr< std::filebuf > buf_ ;
		std::ostream out_ ;

	public:
		FastqWriter( const char* fn ) : buf_( new std::filebuf ), out_( buf_.get() )
		{
			buf_->open( fn, std::ios_base::binary | std::ios_base::out | std::ios_base::trunc ) ;
		}

		FastqWriter( std::streambuf *s ) : buf_(), out_( s ) {}
		virtual ~FastqWriter() {}

		virtual void put_result( const Result& ) ;
} ;

class TableWriter : public Stream
{
	private:
		std::auto_ptr< std::filebuf > buf_ ;
		std::ostream out_ ;
		const char *g_ ;

	public:
		TableWriter( const char* fn, const char* g ) : buf_( new std::filebuf ), out_( buf_.get() ), g_(g)
		{
			buf_->open( fn, std::ios_base::binary | std::ios_base::out | std::ios_base::trunc ) ;
		}

		TableWriter( std::streambuf *s, const char* g ) : buf_(), out_( s ), g_(g) { }
		virtual ~TableWriter() {}

		virtual void put_result( const Result& ) ;
} ;

} // namespace

#endif

