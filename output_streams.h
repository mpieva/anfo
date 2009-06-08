#ifndef INCLUDED_OUTPUT_STREAMS_H
#define INCLUDED_OUTPUT_STREAMS_H

#include "index.h"
#include "stream.h"

#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <fstream>
#include <ios>
#include <map>

namespace output { class Hit ; }
namespace config { class Config ; }

namespace streams
{

class TextWriter : public Stream
{
	private:
		google::protobuf::io::FileOutputStream fos_ ;

		void print_msg( const google::protobuf::Message& ) ;
		void add_alignment( std::string::const_iterator qry, output::Hit &h, const config::Config &conf ) ;
		
		typedef std::map< std::string, CompactGenome > Genomes ; 
		Genomes genomes_ ;

	public:
		TextWriter( int fd ) : fos_( fd ) {}
		TextWriter( const char* fn ) 
			: fos_( throw_errno_if_minus1( creat( fn, 0777 ), "creating", fn ) )
		{ fos_.SetCloseOnDelete( true ) ; }
		virtual ~TextWriter() {}

		virtual void put_header( const Header& h ) 
		{
			hdr_ = h ;
			print_msg( h ) ;
			state_ = need_input ;
		}

		virtual void put_result( const Result& r )
		{
			if( r.has_best_to_genome() && r.best_to_genome().has_cigar() )
			{
				Result r_ = r ;
				add_alignment( r_.sequence().begin(), *r_.mutable_best_to_genome(), hdr_.config() ) ;
				print_msg( r_ ) ;
			}
			else print_msg( r ) ;
		}

		virtual void put_footer( const Footer& f )
		{
			print_msg( f ) ;
			state_ = end_of_stream ;
		}
} ;

class SamWriter : public Stream
{
	private:
		enum bad_stuff { goodness = 0, no_hit, multiple_hits, no_seqid, no_seq, bad_cigar, bad_stuff_max } ; 
		static const char *descr[] ;

		std::auto_ptr< std::filebuf > buf_ ;
		std::ostream out_ ;
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
		SamWriter( const char* fn ) : buf_( new std::filebuf ), out_( buf_.get() )
		{
			buf_->open( fn, std::ios_base::binary | std::ios_base::out | std::ios_base::trunc ) ;
			memset( discarded, 0, sizeof(discarded) ) ;
		}

		SamWriter( std::streambuf *s ) : buf_(), out_( s ) 
		{
			memset( discarded, 0, sizeof(discarded) ) ;
		}

		virtual ~SamWriter() {}

		virtual void put_header( const Header& h ) 
		{
			out_ << "@HD\tVN:1.0" ;
			if( h.is_sorted_by_coordinate() ) out_ << "\tSO:coordinate" ;
			else if( h.is_sorted_by_name() ) out_ << "\tSO:queryname" ;
			out_ << "\n@PG\tID:ANFO\tVN:" << h.version() << '\n' ;
			state_ = need_input ;
		}

		virtual void put_result( const Result& res )
		{
			if (bad_stuff r = protoHit_2_bam_Hit( res )) discarded[r]++;
		}

		virtual void put_footer( const Footer& ) ;
} ;
} // namespace

#endif

