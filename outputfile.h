#ifndef INCLUDED_OUTPUTFILE_H
#define INCLUDED_OUTPUTFILE_H

#include "output.pb.h"
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <fstream>

//! \defgroup outputfile Convenience functions to handle output file
//! The output file is a series of protocol buffer messages.  This makes
//! it compact and extensible.  Instead of writing a serial file, the
//! very same messages could be entered into a Berkeley-DB type
//! database.  Since conversion to all sorts of text files is a no
//! brainer, the binary format is not considered undue hardship on the
//! users of scripting languages.
//!
//! (If you are a user of a glorified logfile munching language and want
//! an ASCII-based table, define a sensible layout and ask for it.  Or
//! write the single function yourself, it isn't that hard.)
//!
//! The protobuf format itself is not self-delimiting, therefore we
//! prepend each message with the number of bytes it occupies, in
//! variable protobuf format.  The first message in a file
//! one is of type output::Header and should contain all sorts of
//! meta information, maybe even configuration of the mapper that was
//! used to create it (an option would be to include the full set of
//! policies).  Subsequent messages are of type output::Result.
//!
//! \todo The output is delimited, but not tagged.  We should do both,
//!       in essence incrementally writing a protobuf message.
//!
//! @{

template< typename Msg >
void write_delimited_message( google::protobuf::io::CodedOutputStream& os, int tag, const Msg& m )
{
	if(
			!os.WriteTag( (tag << 3) | 2 ) ||
			!os.WriteVarint32( m.ByteSize() ) ||
			!m.SerializeToCodedStream( &os )
	  )
		throw "error while serializing" ;
}

template< typename Msg >
bool read_delimited_message( google::protobuf::io::CodedInputStream& is, Msg &m )
{
	uint32_t size ;
	std::string code ;
	if( !is.ReadVarint32( &size ) ) return false ;
	int lim = is.PushLimit( size ) ;
	if( !m.ParseFromCodedStream( &is ) ) throw "error while deserializing" ;
	is.PopLimit( lim ) ;
	return true ;
}

template< typename Hdr, typename Fun, typename Foot >
void scan_output_file( google::protobuf::io::CodedInputStream& is, Hdr h, Fun f, Foot t )
{
	std::string buf ;
	if( !is.ReadString( &buf, 4 ) || buf != "ANFO" ) throw "not an ANFO output file" ;
	while( uint32_t tag = is.ReadTag() )
	{
		if( (tag>>3) == 1 ) {
			output::Header hdr ;
			read_delimited_message( is, hdr ) ;
			h( hdr ) ;
		}
		else if( (tag>>3) == 2 ) {
			output::Result res ;
			read_delimited_message( is, res ) ;
			f( res ) ;
		}
		else if( (tag>>3) == 3 ) {
			output::Footer foot ;
			read_delimited_message( is, foot ) ;
			t( foot ) ;
		}
	}
}

//! @}

class AnfoFile
{
    private:
	const char *name_ ;
	std::ifstream ifs_ ;
	google::protobuf::io::IstreamInputStream iis_ ;
	google::protobuf::io::CodedInputStream cis_ ;
	bool legacy_ ;
	bool error_ ;
	output::Footer foot_ ;

    public: 
	AnfoFile( const char *name ) ;

	output::Header read_header() ;
	output::Result read_result() ;
	output::Footer read_footer() { if( error_ ) foot_.set_exit_code( 1 | foot_.exit_code() ) ; return foot_ ; }
} ;


#endif
