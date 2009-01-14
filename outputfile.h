#ifndef INCLUDED_OUTPUTFILE_H
#define INCLUDED_OUTPUTFILE_H

#include "output.pb.h"
#include <google/protobuf/io/coded_stream.h>

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
//! @{

template< typename Msg >
void write_delimited_message( google::protobuf::io::CodedOutputStream& os, const Msg& m )
{
	std::string code ;
	if( !m.SerializeToString( &code ) ||
			!os.WriteVarint32( code.size() ) ||
			!os.WriteString( code ) )
		throw "error while serializing" ;
}

template< typename Msg >
bool read_delimited_message( google::protobuf::io::CodedInputStream& is, Msg &m )
{
	uint32_t size ;
	std::string code ;
	if( !is.ReadVarint32( &size ) ) return false ;
	if( !is.ReadString( &code, size ) || !m.ParseFromString( code ) )
		throw "error while deserializing" ;
	return true ;
}

template< typename Hdr, typename Fun >
void reduce_output_file( google::protobuf::io::CodedInputStream& is, Hdr h, Fun f )
{
	output::Header hdr ;
	output::Result res ;
	if( !read_delimited_message( is, hdr ) ) throw "cannot read header" ;
	h( hdr ) ;
	while( read_delimited_message( is, res ) ) f( hdr, res ) ;
}

//! @}

#endif
