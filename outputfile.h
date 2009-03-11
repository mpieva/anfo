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
//! \note The output file is essentially a single, giant protobuf message.
//!       However, it is impossible to read it in one go (due to memory
//!       constraints), but it's also impossible to read it through a single
//!       CodedInputStream.  That's why methods in here tend to construct and
//!       destruct a fresh CodedInputStream on each invocation.
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

//! @}

class AnfoFile
{
    private:
	std::string name_ ;
	int fd_ ;
	google::protobuf::io::FileInputStream iis_ ;
	bool legacy_ ;
	bool error_ ;
	bool unlink_on_delete_ ;
	output::Footer foot_ ;

    public: 
	AnfoFile( const std::string& name, bool unlink_on_delete = false ) ;
	~AnfoFile() ;

	output::Header read_header() ;
	output::Result read_result() ;
	output::Footer read_footer() { if( error_ ) foot_.set_exit_code( 1 | foot_.exit_code() ) ; return foot_ ; }
} ;


#endif
