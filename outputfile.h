#ifndef INCLUDED_OUTPUTFILE_H
#define INCLUDED_OUTPUTFILE_H

#include "output.pb.h"
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <fstream>
#include <memory>

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

//! \brief presents ANFO files as series of messages
class AnfoFile
{
	private:
		google::protobuf::io::FileInputStream iis_ ;
		std::auto_ptr<google::protobuf::io::ZeroCopyInputStream> zis_ ;
		bool error_ ;
		std::string name_ ;
		output::Footer foot_ ;

		void check_valid_file() ;

	public: 
		//! \brief opens the named file
		AnfoFile( const std::string& name ) ;

		//! \brief uses a given name and filedescriptor
		//! No actual file is touched, the name for informational
		//! purposes only.
		AnfoFile( int fd, const std::string& name ) ;

		//! \brief reads the header messages
		//! This must be called first after constructing the object.  If
		//! the first message is not a header, an error results.
		output::Header read_header() ;

		//! \brief reads a result message
		//! This should be called repeatedly after having read the
		//! header.  If no more messages follow, an empty Result is
		//! returned (and the footer is stored).  Valid Results from
		//! ANFO files will always have a seq_id.
		output::Result read_result() ;

		//! \brief reads the footer message
		//! This can only be called after all Result messages have been
		//! consumed.  If an error occured during parsing, the exit_code
		//! inside the footer message will have its LSB set.
		output::Footer read_footer() { if( error_ ) foot_.set_exit_code( 1 | foot_.exit_code() ) ; return foot_ ; }
} ;

void merge_sensibly( output::Header& lhs, const output::Header& rhs ) ;
void merge_sensibly( output::Footer& lhs, const output::Footer& rhs ) ;
void merge_sensibly( output::Result& lhs, const output::Result& rhs ) ;

#endif
