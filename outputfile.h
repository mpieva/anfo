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

//! \brief stream of result messages
//! Each stream is a header, followed by many results, followed by a
//! single footer.  The header will be cached internally, so it can be
//! asked for repeatedly.  Results are forgotten once read, the footer
//! is only available after all results have been read, and then it is
//! stored and can be read repeatedly.  It is undefined behaviour to
//! request the footer without first consuming all results.
class Stream
{
	public:
		virtual ~Stream() {}

		//! \brief returns the header
		//! The header can be requested repeatedly, once the stream was
		//! correctly intialized.
		virtual const output::Header& get_header() = 0 ;

		//! \brief returns the next result
		//! Every result can only be read once, internal iterator style.
		//! \param r place to store the next result in
		//! \return true iff another result was available
		virtual bool read_result( output::Result& r ) = 0 ;

		//! \brief returns the footer
		//! Only after all results have been consumed is the footer
		//! available.
		virtual const output::Footer& get_footer() = 0 ;
} ;

//! \brief presents ANFO files as series of messages
//! This class will read a possibly compressed result stream and present
//! it using an iteration interface.  Normally, gzip and bzip2
//! decompression will transparently be done.
class AnfoFile : public Stream
{
	private:
		google::protobuf::io::FileInputStream iis_ ;
		std::auto_ptr<google::protobuf::io::ZeroCopyInputStream> zis_ ;
		std::string name_ ;

		output::Header hdr_ ;
		output::Footer foot_ ;

		void initialize() ;
		static int num_files_ ; // tracked to avoid bumping into the file descriptor limit

	public: 
		//! \brief opens the named file
		AnfoFile( const std::string& name ) ;
		virtual ~AnfoFile() { --num_files_ ; }

		//! \brief uses a given name and filedescriptor
		//! No actual file is touched, the name is for informational
		//! purposes only.
		AnfoFile( int fd, const std::string& name ) ;

		virtual const output::Header& get_header() { return hdr_ ; }
		virtual bool read_result( output::Result& ) ;

		//! \brief reads the footer message
		//! If an error occured during earlier processing, the exit_code
		//! inside the footer message will have its LSB set.
		virtual const output::Footer& get_footer() { return foot_ ; }

		//! \internal
		static unsigned num_open_files() { return num_files_ ; }
} ;

void merge_sensibly( output::Header& lhs, const output::Header& rhs ) ;
void merge_sensibly( output::Footer& lhs, const output::Footer& rhs ) ;
void merge_sensibly( output::Result& lhs, const output::Result& rhs ) ;

//! \brief writes a stream of results to a file
//! The file will be in a format that can be read in by ::AnfoFile.
//! \param fd file descriptor to write to
//! \param s stream to copy to the file
//! \param expensive if set, computationally expensive compression is
//!                  used (normally meaning bzip2 instead of gzip)
int write_stream_to_file( int fd, Stream &s, bool expensive = false ) ;

//! \brief writes a stream of results to an output stream
//! The byte stream will be in a format that can be read in by ::AnfoFile.
//! \param zos Google-style output stream to write to
//! \param s stream to copy to the file
int write_stream_to_file( google::protobuf::io::ZeroCopyOutputStream *zos, Stream &s ) ;

#endif
