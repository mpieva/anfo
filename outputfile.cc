#include "outputfile.h"
#include "compress_stream.h"
#include "util.h"

#include <iostream>
#include <fcntl.h>

using namespace google::protobuf::io ;
using namespace output ;
using namespace std ;

AnfoFile::AnfoFile( const std::string& name, bool unlink_on_delete ) 
	: name_( name ), fd_( throw_errno_if_minus1( open( name.c_str(), O_RDONLY ), "opening ", name.c_str() ) )
        , iis_( fd_ ), zis_( decompress( &iis_ ) ), legacy_(false), error_(false), unlink_on_delete_(unlink_on_delete)
{
	std::string tag ;
	CodedInputStream cis( zis_.get() ) ;

	if( !cis.ReadString( &tag, 4 ) || tag != "ANFO" ) {
		clog << "\033[K" << name_ << ": not an ANFO file" << endl ;
		error_ = true ;
	}
	foot_.set_exit_code(0) ;
}

AnfoFile::~AnfoFile()
{
	close( fd_ ) ;
	if( unlink_on_delete_ ) unlink( name_.c_str() ) ;
}

Header AnfoFile::read_header()
{
	uint32_t tag ;
	Header hdr ;
	CodedInputStream cis( zis_.get() ) ;
	if( cis.ReadVarint32( &tag ) )
	{
		if( tag != 10 )
		{
			legacy_ = true ;
			clog << "\033[K" << name_ << ": legacy file (" << tag << ")" << endl ;
		}
		if( tag != 10 || cis.ReadVarint32( &tag ) )
		{
			int lim = cis.PushLimit( tag ) ;
			if( hdr.ParseFromCodedStream( &cis ) ) {
				cis.PopLimit( lim ) ;
				return hdr ;
			}
		}
	}

	clog << "\033[K" << name_ << ": deserialization error" << endl ;
	error_ = true ;
	return hdr ;
}

Result AnfoFile::read_result()
{
	uint32_t tag = 0 ;
	CodedInputStream cis( zis_.get() ) ;
	if( !error_ ) {
		if( cis.ExpectAtEnd() ) {
			clog << "\033[K" << name_ << ": end of stream" << endl ;
			return Result() ;
		}
		if( legacy_ || (tag = cis.ReadTag()) ) {
			uint32_t size ;
			std::string buf ;
			if( cis.ReadVarint32( &size ) && cis.ReadString( &buf, size ) ) {
				if( legacy_ || tag == 18 ) {
					Result res ;
					if( res.ParseFromString( buf ) && res.has_seqid() ) {
						return res ;
					}
				}
				if( legacy_ || tag == 26 ) {
					if( foot_.ParseFromString( buf ) ) return Result() ;
				}
				clog << "\033[K" << name_ << ": deserialization error" << endl ;
				error_ = true ;
			}
		}
	}
	return Result() ;
}
