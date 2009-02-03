#include "outputfile.h"
#include <iostream>

using namespace google::protobuf::io ;
using namespace output ;
using namespace std ;

AnfoFile::AnfoFile( const char *name ) 
	: name_( name ), ifs_( name ), iis_( &ifs_ ), cis_( &iis_ ), legacy_(false), error_(false)
{
	std::string tag ;
	cis_.SetTotalBytesLimit( INT_MAX, INT_MAX ) ;
	if( !cis_.ReadString( &tag, 4 ) || tag != "ANFO" ) {
		clog << "\033[K" << name_ << ": not an ANFO file" << endl ;
		error_ = true ;
	}
	foot_.set_exit_code(0) ;
}

Header AnfoFile::read_header()
{
	uint32_t tag ;
	Header hdr ;
	if( cis_.ReadVarint32( &tag ) )
	{
		if( tag != 10 )
		{
			legacy_ = true ;
			clog << "\033[K" << name_ << ": legacy file (" << tag << ")" << endl ;
		}
		if( tag != 10 || cis_.ReadVarint32( &tag ) )
		{
			int lim = cis_.PushLimit( tag ) ;
			if( hdr.ParseFromCodedStream( &cis_ ) ) {
				cis_.PopLimit( lim ) ;
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
	if( !error_ ) {
		if( cis_.ExpectAtEnd() ) {
			clog << "\033[K" << name_ << ": end of stream" << endl ;
			return Result() ;
		}
		if( legacy_ || (tag = cis_.ReadTag()) ) {
			uint32_t size ;
			std::string buf ;
			if( cis_.ReadVarint32( &size ) && cis_.ReadString( &buf, size ) ) {
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
