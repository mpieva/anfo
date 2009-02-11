#include "output.pb.h"
#include "outputfile.h"

#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/text_format.h>

using namespace google::protobuf::io ;
using namespace google::protobuf ;
using namespace config ;
using namespace output ;

FileOutputStream cout(1) ;

template< typename Msg > void print_msg( const Msg& m )
{
	TextFormat::Print( m, &cout ) ;
	void* data ; int size ;
	cout.Next( &data, &size ) ;
	*(char*)data = '\n' ;
	data = (char*)data + 1 ;
	--size ;
	if( !size ) cout.Next( &data, &size ) ;
	*(char*)data = '\n' ;
	data = (char*)data + 1 ;
	--size ;
	if( size ) cout.BackUp( size ) ;
}

int main_( int argc, const char** argv )
{
	GOOGLE_PROTOBUF_VERIFY_VERSION ;
	for( int argi = 1 ; argi != argc ; ++argi )
	{
		AnfoFile af( argv[argi] ) ;
		print_msg( af.read_header() ) ;
		for(;;) {
			Result r = af.read_result() ;
			if( !r.has_seqid() ) break ;
			print_msg( r ) ;
		}
		print_msg( af.read_footer() ) ;
	}
	return 0 ;
}

