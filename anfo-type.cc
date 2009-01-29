#include "output.pb.h"
#include "outputfile.h"

#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/text_format.h>
#include <fstream>
#include <iostream>

using namespace google::protobuf::io ;
using namespace google::protobuf ;
using namespace config ;
using namespace output ;

template< typename Msg > void print_msg( const Msg& m )
{
	{
		OstreamOutputStream oos( &std::cout ) ;
		TextFormat::Print( m, &oos ) ;
	}
	std::cout << '\n' << std::endl ;
}

void print_hdr( const Header& h ) { print_msg( h ) ; }
void print_foot( const Footer& h ) { print_msg( h ) ; }
void print_res( const Header&, const Result& r ) { print_msg( r ) ; }

int main_( int argc, const char** argv )
{
	GOOGLE_PROTOBUF_VERIFY_VERSION ;
	for( int argi = 1 ; argi != argc ; ++argi )
	{
		std::ifstream is( argv[argi] ) ;
		IstreamInputStream iis( &is ) ;
		CodedInputStream cis( &iis ) ;
		reduce_output_file( cis, &print_hdr, &print_foot, &print_res ) ;
	}
	return 0 ;
}

