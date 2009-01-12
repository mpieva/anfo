#include "config.pb.h"
#include "index.h"
#include "util.h"

#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/text_format.h>

#include <fstream>
#include <iostream>

int main_( int argc, const char**argv )
{
	for( int argi = 1 ; argi != argc ; ++argi )
	{
		std::cout << argv[argi] << ": " ;
		std::ifstream fd( argv[argi] ) ;
		uint32_t sig ;
		fd.read( (char*)&sig, 4 ) ;
		if( sig == CompactGenome::signature )
		{
			std::cout << "DNA file" << std::endl ;
			uint32_t off = 0, len = 0 ;
			fd.read( (char*)&off, 4 ) ;
			fd.read( (char*)&len, 4 ) ;
			char buf[ len ] ;
			fd.seekg( off ) ;
			fd.read( (char*)buf, len ) ;
			config::Genome g ;
			g.ParseFromArray( buf, len ) ;
			google::protobuf::io::OstreamOutputStream os( &std::cout ) ;
			google::protobuf::TextFormat::Print( g, &os ) ;
		}
		else if( sig == FixedIndex::signature )
		{
			std::cout << "inverted list index" << std::endl ;
			uint32_t len = 0 ;
			fd.read( (char*)&len, 4 ) ;
			char buf[ len ] ;
			fd.read( (char*)buf, len ) ;
			config::CompactIndex ci ;
			ci.ParseFromArray( buf, len ) ;
			google::protobuf::io::OstreamOutputStream os( &std::cout ) ;
			google::protobuf::TextFormat::Print( ci, &os ) ;
		}
		else std::cout << "unknown\n" ;
		std::cout << std::endl ;
	}
	return 0 ;
}
