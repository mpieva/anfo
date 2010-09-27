//    Copyright 2009 Udo Stenzel
//    This file is part of ANFO
//
//    ANFO is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    Anfo is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with Anfo.  If not, see <http://www.gnu.org/licenses/>.

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include "align_fwd.h"
#include "conffile.h"
#include "config.pb.h"
#include "index.h"
#include "stream.h"
#include "util.h"

#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/text_format.h>

#include <fstream>
#include <iostream>

WRAPPED_MAIN
{
	GOOGLE_PROTOBUF_VERIFY_VERSION ;
	int ret = 0 ;
	bool quiet = false ;
	const char** arg = argv+1 ;
	if( arg != argv+argc && strcmp(*arg,"-q") == 0 ) { quiet = true ; ++arg ; }

	for( ; arg != argv+argc ; ++arg )
	{
		if( !quiet ) std::cout << *arg << ": " ;
		std::ifstream fd( *arg ) ;
		uint32_t sig ;
		fd.read( (char*)&sig, 4 ) ;
		if( sig == CompactGenome::signature )
		{
			if( !quiet ) std::cout << "DNA file" << std::endl ;
			uint32_t off = 0, len = 0 ;
			fd.read( (char*)&off, 4 ) ;
			fd.read( (char*)&len, 4 ) ;
			char buf[ len ] ;
			fd.seekg( off ) ;
			fd.read( (char*)buf, len ) ;
			config::Genome g ;
			g.ParseFromArray( buf, len ) ;
			if( !quiet ) {
				google::protobuf::io::OstreamOutputStream os( &std::cout ) ;
				google::protobuf::TextFormat::Print( g, &os ) ;
			}
		}
		else if( sig == FixedIndex::signature || (sig & 0xffffff) == (FixedIndex::signature & 0xffffff) )
		{
			if( !quiet ) std::cout << "inverted list index"
				<< ( sig == FixedIndex::signature ? "" : " (legacy format)" ) 
				<< std::endl ;
			uint32_t len = 0 ;
			fd.read( (char*)&len, 4 ) ;
			char buf[ len ] ;
			fd.read( (char*)buf, len ) ;
			config::CompactIndex ci ;
			ci.ParseFromArray( buf, len ) ;
			if( !quiet ) {
				google::protobuf::io::OstreamOutputStream os( &std::cout ) ;
				google::protobuf::TextFormat::Print( ci, &os ) ;
			}
		}
		else 
		{
			output::Header header ;
			Holder< streams::Stream > af( new streams::UniversalReader( *arg ) ) ;

			try { header = af->fetch_header() ; } catch(...) {}
			if( header.has_config() ) 
			{
				if( !quiet ) std::cout << "ANFO stream file" << std::endl ;
				if( !quiet ) std::cout << adna_parblock( header.config().aligner() ) ;

				try {
					output::Footer footer = af->lookahead_footer() ;
					if( !quiet ) {
						google::protobuf::io::OstreamOutputStream os( &std::cout ) ;
						google::protobuf::TextFormat::Print( header, &os ) ;
						google::protobuf::TextFormat::Print( footer, &os ) ;
					}
					ret |= footer.IsInitialized() ? footer.exit_code() : 255 ;
				} 
				catch(...) { ret |= 255 ; }
			}
			else
			{
				try {
					config::Config conf = parse_text_config( *arg ) ;
					if( !quiet ) {
						std::cout << "ANFO configuration file" << std::endl ;
						google::protobuf::io::OstreamOutputStream os( &std::cout ) ;
						google::protobuf::TextFormat::Print( conf, &os ) ;
						std::cout << adna_parblock( conf.aligner() ) ;
					}
				}
				catch(...)
				{
					if( !quiet ) {
						if( af->get_state() == streams::Stream::have_output )
							std::cout << "FastA/FastQ sequence file(?)" << std::endl ;
						else std::cout << "unknown\n" ;
					}
				}
			}
		}
		if( !quiet ) std::cout << std::endl ;
	}
	return ret ;
}
