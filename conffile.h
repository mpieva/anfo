//    Copyright 2009 Udo Stenzel
//    This file is part of ANFO
//
//    ANFO is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    Foobar is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

#ifndef INCLUDED_CONFFILE_H
#define INCLUDED_CONFFILE_H

#include "config.pb.h"
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/text_format.h>
#include <fstream>
#include <string>

/*! \page configuration Configuration Files, Metadata, Output Format
 *
 * In order to keep configuration files and somewhat complicated,
 * persistent data structures in a sensible format and extensible, we
 * abuse Google's protocol buffers (http://code.google.com/p/protobuf/)
 * for this.  Code to handle protocol buffers can be generated for C++,
 * Java and Pyhton by the Google toolchain, but also for Haskell
 * (http://hackage.haskell.org/cgi-bin/hackage-scripts/package/hprotoc)
 * and possibly for more languages.  The definition of the possible
 * parts of the configuration resides in \c config.proto, we store the
 * configuration in a single file using the protobuf text format.
 * Metadata for genomes and indices is also defined there, those files
 * incorporate their metadata in binary form.  If the need arises, we
 * can easily move to multiple files, binary files, or a combination of
 * both.  Text and binary files are interconvertible using \c protoc.
 */

//! \brief reads the configuration from a text file.
//! A text format configuration file is read and parsed into a config
//! object.
//!
//! \param filename filename of configuration file
//! \return parsed configuration 
inline config::Config parse_text_config( const std::string& filename )
{
	config::Config c ;
	std::ifstream config_in( filename.c_str() ) ;
	google::protobuf::io::IstreamInputStream fis( &config_in ) ;
	if( !google::protobuf::TextFormat::Parse( &fis, &c ) )
		throw "parse error reading " + filename ;
	return c ;
}

inline config::Config get_default_config( const char* config_file = 0 ) 
{
	if( config_file ) return parse_text_config( config_file ) ;
	else if( !access( "anfo.cfg", F_OK ) ) return parse_text_config( "anfo.cfg" ) ;
	else if( !access( ".anfo.cfg", F_OK ) ) return parse_text_config( ".anfo.cfg" ) ;
	else {
		std::string f = getenv( "HOME" ) + std::string( ".anfo.cfg" ) ;
		if( !access( f.c_str(), F_OK ) ) return parse_text_config( f.c_str() ) ;
		else throw "no config file found" ;
	}
}

#endif

