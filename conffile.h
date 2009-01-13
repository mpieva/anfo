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

#endif

