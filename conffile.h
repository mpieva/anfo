#ifndef INCLUDED_CONFFILE_H
#define INCLUDED_CONFFILE_H

#include "config.pb.h"

/*! \defgroup configuration Convenient access to configuration files
 *
 * Configuration files contain mostly metadata, e.g. which genomes are
 * known, what contigs do they consist of, how were they indexed.  It
 * might be useful to add other stuff, such as seeding policies and
 * alignment parameters.
 *
 * In order to keep configuration files in a sensible format and
 * extensible, we abuse Google's protocol buffers
 * (http://code.google.com/p/protobuf/) for this.  Code to handle
 * protocol buffers can be generated for C++, Java and Pyhton by the
 * Google toolchain, but also for Haskell
 * (http://hackage.haskell.org/cgi-bin/hackage-scripts/package/hprotoc)
 * and possibly for more languages.  The definition of the possible
 * parts of the configuration resides in \c config.proto, we
 * currently store the configuration in a single file using the protobuf
 * text format.  If the need arises, we can easily move to multiple
 * files, binary files, or a combination of both.  Text and binary files
 * are interconvertible using \c protoc.
 *
 * Our main configuration is of type config::Config, pointers to
 * this type are expected by the functions contained in here.
 *
 * @{ */

//! \brief merges the configuration from a text file.
//! A text format configuration file is read and the information
//! incorporated into an existing (possibly empty) configuration.
//! Multiple configurations can be merged by repeated calls.  If the file
//! doesn't exist, nothing is merged.  If anything else goes wrong, an
//! exception is thrown.
//!
//! \param filename filename of configuration file
//! \param mi configuration to be merged into
void merge_text_config( const std::string& filename, config::Config& mi ) ;

//! \brief merges the configuration from a binary file.
//! A binary format configuration file is read and the information
//! incorporated into an existing (possibly empty) configuration.
//! Multiple configurations can be merged by repeated calls.  If the file
//! doesn't exist, nothing is merged.  If anything else goes wrong, an
//! exception is thrown.
//!
//! \param filename filename of configuration file
//! \param mi configuration to be merged into
void merge_binary_config( const std::string& filename, config::Config& mi ) ;

//! \brief writes a configuration in text format.
//! The configuration is written to a file in Google's text format.  A
//! new file is always created and then moved to the correct name,
//! overwriting any preexisting file.  If anything goes wrong, an
//! exception is thrown.  Note that this not really the inverse of
//! merge_text_config() or merge_binary_config(), since reading a config
//! and then writing it unfolds the include options!
//!
//! \param filename filename of configuration file
//! \param mi pointer to configuration
void write_text_config( const std::string& filename, const config::Config& mi ) ;

//! \brief writes a configuration in binary format.
//! The configuration is written to a file in Google's native binary
//! format.  A new file is always created and then moved to the correct
//! name, overwriting any preexisting file.  If anything goes wrong, an
//! exception is thrown.  Note that this not really the inverse of
//! merge_text_config() or merge_binary_config(), since reading a config
//! and then writing it unfolds the include options!
//!
//! \param filename filename of configuration file
//! \param mi pointer to configuration
void write_binary_config( const std::string& filename, const config::Config& mi ) ;

//! @}
#endif

