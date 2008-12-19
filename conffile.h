#ifndef INCLUDED_CONFFILE_H
#define INCLUDED_CONFFILE_H

#include "metaindex.pb.h"

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
 * parts of the configuration resides in \c metaindex.proto, we
 * currently store the configuration in a single file using the protobuf
 * text format.  If the need arises, we can easily move to multiple
 * files, binary files, or a combination of both.  Text and binary files
 * are interconvertible using \c protoc.
 *
 * Our main configuration is of type metaindex::Config, pointers to
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
void merge_text_config( const std::string& filename, metaindex::Config& mi ) ;

//! \brief merges the configuration from a binary file.
//! A binary format configuration file is read and the information
//! incorporated into an existing (possibly empty) configuration.
//! Multiple configurations can be merged by repeated calls.  If the file
//! doesn't exist, nothing is merged.  If anything else goes wrong, an
//! exception is thrown.
//!
//! \param filename filename of configuration file
//! \param mi configuration to be merged into
void merge_binary_config( const std::string& filename, metaindex::Config& mi ) ;

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
void write_text_config( const std::string& filename, const metaindex::Config& mi ) ;

//! \brief finds a named genome.
//! The metaindex is scanned for a genome exactly matching the supplied
//! name.  If no genome is found, an exception is thrown.
//! \param mi the configuration
//! \param genome_name name to look for
//! \return A reference to the genome config.
const metaindex::Genome &find_genome( const metaindex::Config& mi, const std::string genome_name ) ;

//! \brief finds or creates a named genome.
//! The metaindex is scanned for a genome exactly matching the supplied
//! name.  If no genome is found, an empty one is generated instead.
//! \param mi the configuration
//! \param genome_name name to look for
//! \return A reference to the genome config.
metaindex::Genome& find_or_create_genome( metaindex::Config& mi, std::string genome_name ) ;

//! \brief finds or creates a compact index.
//! Looks for a compact index for the given named genome with the given
//! wordsize.  If none is found, a new one is created, but it still has
//! to be filled in.
//! \param mi the configuration
//! \param genome_name name of the indexed genome
//! \param wordsize length of indexed words
//! \return reference to the index configuration
metaindex::CompactIndex& find_or_create_compact_index(
		metaindex::Config& mi, const std::string genome_name, unsigned wordsize ) ;

//! \brief finds any index for a given genome
//! Looks for a compact index for the named genome.  The first index
//! with the right word size is returned.  If the requested word size is
//! zero, the first index is returned regardless of its word size.  An
//! exception is thrown if nothing appropriate is found.
//! \param mi the configuration
//! \param genome_name name of the indexed genome
//! \param wordsize desired length of indexed words
//! \return reference to the index configuration
const metaindex::CompactIndex& find_compact_index(
		const metaindex::Config& mi, const std::string genome_name, unsigned wordsize = 0 ) ; 

// }@
#endif

