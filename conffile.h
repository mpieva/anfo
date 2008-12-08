#ifndef INCLUDED_CONFFILE_H
#define INCLUDED_CONFFILE_H

#include "metaindex.pb.h"

/*! \brief Convenient access to configuration files.
 *
 * Configuration files contain mostly metadata, e.g. which genomes are
 * known, what contigs do they consist of, how were they indexed.  It
 * might be useful to add other stuff, such as seeding policies and
 * alignment parameters.
 *
 * In order to keep configuration files in a sensible format and
 * extensible, we abuse Google's protocol buffers
 * (http://code.google.com/p/protobuf/) for this.  The definition of the
 * possible parts of the configuration resides in \c metaindex.proto, we
 * currently store the configuration in a single file using the protobuf
 * text format.  If the need arises, we can easily move to multiple
 * files, binary files, or a combination of both.  Text and binary files
 * are interconvertible using \c protoc.
 */
class Config
{
	private:
		std::string filename_ ;
		metaindex::MetaIndex mi_ ;

	public:
		Config( const std::string& filename ) ;

		const metaindex::MetaIndex *operator -> () const { return &mi_ ; }
		metaindex::MetaIndex *operator -> () { return &mi_ ; }

		const metaindex::Genome& find_genome( const std::string genome_name ) const ;
		metaindex::Genome *find_or_create_genome( const std::string genome_name ) ;
		metaindex::CompactIndex *find_or_create_compact_index( const std::string genome_name, unsigned wordsize ) ;
		const metaindex::CompactIndex& find_compact_index( const std::string genome_name ) const ;

		void write() const ;
} ;

#endif

