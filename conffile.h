#ifndef INCLUDED_CONFFILE_H
#define INCLUDED_CONFFILE_H

#include "metaindex.pb.h"

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

