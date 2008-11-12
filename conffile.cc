#include "conffile.h"
#include "util.h"

#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/text_format.h>

#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>

Config::Config( const std::string& filename ) : filename_( filename )
{
	int config_in = open( filename.c_str(), O_RDONLY ) ;
	if( config_in != -1 || errno != ENOENT )
	{
		throw_errno_if_minus1( config_in, "reading config" ) ;
		google::protobuf::io::FileInputStream fis( config_in ) ;
		google::protobuf::TextFormat::Parse( &fis, &mi_ ) ;
		close( config_in ) ;
	}
}

void Config::write() const
{
	std::string new_file_name = filename_ + '#' ;
	int config_out = throw_errno_if_minus1(
			creat( new_file_name.c_str(), 0644 ), "writing config" ) ;
	google::protobuf::io::FileOutputStream fos( config_out ) ;
	fos.SetCloseOnDelete( true ) ;
	google::protobuf::TextFormat::Print( mi_, &fos ) ;

	throw_errno_if_minus1( 
			rename( new_file_name.c_str(), filename_.c_str() ),
			"renaming config" ) ;
}

const metaindex::Genome& Config::find_genome( const std::string genome_name ) const
{
	for( int i = 0 ; i != mi_.genome_size() ; ++i )
		if( mi_.genome(i).name() == genome_name )
			return mi_.genome(i) ;
	throw "don't know about genome " + genome_name ;
}

metaindex::Genome *Config::find_or_create_genome( const std::string genome_name )
{
	for( int i = 0 ; i != mi_.genome_size() ; ++i )
		if( mi_.genome(i).name() == genome_name )
			return mi_.mutable_genome(i) ;
	return mi_.add_genome() ;
}

metaindex::CompactIndex *Config::find_or_create_compact_index(
		const std::string genome_name, unsigned wordsize )
{
	for( int i = 0 ; i != mi_.compact_index_size() ; ++i )
		if( mi_.compact_index(i).genome_name() == genome_name 
				&& mi_.compact_index(i).wordsize() == wordsize )
			return mi_.mutable_compact_index(i) ;
	return mi_.add_compact_index() ;
}

const metaindex::CompactIndex& Config::find_compact_index(
		const std::string genome_name ) const
{
	for( int i = 0 ; i != mi_.compact_index_size() ; ++i )
		if( mi_.compact_index(i).genome_name() == genome_name )
			return mi_.compact_index(i) ;
	throw "don't have an index for genome " + genome_name ;
}

