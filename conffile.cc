#include "conffile.h"
#include "util.h"

#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/text_format.h>

#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>

using namespace std ;
using namespace metaindex ;

static void unfold_includes( Config &mi )
{
	for( int i = 0 ; i != mi.load_size() ; ++i )
		merge_binary_config( mi.load(i), mi ) ;

	for( int i = 0 ; i != mi.include_size() ; ++i )
		merge_text_config( mi.include(i), mi ) ;
}

void merge_text_config( const string& filename, Config &mi )
{
	Config c ;
	int config_in = open( filename.c_str(), O_RDONLY ) ;
	if( config_in != -1 || errno != ENOENT )
	{
		throw_errno_if_minus1( config_in, "reading config" ) ;
		google::protobuf::io::FileInputStream fis( config_in ) ;
		if( !google::protobuf::TextFormat::Parse( &fis, &c ) )
			throw "parse error reading " + filename ;
		close( config_in ) ;
	}
	unfold_includes( c ) ;
	mi.MergeFrom( c ) ;
}

void merge_binary_config( const string& filename, Config &mi )
{
	Config c ;
	int config_in = open( filename.c_str(), O_RDONLY ) ;
	if( config_in != -1 || errno != ENOENT )
	{
		throw_errno_if_minus1( config_in, "reading config" ) ;
		if( !c.ParseFromFileDescriptor( config_in ) )
			throw "unmarshal error reading " + filename ;
		close( config_in ) ;
	}
	unfold_includes( c ) ;
	mi.MergeFrom( c ) ;
}

void write_text_config( const string& filename, const Config &mi )
{
	string new_file_name = '#' + filename + '#' ;
	int config_out = throw_errno_if_minus1(
			creat( new_file_name.c_str(), 0644 ), "writing config" ) ;
	google::protobuf::io::FileOutputStream fos( config_out ) ;
	fos.SetCloseOnDelete( true ) ;
	google::protobuf::TextFormat::Print( mi, &fos ) ;

	throw_errno_if_minus1( 
			rename( new_file_name.c_str(), filename.c_str() ),
			"renaming config" ) ;
}

const Genome& find_genome( const Config &mi, const string genome_name ) 
{
	for( int i = 0 ; i != mi.genome_size() ; ++i )
		if( mi.genome(i).name() == genome_name )
			return mi.genome(i) ;
	throw "don't know about genome " + genome_name ;
}

Genome& find_or_create_genome( Config &mi, const string genome_name )
{
	for( int i = 0 ; i != mi.genome_size() ; ++i )
		if( mi.genome(i).name() == genome_name )
			return *mi.mutable_genome(i) ;
	return *mi.add_genome() ;
}

CompactIndex& find_or_create_compact_index(
		Config &mi, const string genome_name, unsigned wordsize )
{
	for( int i = 0 ; i != mi.compact_index_size() ; ++i )
		if( mi.compact_index(i).genome_name() == genome_name 
				&& mi.compact_index(i).wordsize() == wordsize )
			return *mi.mutable_compact_index(i) ;
	return *mi.add_compact_index() ;
}

const CompactIndex& find_compact_index( const
		Config &mi, const string genome_name, unsigned wordsize )
{
	for( int i = 0 ; i != mi.compact_index_size() ; ++i )
		if( mi.compact_index(i).genome_name() == genome_name
				&& ( !wordsize || mi.compact_index(i).wordsize() == wordsize) )
			return mi.compact_index(i) ;
	throw "don't have an index for genome " + genome_name ;
}

