#include "output.pb.h"
#include "outputfile.h"
#include "conffile.h"

#include "index.h"

#include <google/protobuf/io/zero_copy_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/text_format.h>

using namespace google::protobuf::io ;
using namespace google::protobuf ;
using namespace config ;
using namespace output ;

FileOutputStream out(1) ;

template< typename Msg > void print_msg( const Msg& m )
{
	TextFormat::Print( m, &out ) ;
	void* data ; int size ;
	out.Next( &data, &size ) ;
	*(char*)data = '\n' ;
	data = (char*)data + 1 ;
	--size ;
	if( !size ) out.Next( &data, &size ) ;
	*(char*)data = '\n' ;
	data = (char*)data + 1 ;
	--size ;
	if( size ) out.BackUp( size ) ;
}

typedef map< string, CompactGenome > Genomes ; 

void print_str( const std::string& s )
{
	int i = 0, j = s.size() ;
	while( i != j )
	{
		void *p0 ;
		int sz, k = 0 ;
		out.Next( &p0, &sz ) ;
		for( char *p = (char*)p0 ; k != sz && i != j ; ++k ) *p++ = s[i++] ;
		if( k != sz ) out.BackUp( sz - k ) ;
	}
}

void add_alignment( std::string::const_iterator qry, output::Hit &h, const config::Config &conf, Genomes &genomes )
{
	CompactGenome &g = genomes[ h.genome() ] ;
	if( !g.get_base() ) 
	{
		try { CompactGenome( h.genome(), conf, MADV_RANDOM ).swap( g ) ; }
		catch(...) { /* too bad, we'll make do without a genome */ }
	}

	if( g.get_base() ) 
	{
		std::string &r = *h.mutable_ref(), &q = *h.mutable_qry(), &c = *h.mutable_con() ;
		r.clear() ; q.clear() ; c.clear() ;
		DnaP ref = g.find_pos( h.sequence(), h.start_pos() ) ; 
		if( h.aln_length() < 0 ) ref = ref.reverse() + h.aln_length() + 1 ;

		for( size_t i = 0 ; i != h.cigar().size() ; ++i )
		{
			if( (uint8_t)h.cigar()[i] == 0 ) {
				r.push_back('~') ;
				q.push_back('~') ;
				c.push_back('~') ;
			}
			else if( (uint8_t)h.cigar()[i] < 128 ) for( size_t j = 0 ; j != (uint8_t)h.cigar()[i] ; ++j ) {
				r.push_back( from_ambicode( *ref ) ) ;
				q.push_back( *qry ) ;
				c.push_back( from_ambicode( *ref ) == *qry ? '*' : ' ' ) ;
				++ref; ++qry ;
			}
			else if( (uint8_t)h.cigar()[i] < 192 ) for( size_t j = 128 ; j != (uint8_t)h.cigar()[i] ; ++j, ++qry ) {
				r.push_back( '-' ) ;
				q.push_back( *qry ) ;
				c.push_back( ' ' ) ;
			}
			else for( size_t j = 192 ; j != (uint8_t)h.cigar()[i] ; ++j, ++ref ) {
				r.push_back( from_ambicode( *ref ) ) ;
				q.push_back( '-' ) ;
				c.push_back( ' ' ) ;
			}
		}
	}
	else
	{
		std::string &q = *h.mutable_qry() ;
		q.clear() ;

		for( size_t i = 0 ; i != h.cigar().size() ; ++i )
		{
			if( (uint8_t)h.cigar()[i] == 0 ) {
				q.push_back('~') ;
			}
			else if( (uint8_t)h.cigar()[i] < 128 ) for( size_t j = 0 ; j != (uint8_t)h.cigar()[i] ; ++j, ++qry ) {
				q.push_back( *qry ) ;
			}
			else if( (uint8_t)h.cigar()[i] < 192 ) for( size_t j = 128 ; j != (uint8_t)h.cigar()[i] ; ++j, ++qry ) {
				q.push_back( *qry ) ;
			}
			else for( size_t j = 192 ; j != (uint8_t)h.cigar()[i] ; ++j ) {
				q.push_back( '-' ) ;
			}
		}
	}
}

int main_( int argc, const char** argv )
{
	GOOGLE_PROTOBUF_VERIFY_VERSION ;
	config::Config conf ;
	Genomes genomes ;

	for( int argi = 1 ; argi != argc ; ++argi )
	{
		AnfoFile af( argv[argi] ) ;
		output::Header hdr = af.get_header() ;
		print_msg( hdr ) ;
		for( Result r ; af.read_result( r ) ; ) {
			if( r.has_best_to_genome() && r.best_to_genome().has_cigar() )
				add_alignment( r.sequence().begin(), *r.mutable_best_to_genome(), hdr.config(), genomes ) ;
			print_msg( r ) ;
		}
		print_msg( af.get_footer() ) ;
	}
	return 0 ;
}

