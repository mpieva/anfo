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

typedef map< string, CompactGenome > Genomes ; 

void add_alignment( std::string::const_iterator qry, output::Hit &h, const config::Config &conf, Genomes &genomes )
{
	CompactGenome &g = genomes[ h.genome_file() ] ;
	if( !g.get_base() ) CompactGenome( h.genome_file(), conf, MADV_RANDOM ).swap( g ) ;

	if( g.get_base() ) 
	{
		std::string &r = *h.mutable_ref(), &q = *h.mutable_qry(), &c = *h.mutable_con() ;
		r.clear() ; q.clear() ; c.clear() ;
		DnaP ref = g.find_pos( h.sequence(), h.start_pos() ) ; 
		if( h.aln_length() < 0 ) ref = ref.reverse() + h.aln_length() + 1 ;

		for( size_t i = 0 ; i != h.cigar().size() ; ++i )
		{
			if( (uint8_t)h.cigar()[i] < 128 ) for( size_t j = 0 ; j != (uint8_t)h.cigar()[i] ; ++j ) {
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
			if( (uint8_t)h.cigar()[i] < 128 ) for( size_t j = 0 ; j != (uint8_t)h.cigar()[i] ; ++j, ++qry ) {
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
	Genomes genomes ;

	if( argc <= 3 ) {
		std::clog << "Usage: " << argv[0] << " <slope> <offset> [<file> ...]\n"
			         "  Extracts alignments with score < slope * (length - offset) from ANFO\n"
					 "  files and writes alignments in multi-FASTA format to stdout" << std::endl ; 
		return 1 ;
	}

	float slope = atof( argv[1] ) ;
	int  offset = atoi( argv[2] ) ;

	for( const char** arg = argv+3 ; arg != argv+argc ; ++arg )
	{
		AnfoFile af( *arg ) ;
		output::Header hdr = af.get_header() ;
		for( Result r ; af.read_result( r ) ; ) {
			if( r.has_best_to_genome() && r.best_to_genome().has_cigar() )
			{
				int leff = r.has_trim_right() ? r.trim_right() : r.sequence().size() ;
				if( r.best_to_genome().score() <= slope * (leff - offset) ) 
				{
					add_alignment( r.sequence().begin(), *r.mutable_best_to_genome(), hdr.config(), genomes ) ;
					std::cout << '>' << r.best_to_genome().sequence() << ' '
						<< r.best_to_genome().start_pos()
						<< (r.best_to_genome().aln_length() > 0 ? '+' : '-')
						<< r.best_to_genome().start_pos() + abs(r.best_to_genome().aln_length()) - 1
						<< '\n' << r.best_to_genome().ref() << '\n' 
						<< '>' << r.seqid() 
						<< ( r.has_trim_right() ? " adapter cut off\n" : "\n" ) 
						<< r.best_to_genome().qry() << std::endl ;
				}
			}
		}
	}
	return 0 ;
}
