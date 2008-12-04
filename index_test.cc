#include "util.h"
#include "Index.h"
#include "conffile.h"
#include "Align.h"

using namespace std ;

int main_( int argc, const char * argv[] )
{
	if( argc < 3 ) 
	{
		std::clog << "Usage: " << argv[0] << " <conffile> <genome> [<sequence>]" << std::endl ;
		return 1 ;
	}

	Config cfg( argv[1] ) ;
	metaindex::CompactIndex cix = cfg.find_compact_index( argv[2] ) ;
	std::cout << "Found index " << cix.filename() << " with wordsize " 
		<< cix.wordsize() << " and " ;
	if( cix.has_cutoff() ) std::cout << "cutoff " << cix.cutoff() << '.' << std::endl ;
	else std::cout << "no cutoff." << std::endl ;

	FixedIndex ix( cix.filename().c_str(), cix.wordsize() ) ;

	vector<Seed> seeds ;
	ix.lookup( argc < 4 ? "TAGGTCTTTTCCCAGGCCCAGTATCTGTGATTTGCTGTACATAACAGCTG" : argv[3] , seeds ) ;
	cout << "got " << seeds.size() << " seeds, combined into " ;
	combine_seeds( seeds ) ;
	cout << seeds.size() << " larger ones, clumped into " ;
	select_seeds( seeds, /* ±d */ 2, /* ±o */ 16, /* m */ 24 ) ;
	cout << seeds.size() << " clumps." << endl ;

	for( unsigned i = 0 ; i != seeds.size() ; ++i )
		cout << seeds[i].offset << ' '
			 << seeds[i].diagonal << ' ' 
			 << seeds[i].size << endl ;

	// quick hack to init alignments... XXX
	return 0 ;
}

