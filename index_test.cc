#include "util.h"
#include "Index.h"
#include "conffile.h"

using namespace std ;

int main_( int argc, const char * argv[] )
{
	if( argc < 2 ) return 1 ;
	Config cfg( argv[1] ) ;
	metaindex::CompactIndex cix = cfg.find_compact_index( argv[2] ) ;

	FixedIndex ix( cix.filename().c_str(), cix.wordsize() ) ;

	vector<Seed> seeds ;
	ix.lookup( "TAGGTCTTTTCCCAGGCCCAGTATCTGTGATTTGCTGTACATAACAGCTG", seeds ) ;
	cout << "got " << seeds.size() << " seeds, combined into " ;
	combine_seeds( seeds ) ;
	cout << seeds.size() << " larger ones, clumped into " ;
	select_seeds( seeds, 2, 16, 20 ) ;
	cout << seeds.size() << " clumps." << endl ;

	for( unsigned i = 0 ; i != seeds.size() ; ++i )
		cout << seeds[i].offset << ' '
			 << seeds[i].diagonal << ' ' 
			 << seeds[i].size << endl ;

	return 0 ;
}

