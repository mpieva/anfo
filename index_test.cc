#include "util.h"
#include "Index.h"

using namespace std ;

int main_( int argc, const char * const argv[] )
{
	if( argc < 2 ) return 1 ;
	FixedIndex ix( argv[1], 10 ) ;

	vector<Seed> seeds ;
	ix.lookup( "TAGGTCTTTTCCCAGGCCCAGTATCTGTGATTTGCTGTACATAACAGCTG", seeds ) ;
	cout << "got " << seeds.size() << " seeds, combined into " ;
	combine_seeds( seeds ) ;
	cout << seeds.size() << " larger ones, clumped into " ;
	select_seeds( seeds, 2, 2, 20 ) ;
	cout << seeds.size() << " clumps." << endl ;

	for( unsigned i = 0 ; i != seeds.size() ; ++i )
		cout << seeds[i].offset << ' '
			 << seeds[i].diagonal << ' ' 
			 << seeds[i].size << endl ;

	return 0 ;
}

