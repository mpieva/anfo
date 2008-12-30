#include "align.h"
#include "conffile.h"
#include "index.h"
#include "util.h"

using namespace std ;

//! \brief Rev-complements ASCII-encoded sequence.
//! This shouldn't be needed, it's only here for debugging.
//! \internal
void revcom( string &s )
{
	for( size_t i=0, j=s.size()-1 ; i<j ; ++i, --j )
	{
		char c = s[i] ;
		s[i] = s[j] ;
		s[j] = c ;
	}

	for( size_t i=0 ; i != s.length() ; ++i )
	{
		switch( s[i] )
		{
			case 'A': s[i] = 'T' ; break ;
			case 'C': s[i] = 'G' ; break ;
			case 'G': s[i] = 'C' ; break ;
			case 'T': s[i] = 'A' ; break ;
			case 'M': s[i] = 'K' ; break ;
			case 'K': s[i] = 'M' ; break ;
			case 'Y': s[i] = 'R' ; break ;
			case 'R': s[i] = 'Y' ; break ;
			case 'B': s[i] = 'V' ; break ;
			case 'H': s[i] = 'D' ; break ;
			case 'D': s[i] = 'H' ; break ;
			case 'V': s[i] = 'B' ; break ;
		}
	}
}


int main_( int argc, const char * argv[] )
{
	metaindex::Config mi ;
	merge_text_config( argc < 2 || !strcmp( argv[1], "R" ) ? "config.txt" : argv[1], mi ) ;
	metaindex::CompactIndex cix = find_compact_index( mi, argc < 3 ? "chr21" : argv[2] ) ;
	metaindex::Genome gdef = find_genome( mi, cix.genome_name() ) ;

	cout << "Found index " << cix.filename() << " with wordsize " << cix.wordsize() << " and " ;
	if( cix.has_cutoff() ) cout << "cutoff " << cix.cutoff() << '.' ;
	                  else cout << "no cutoff." ;
	cout << "\nFound genome " << gdef.filename() << " of length " << gdef.total_size() << endl ;

	FixedIndex ix( cix.filename().c_str(), cix.wordsize() ) ;
	CompactGenome g( gdef ) ;

	// string sq = argc < 4 ? "TAGGTCTTTTCCCAGGCCCAGTATCTGTGATTTGCTGTACATAACAGCTG" : argv[3] ;
	string sq = argc < 4 ? "AGVTMTTTTACCCAGGCCCAGTATCTGTGATTTGCTGTAGATAACGCTG" : argv[3] ;
	if( argc >= 2 && strcmp( argv[1], "R" ) == 0 ) revcom(sq) ;
	QSequence ps( sq.c_str() ) ;

	vector<Seed> seeds ;
	int num_raw = ix.lookup( ps, seeds ) ;
	int num_comb = seeds.size() ;
	select_seeds( seeds, /* ±d */ 2, /* ±o */ 16, /* m */ 12, g.get_contig_map() ) ;
	int num_clumps = seeds.size() ;

	cout << "got " << num_raw << " seeds, combined into " << num_comb
		 << " larger ones, clumped into " << num_clumps << " clumps."
		 << endl ;

	deque<flat_alignment> ol ;
	setup_alignments( g, ps, seeds.begin(), seeds.end(), ol ) ;
	flat_alignment best = find_cheapest( ol, INT_MAX, true ) ;

	cout << "Done near " << best.reference.abs() - g.get_base() << " costing " << best.penalty << ':' << endl ;

	deque< pair< flat_alignment, const flat_alignment* > > ol_ ;
	reset( best ) ;
	greedy( best ) ;
	(enter_bt<flat_alignment>( ol_ ))( best ) ;
	cout << find_cheapest( ol_ ) << endl ;

	return 0 ;
}

