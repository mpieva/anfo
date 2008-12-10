#include "align.h"
#include "conffile.h"
#include "index.h"
#include "util.h"

using namespace std ;

void revcom( std::string &s )
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
		}
	}
}


int main_( int argc, const char * argv[] )
{
	Config cfg( argc < 2 ? "index.txt" : strcmp( argv[1], "R" ) ? argv[1] : "index.txt" ) ;
	metaindex::CompactIndex cix = cfg.find_compact_index( argc < 3 ? "chr21" : argv[2] ) ;
	std::cout << "Found index " << cix.filename() << " with wordsize " 
		<< cix.wordsize() << " and " ;
	if( cix.has_cutoff() ) std::cout << "cutoff " << cix.cutoff() << '.' << std::endl ;
	else std::cout << "no cutoff." << std::endl ;

	FixedIndex ix( cix.filename().c_str(), cix.wordsize() ) ;

	std::string s = argc < 4 ? "TAGGTCTTTTCCCAGGCCCAGTATCTGTGATTTGCTGTACATAACAGCTG" : argv[3] ;
	// std::string s = argc < 4 ? "AGVTMTTTTACCCAGGCCCAGTATCTGTGATTTGCTGTAGATAACGCTG" : argv[3] ;
	if( argc >= 2 && strcmp( argv[1], "R" ) == 0 ) revcom(s) ;

	vector<Seed> seeds ;
	ix.lookup( s, seeds ) ;

	cout << "got " << seeds.size() << " seeds, combined into " ;
	combine_seeds( seeds ) ;
	cout << seeds.size() << " larger ones, clumped into " ;
	select_seeds( seeds, /* ±d */ 2, /* ±o */ 16, /* m */ 24 ) ;
	cout << seeds.size() << " clumps." << endl ;

	metaindex::Genome gdef = cfg.find_genome( cix.genome_name() ) ;
	std::cout << "Found genome " << gdef.filename() << " of length "
		<< gdef.total_size() << std::endl ;

	CompactGenome g( gdef.filename().c_str() ) ;
	std::cout << "Loaded genome." << std::endl ;

	PreparedSequence ps( s.c_str() ) ;
	std::deque< flat_alignment > ol ;

	// quick hack to init alignments... XXX: move this somewhere
	// sensible... as soon as it works.
	for( std::vector<Seed>::const_iterator s = seeds.begin() ; s != seeds.end() ; ++s )
	{
		std::clog << *s << std::endl ;

		flat_alignment fa ;
		if( s->offset >= 0 ) {
			fa.reference = g.get_base() + s->offset + s->diagonal + s->size / 2 ;
			fa.query = ps.forward() + s->offset + s->size / 2 ;
		}
		else
		{
			// XXX I want to start in the middle of the seed.  Why
			// does this not work on the reverse-complement strand?
			fa.reference = g.get_base() + s->offset + s->diagonal ; // + s->size / 2 ;
			fa.query = ps.reverse() - s->offset ; // + s->size / 2 ;
		}
		reset( fa ) ;
		greedy( fa ) ;
		ol.push_back( fa ) ;
	}
	std::make_heap( ol.begin(), ol.end() ) ;
	flat_alignment best = find_cheapest( ol ) ;

	std::cout << "Done near " << best.reference - g.get_base()
		      << " costing " << best.penalty << std::endl ;

	std::deque< std::pair< flat_alignment, const flat_alignment* > > ol_ ;
	reset( best ) ;
	greedy( best ) ;
	ol_.push_back( std::make_pair( best, (const flat_alignment*)0 ) ) ;
	std::cout << find_cheapest( ol_ ) << std::endl ;

	return 0 ;
}

