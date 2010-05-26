#ifndef INCLUDED_CHAIN_H
#define INCLUDED_CHAIN_H

#include <istream>
#include <map>
#include <string>

// handling of UCSC liftover chains

// What we're interested in is actually a back-and-forth mapping between
// two genomes.  We compile the two chains down to a single table of
// blocks that pass all filters.  Here's that list:

struct Block
{
	// string p_chrom 		-- implicit, chromosome in primary genome
	// int start 			-- implicit, start position in primary genome
	int size : 30 ;			// matched size
	int multiple : 1 ;		// set if multiple blocks overlap
	int s_flip : 1 ;		// indicates reverse-complement
	int s_start ;			// start on secondary (negative for RC strand)
	std::string s_chromosome	;	// chromosome in secondary
	std::string chain_id ;
} ;

typedef std::map< int, Block > Blocks ;
typedef std::map< std::string, Blocks > Chains ;

struct DoubleCrossChain 
{
	std::string primary_genome ;
	std::string secondary_genome ;
	Chains chains ;
} ;


void fill_cross_chains(
		DoubleCrossChain &chain,
		std::string primary_genome, std::string secondary_genome,
		std::istream& primary_chain, const std::string& pfile,
		std::istream& secondary_chain, const std::string& sfile ) ;

#endif
