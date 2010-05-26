#include "../config.h"

#include "chain.h"
#include "util.h"

#include <cassert>
#include <sstream>

static void put_block( Blocks& blocks, int tStart, int match, int qStart, std::string qName, std::string cid )
{
	bool flip = 0 ;

	// always operate on forward target strand
	if( tStart < 0 ) {
		tStart = -tStart - match + 1 ;
		flip = 1 ;
	}

	if( qStart < 0 ) {
		qStart = -qStart - match + 1 ;
		flip ^= 1 ;
	}

	// iterate over potentially overlapping block.
	// start by finding the earliest block that doesn't overlap:
	for( Blocks::iterator b = blocks.lower_bound( tStart + match ) ; match ; )
	{
		// where are we now?  
		// 1) if we don't have an appropriate block or the one we have
		// doesn't overlap the new block, we can store it and are
		// finished.
		if( b == blocks.begin() ||
				( --b, b->first + b->second.size <= tStart ) )
		{
			Block& nb = blocks[ tStart ] ;
			nb.size = match ;
			nb.s_start = qStart ;
			nb.multiple = 0 ;
			nb.s_flip = flip ;
			nb.s_chromosome = qName ;
			nb.chain_id = cid ;
			break ;
		}

		--b ;
		// we have a block now, and it is guaranteed to overlap
		// 2) if the old block is contained in the new one, we mark it
		//    as duplicate, store a remainder of the new block (if
		//    any), reduce the new block's size and continue
		if( b->first >= tStart && b->first + b->second.size <= tStart + match )
		{
			console.output( Console::notice, "case 2" ) ;
			b->second.multiple = 1 ;
			if( b->first + b->second.size < tStart + match )
			{
				Block& nb = blocks[ b->first + b->second.size ] ;
				nb = b->second ;
				nb.size = (tStart + match) - (b->first + b->second.size) ;
				nb.s_start += b->second.size ;
			} 
			match = b->first - tStart ;
		}

		// 3) the old block overlaps us to the right.  we split it, mark
		//    one part as multiple, reduce the new block's size and go on.
		else if( b->first >= tStart )
		{
			console.output( Console::notice, "case 3" ) ;
			int d = tStart + match - b->first ;
			Block& nb = blocks[ tStart + match ] ;
			nb = b->second ;
			nb.size -= d ;
			nb.s_start += d ;

			b->second.size = d ;
			b->second.multiple = 1 ;

			match -= d ;
		}

		// 4) the old block overlaps us to the left.  we split it, store
		//    what we have left and finish.
		else if( b->first + b->second.size <= tStart + match ) 
		{
			console.output( Console::notice, "case 4" ) ;
			int d = (tStart + match) - tStart ;

			Block& nb = blocks[ tStart ] ;
			nb = b->second ;
			nb.size = d ;
			nb.s_start += b->second.size - d ;
			nb.multiple = 1 ;

			Block& mb = blocks[ b->first + b->second.size ] ;
			mb.size = match - d ;
			mb.s_start = qStart + d ;
			mb.multiple = 0 ;
			mb.s_flip = flip ;
			mb.s_chromosome = qName ;
			mb.chain_id = cid ;

			b->second.size -= d ;
			break ;
		}

		// 5) we are contained in the old block.  we split it in three,
		//    mark the middle one duplicated and finish.
		else
		{
			console.output( Console::notice, "case 5" ) ;
			int r = (b->first + b->second.size) - (tStart + match) ;

			Block& rb = blocks[ tStart + match ] ;
			rb = b->second ;
			rb.size = r ;
			rb.s_start += b->second.size - r ;

			Block& mb = blocks[ tStart ] ;
			mb.size = match ;
			mb.s_start = qStart ;
			mb.multiple = 1 ;
			mb.s_flip = flip ;
			mb.s_chromosome = qName ;
			mb.chain_id = cid ;

			b->second.size = tStart - b->first ;
			break ;
		}
	}
}

static void read_chain_file( std::istream &is, const std::string& name, Chains& chains )
{
    console.output( Console::notice, "reading " + name ) ;
	std::string line, score, key, tName, qName, tStrand, qStrand, cid ;
	int tSize, tStart, tEnd, qSize, qStart, qEnd, nchains = 0, nblocks = 0 ;

	if( getline( is, line ) ) for(;;)
	{
		std::stringstream ss( line ) ;
		// read a line, assuming (and then checking) that it is a chain
		// header
		if( ss >> key >> score >> tName >> tSize >> tStrand >> tStart
				>> tEnd >> qName >> qSize >> qStrand >> qStart >> qEnd >> cid
				&& key == "chain" )
		{
			++nchains ;
			assert( 0 <= tStart && tStart <= tEnd ) ;
			assert( 0 <= qStart && qStart <= qEnd ) ;
			assert( tStrand == "+" || tStrand == "-" ) ;
			assert( qStrand == "+" || qStrand == "-" ) ;
			if( tStrand == "-" ) tStart = -tEnd + 1 ;
			if( qStrand == "-" ) qStart = -qEnd + 1 ;

			Blocks& blocks = chains[ tName ] ;

			// read a line, it should define a block
			while( getline( is, line ) ) 
			{
				// try parsing as block definition
				int match, tGap = 0, qGap = 0 ;
				std::stringstream ss( line ) ;
				if( !(ss >> match) ) break ; // get out of here, but keep the line
				ss >> tGap >> qGap ;
				++nblocks ;

				put_block( blocks, tStart, match, qStart, qName, cid ) ;

				tStart += match + tGap ;
				qStart += match + qGap ;
			}
		}
		else if( !getline( is, line ) ) break ;
	}

	int gblocks = 0 ;
	for( Chains::iterator c = chains.begin() ; c !=chains.end() ; )
	{
		Blocks& blocks = c->second ;
		for( Blocks::iterator b = blocks.begin() ; b != blocks.end() ; )
		{
			Blocks::iterator b_ = b ; 
			++b_ ;

			if( b->second.multiple ) blocks.erase( b ) ;
			else ++gblocks ;

			b = b_ ;
		}

		Chains::iterator c_ = c ;
		++c_ ;
		if( blocks.empty() ) chains.erase( c ) ;
		c = c_ ;
	}

	std::stringstream ss ;
	ss  << nchains << " chains w/ " << nblocks << " blocks, "
		<< gblocks << " of which are good" ;
    console.output( Console::notice, ss.str() ) ;
}

// XXX hrmpf.
static void make_invertible( Chains& cs, Chains& ics )
{
	int chrom = 0, no_block = 0, perfect = 0, shrug = 0 ;
	for( Chains::iterator c = cs.begin() ; c != cs.end() ; ++c )
	{
		Blocks& blocks = c->second ;
		for( Blocks::iterator b = blocks.begin() ; b != blocks.end() ; )
		{
			Blocks& ic = ics[ b->second.s_chromosome ] ; 
			Blocks::iterator bh = ic.lower_bound( b->second.s_start + b->second.size ) ;
			Blocks::iterator b_ = b ;
			++b_ ;

			bool good = false ;
			if( bh == ic.begin() ) ++no_block ;
			else {
				--bh ;
				if( bh->first + bh->second.size <= b->first ) ++no_block ;
				else if( bh->second.s_chromosome != c->first ) ++chrom ;
				else if( bh->second.s_start == b->first && bh->second.size == b->second.size
						&& bh->second.s_flip == b->second.s_flip ) { good = true ; ++perfect ; }
				else {
					++shrug ;
					std::ostringstream ss ;
					ss << "Weird case: (hg18)" << c->first << ':' << b->first << '+' << b->second.size
						<< " --> (pt2)" << b->second.s_chromosome << ':' << b->second.s_start 
						<< " --> (pt2)" << bh->first << '+' << bh->second.size 
						<< " --> (hg18)" << bh->second.s_chromosome << ':' << bh->second.s_start 
						<< '+' << bh->second.size ;
					console.output( Console::notice, ss.str() ) ;
				}
			}

			if( !good ) blocks.erase( b ) ;
			b = b_ ;
		}
	}

	std::ostringstream ss ;
	ss << no_block << " without equivalent, " << chrom << " blocks on wrong chromosome, "
		<< perfect << " perfect matches, " << shrug << " other cases." ;
	console.output( Console::notice, ss.str() ) ;
}

void fill_cross_chains(
		DoubleCrossChain &chain,
		std::string primary_genome, std::string secondary_genome,
		std::istream& primary_chain, const std::string& pfile,
		std::istream& secondary_chain, const std::string& sfile )
{
	chain.primary_genome = primary_genome ;
	chain.secondary_genome = secondary_genome ;

	Chains ichains ;
	read_chain_file( primary_chain, pfile, chain.chains ) ;
	read_chain_file( secondary_chain, sfile, ichains ) ;

	make_invertible( chain.chains, ichains ) ;
}

	
