#include "../config.h"

#include "chain.h"
#include "util.h"

#include <fstream>

const char hg_pt[] = "/mnt/sequencedb/ucsc/goldenPath/hg18/liftOver/hg18ToPanTro2.over.chain" ;
const char pt_hg[] = "/mnt/sequencedb/ucsc/goldenPath/panTro2/liftOver/panTro2ToHg18.over.chain" ;

int main()
{
	console.more_verbose() ;
	DoubleCrossChain c ;
	std::ifstream f1( hg_pt ), f2( pt_hg ) ;
	fill_cross_chains( c, "hg18", "pt2", f1, hg_pt, f2, pt_hg ) ;
}
