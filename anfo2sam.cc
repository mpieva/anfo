/* 
 * File:   anfo2sam.cc
 * Author: michael_siebauer
 *
 * Created on March 13, 2009, 1:58 PM
 */

#include <cmath>
#include <stdlib.h>
#include "output.pb.h"
#include <iostream>
#include <stdint.h>
#include "outputfile.h"
#include <string>
#include "util.h"

// BAM FLAGS
#define BAM_FPAIRED        1   // read is paired in sequencing
#define BAM_FPROPER_PAIR   2   // read is mapped in proper pair
#define BAM_FUNMAP         4   // query seq. is unmapped
#define BAM_FMUNMAP        8   // mate is umapped
#define BAM_FREVERSE      16   // strand of query (0 - fwd, 1 - rev)
#define BAM_FMREVERSE     32   // strand of mate
#define BAM_FREAD1        64   // read is 1. in pair
#define BAM_FREAD2       128   // read is 2. in pair
#define BAM_FSECONDARY   256   // alignment is NOT primary
#define BAM_FQCFAIL      512   // read fails due low quality
#define BAM_FDUP        1024   // read is duplicate


std::ostream& decode_binCigar(std::ostream& s, const std::string &cigar) {
	for( size_t i = 0 ; i != cigar.size() ; ++i )
	{
		if( (uint8_t)cigar[i] == 0 ) continue;
		else if( (uint8_t)cigar[i] < 128 ) s << (unsigned)cigar[i]       << 'M' ;
		else if( (uint8_t)cigar[i] < 192 ) s << (unsigned)cigar[i] - 128 << 'I' ;
		else                               s << (unsigned)cigar[i] - 192 << 'D' ;
	}
    return s ;
}

// About MAPQ: this is supposed to be the "quality of the mapping",
// which we interpret as confidence in the aligned position.  We can
// infer this from the score of the best and the second best alignment:
// if the difference is D, then the best alignment is exp(D/10) times
// more likely than the second best.  The SAM documentation hints that
// this is exactly the quantity they want to encode in MAPQ, so all we
// need to do is rescale it (from natural to decadic logarithm).  If
// there is no second hit, we'll assume a "perfect" mapping by writing
// out a 255.
int protoHit_2_bam_Hit(output::Result &result){

    // I need exactly one hit   (!XOR)
    if (result.has_best_hit() == result.has_best_to_genome())
        return (EXIT_FAILURE);

    if (!result.has_seqid())
        return (EXIT_FAILURE);

    if (!result.has_sequence())
        return (EXIT_FAILURE);

	// XXX 
	// Either one of those two is fine, depending on what we actually
	// want.  mixing them is probably wrong.  (But right now it's okay,
	// since only best_to_genome will ever be present.)
    output::Hit hit = (result.has_best_hit())?result.best_hit():result.best_to_genome();

/*QNAME*/   std::cout << result.seqid() << "\t";
/*FLAG */   std::cout << ((hit.aln_length() < 0)?(BAM_FREVERSE):(0)) << "\t";  // TODO: calc flag
/*RNAME*/   std::cout << hit.sequence() << "\t";
/*POS*/     std::cout << ( (hit.aln_length() >= 0)?hit.start_pos():(hit.start_pos() + hit.aln_length() + 1) ) << "\t";
/*MAPQ*/   	std::cout << ( result.has_diff_to_next() 
					? (int)( 0.5 + result.diff_to_next() / std::log(10.0) ) : 255 ) << '\t' ;
/*CIGAR*/   decode_binCigar(std::cout, hit.cigar()) << "\t";
          // We don't have paired end reads
/*MRNM*/    std::cout << "*" << "\t";
/*MPOS*/    std::cout << "0" << "\t";
/*ISIZE*/   std::cout << "0" << "\t";

/*SEQ*/     std::cout << result.sequence() << "\t";
/*QUAL*/    if( result.has_quality() )
				for (size_t i = 0; i < result.quality().size(); i++)
					std::cout << char((uint8_t)result.quality()[i] + 33);
			else
				std::cout << '*' ;

/*[TAGS]*/

    std::cout << std::endl;
    return (EXIT_SUCCESS);
}



int main_( int argc, const char * argv[] ){
    GOOGLE_PROTOBUF_VERIFY_VERSION;

    if (argc <= 1){
       std::cerr << "Convert .anfo alignments to SAM/BAM format" << std::endl;
       std::cerr << "Usage: " << argv[0] << " <anfo-alignment> " << std::endl;
       return (EXIT_FAILURE);
    }

    AnfoFile f_anfo(argv[1]);
    f_anfo.read_header();
    output::Result res = f_anfo.read_result();

    int discarded = 0;
    while (res.has_seqid()) {
        if (protoHit_2_bam_Hit(res) != EXIT_SUCCESS)
            discarded++;

        res = f_anfo.read_result();
    }

    if (discarded)
        std::cerr << discarded << " reads could not been converted!" << std::endl;

    return (EXIT_SUCCESS);
}




