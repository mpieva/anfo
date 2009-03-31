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
		else if( (uint8_t)cigar[i] < 128 ) s << (unsigned)(uint8_t)cigar[i]       << 'M' ;
		else if( (uint8_t)cigar[i] < 192 ) s << (unsigned)(uint8_t)cigar[i] - 128 << 'I' ;
		else                               s << (unsigned)(uint8_t)cigar[i] - 192 << 'D' ;
	}
    return s ;
}

// About len: since we're doing semi-global alignments, it is actually
// an error if the sequence length and the effective len from the cigar
// do not match.  However, files with such errors exist, so we need to
// filter out those faulty alignments.
int len_from_bin_cigar( const std::string& cigar ) {
	int l = 0 ;
	for( size_t i = 0 ; i != cigar.size() ; ++i )
	{
		if( (uint8_t)cigar[i] < 128 ) l += (uint8_t)cigar[i] ;
		else if( (uint8_t)cigar[i] < 192 ) l += (unsigned)(uint8_t)cigar[i] - 128 ;
	}
	return l ;
}

enum bad_stuff { goodness = 0, no_hit, multiple_hits, no_seqid, no_seq, bad_cigar, bad_stuff_max } ; 
const char *descr[] = { 0, "had no hit", "had multiple hits", "missed the sequence id"
	                  , "missed the sequence", "had a bad CIGAR" } ;

// About MAPQ: this is supposed to be the "quality of the mapping",
// which we interpret as confidence in the aligned position.  We can
// infer this from the score of the best and the second best alignment:
// if the difference is D, then the best alignment is exp(D/10) times
// more likely than the second best.  The SAM documentation hints that
// this is exactly the quantity they want to encode in MAPQ, so all we
// need to do is rescale it (from natural to decadic logarithm).  If
// there is no second hit, we'll assume a "perfect" mapping by writing
// out a 255.
bad_stuff protoHit_2_bam_Hit(output::Result &result){

    if (!result.has_best_hit() && !result.has_best_to_genome()) return no_hit ;
    if (result.has_best_hit() && result.has_best_to_genome()) return multiple_hits ;
    if (!result.has_seqid()) return no_seqid;
    if (!result.has_sequence()) return no_seq;

	// XXX 
	// Either one of those two is fine, depending on what we actually
	// want.  mixing them is probably wrong.  (But right now it's okay,
	// since only best_to_genome will ever be present.)
    output::Hit hit = (result.has_best_hit())?result.best_hit():result.best_to_genome();

	if (len_from_bin_cigar(hit.cigar()) != result.sequence().length()) return bad_cigar ;

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
    return goodness ;
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

    int discarded[bad_stuff_max] = {0};
    while (res.has_seqid()) {
        if (bad_stuff r = protoHit_2_bam_Hit(res))
            discarded[r]++;

        res = f_anfo.read_result();
    }

	for( int b = 1 ; b != bad_stuff_max ; ++b )
		if (discarded[b])
			std::cerr << discarded[b] << " reads " << descr[b] << std::endl;

    return (EXIT_SUCCESS);
}




