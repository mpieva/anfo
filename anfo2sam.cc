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


template <typename Iter> 
std::ostream& decode_binCigar(std::ostream& s, Iter begin, Iter end ) {
	for( ; begin != end ; ++begin )
	{
		if( (uint8_t)(*begin) == 0 ) continue;
		else if( (uint8_t)(*begin) < 128 ) s << (unsigned)(uint8_t)(*begin)       << 'M' ;
		else if( (uint8_t)(*begin) < 192 ) s << (unsigned)(uint8_t)(*begin) - 128 << 'I' ;
		else                               s << (unsigned)(uint8_t)(*begin) - 192 << 'D' ;
	}
    return s ;
}

// About len: since we're doing semi-global alignments, it is actually
// an error if the sequence length and the effective len from the cigar
// do not match.  However, files with such errors exist, so we need to
// filter out those faulty alignments.
unsigned len_from_bin_cigar( const std::string& cigar ) {
	unsigned l = 0 ;
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

//! \page anfo_to_sam Conversion to SAM
//! \c anfo2sam converts a single ANFO file into a single SAM file.  If
//! the ANFO file is sorted by genome coordinate, the resulting SAM file
//! will be sorted, too, and be labelled correctly, so it can directly
//! be converted to BAM and be indexed. A few decisions had to be made
//! on the interpretation of SAM fields.  
//!
//! About MAPQ: this is supposed to be the "quality of the mapping",
//! which should apparently be interpreted as confidence in the aligned
//! position.  We can infer this from the score of the best and the
//! second best alignment: if the difference is D, then the best
//! alignment is exp(D/10) times more likely than the second best.  The
//! SAM documentation hints that this is exactly the quantity they want
//! to encode in MAPQ, so all we need to do is rescale it (from natural
//! to decadic logarithm).  If there is no second hit, we'll assume a
//! "perfect" mapping by writing out a 255.
//!
//! \todo In reality, the map quality would only be the score cutoff,
//!       which is the cutoff from the header (appropriate to this
//!       particular sequence) scaled by the sequence length.  In the
//!       current setup, only sequences of length 40 or above really
//!       achieve a MAPQ of 255.
//!
//! About the hit to convert: ANFO is happy represernting more than one
//! hit per sequence, but SAM is not.  We will convert the \c
//! best_to_genome entry, ignoring everything else.  This is both the
//! most useful alignment and the one best in line with the intended use
//! of SAM.
//!
//! About the alignment score: This optional field of SAM is filled with
//! the raw ANFO alignment score, which is a bit difficult to interpret
//! if you're used to other aligners.  A lower score means a better
//! alignment, and a score of (length-20) * 18 is about the quality that
//! would give a zero score in BLAST.  Since SAM doesn't specify
//! anything about the alignment, we'll leave it that way.
//!
//! \todo calculate other SAM/BAM flags

bad_stuff protoHit_2_bam_Hit(output::Result &result){

    if (!result.has_best_to_genome()) return no_hit ;
    if (!result.has_seqid()) return no_seqid ;
    if (!result.has_sequence()) return no_seq ;

    output::Hit hit = result.best_to_genome() ;

	if (len_from_bin_cigar(hit.cigar()) != result.sequence().length()) return bad_cigar ;

/*QNAME*/   std::cout << result.seqid() << "\t";
/*FLAG */   std::cout << (hit.aln_length() < 0 ? BAM_FREVERSE : 0) << "\t";
/*RNAME*/   std::cout << hit.sequence() << "\t";
/*POS*/     std::cout << 1+hit.start_pos() << "\t";
/*MAPQ*/   	std::cout << ( result.has_diff_to_next() 
					? (int)( 0.5 + result.diff_to_next() / std::log(10.0) ) : 255 ) << '\t' ;
/*CIGAR*/   ( hit.aln_length() >= 0 ? decode_binCigar(std::cout, hit.cigar().begin(), hit.cigar().end() )
			                        : decode_binCigar(std::cout, hit.cigar().rbegin(), hit.cigar().rend() ) ) << "\t";
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
/*SCORE*/	std::cout << "\tAS:i:" << hit.score() << '\n' ;

    return goodness ;
}



int main_( int argc, const char * argv[] ){
    GOOGLE_PROTOBUF_VERIFY_VERSION;

    if (argc <= 1){
       std::cerr << "Convert .anfo alignments to SAM/BAM format" << std::endl;
       std::cerr << "Usage: " << argv[0] << " <anfo-alignment> " << std::endl;
       return (EXIT_FAILURE);
    }

    int discarded[bad_stuff_max] = {0};
    AnfoFile f_anfo(argv[1]);
	output::Header hdr = f_anfo.get_header();
	std::cout << "@HD\tVN:1.0" ;
	if( hdr.is_sorted_by_coordinate() ) std::cout << "\tSO:coordinate" ;
	else if( hdr.is_sorted_by_name() ) std::cout << "\tSO:queryname" ;
	std::cout << "\n@PG\tID:ANFO\tVN:" << hdr.version() << '\n' ;

    for( output::Result res ; f_anfo.read_result( res ) ; ) 
        if (bad_stuff r = protoHit_2_bam_Hit(res))
            discarded[r]++;

	for( int b = 1 ; b != bad_stuff_max ; ++b )
		if (discarded[b])
			std::cerr << discarded[b] << " reads " << descr[b] << std::endl;

    return (EXIT_SUCCESS);
}




