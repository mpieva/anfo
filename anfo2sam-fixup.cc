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

// This is a VERY UGLY fix for a very ugly coordinate mishap.  Versions
// 0.8.0 and before of ANFO miscalculated the position of some genomic
// contigs.  We're trying to fix this here for hg18 by adding offsets to
// the coordinates.  This is a band-aid, will probably not work reliably
// and chr{8,17,X}_random remain FUBAR'ed anyway.  DO NOT USE this
// version of anfo2sam unless you absolutely must.

struct Fixup {
	const char *sequence ;
	uint32_t start_of_fix ;
	int32_t offset ;
} ;

Fixup fixups[] = {
	"chr3",			0,			-2,
	"chr5",			0,			-2,
	"chr6",			0,			-2,
	"chr7",			0,			-2,
	"chr10",		0,			-2,
	"chr11",		0,			-2,
	"chr12",		0,			-2,
	"chr12",		122494037,	-4,
	"chr13",		0,			-2,
	"chr14",		0,			-2,
	"chr15",		0,			-2,
	"chr19",		0,			-2,
	"chr20",		0,			-2,
	"chr21",		0,			-2,
	"chr21",		32078910,	-4,
	"chr21",		32079256,	-6,
	"chr21",		39207822,	-8,
	"chr21",		42908973,	-10,
	"chr22",		0,			-2,
	"chr22",		17558165,	-4,
	"chr1_random",	483265,		-2,
	"chr7_random",	402174,		-2,
	"chr9_random",	92729,		-2,
	"chr9_random",	144828,		-4
} ;

uint32_t fix_position( const std::string& seq, uint32_t raw_pos )
{
	int32_t off = 0 ;
	for( size_t i = 0 ; i != sizeof(fixups)/sizeof(fixups[0]) ; ++i )
		if( seq == fixups[i].sequence && raw_pos >= fixups[i].start_of_fix )
			off = fixups[i].offset ;
	return raw_pos+off ;
}

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
unsigned len_from_bin_cigar( const std::string& cigar ) {
	unsigned l = 0 ;
	for( size_t i = 0 ; i != cigar.size() ; ++i )
	{
		if( (uint8_t)cigar[i] < 128 ) l += (uint8_t)cigar[i] ;
		else if( (uint8_t)cigar[i] < 192 ) l += (unsigned)(uint8_t)cigar[i] - 128 ;
	}
	return l ;
}

bool is_clean_dna( const std::string& s )
{
	for( size_t i = 0 ; i != s.size() ; ++i )
	{
		switch( s[i] & ~32 )
		{
			case 'A': case 'C': case 'G': case 'T': case 'U': case 'N': break ;
			default: return false ;
		}
	}
	return true ;
}

std::string revcompl( std::string s )
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
			case 'U':
			case 'T': s[i] = 'A' ; break ;
		}
	}
	return s ;
}

enum bad_stuff { goodness = 0, no_hit, multiple_hits, no_seqid, no_seq, bad_cigar, parse_accident, bad_stuff_max } ; 
const char *descr[] = { 0, "had no hit", "had multiple hits", "missed the sequence id"
	                  , "missed the sequence", "had a bad CIGAR", "were improperly parsed" } ;

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
	if (!is_clean_dna( result.sequence() )) return parse_accident ;

	// XXX 
	// Either one of those two is fine, depending on what we actually
	// want.  mixing them is probably wrong.  (But right now it's okay,
	// since only best_to_genome will ever be present.)
    output::Hit hit = (result.has_best_hit())?result.best_hit():result.best_to_genome();

	if (len_from_bin_cigar(hit.cigar()) != result.sequence().length()) return bad_cigar ;

/*QNAME*/   std::cout << result.seqid() << "\t";
/*FLAG */   std::cout << (hit.aln_length() < 0 ? BAM_FREVERSE : 0) << "\t";  // TODO: calc flag
/*RNAME*/   std::cout << hit.sequence() << "\t";
	        uint32_t eff_start = fix_position( hit.sequence(), hit.old_start_pos() ) ;
/*POS*/     std::cout << ( hit.aln_length() >= 0 ? eff_start : eff_start+hit.aln_length()+1 ) << "\t";
/*MAPQ*/   	std::cout << ( result.has_diff_to_next() 
					? (int)( 0.5 + result.diff_to_next() / std::log(10.0) ) : 255 ) << '\t' ;
/*CIGAR*/   decode_binCigar(std::cout, hit.cigar()) << "\t";
          // We don't have paired end reads
/*MRNM*/    std::cout << "*" << "\t";
/*MPOS*/    std::cout << "0" << "\t";
/*ISIZE*/   std::cout << "0" << "\t";

/*SEQ*/     if( hit.aln_length() >= 0 ) std::cout << result.sequence() << '\t' ;
            else                        std::cout << revcompl( result.sequence() ) << "\t";
/*QUAL*/    if( result.has_quality() )
			{
				if( hit.aln_length() >= 0 )
					for (size_t i = 0; i != result.quality().size(); ++i)
						std::cout << char((uint8_t)result.quality()[i] + 33);
				else
					for (size_t i = result.quality().size(); i != 0 ; --i)
						std::cout << char((uint8_t)result.quality()[i-1] + 33);
			}
			else
				std::cout << '*' ;

/*[TAGS]*/
/*SCORE*/	std::cout << "\tAS:i:" << hit.score() ;

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

	std::cerr << "WARNING: This program contains an ugly hack.\n"
		         "DO NOT USE this unless you know what you're doing!" << std::endl ;
    int discarded[bad_stuff_max] = {0};
    AnfoFile f_anfo(argv[1]);
    f_anfo.get_header();

    for( output::Result res ; f_anfo.read_result( res ) ; ) 
        if (bad_stuff r = protoHit_2_bam_Hit(res))
            discarded[r]++;

	for( int b = 1 ; b != bad_stuff_max ; ++b )
		if (discarded[b])
			std::cerr << discarded[b] << " reads " << descr[b] << std::endl;

    return (EXIT_SUCCESS);
}

