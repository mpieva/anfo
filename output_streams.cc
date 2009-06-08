#include "output_streams.h"

#include "conffile.h"
#include "output.pb.h"

#include <google/protobuf/text_format.h>
#include <cmath>

namespace streams {

void TextWriter::print_msg( const google::protobuf::Message& m )
{
	google::protobuf::TextFormat::Print( m, &fos_ ) ;
	void* data ; int size ;
	fos_.Next( &data, &size ) ;
	*(char*)data = '\n' ;
	data = (char*)data + 1 ;
	--size ;
	if( !size ) fos_.Next( &data, &size ) ;
	*(char*)data = '\n' ;
	data = (char*)data + 1 ;
	--size ;
	if( size ) fos_.BackUp( size ) ;
}

void TextWriter::add_alignment( std::string::const_iterator qry, output::Hit &h, const config::Config &conf ) 
{
	CompactGenome &g = genomes_[ h.genome_file() ] ;
	if( !g.get_base() ) 
	{
		try { CompactGenome( h.genome_file(), conf, MADV_RANDOM ).swap( g ) ; }
		catch(...) { /* too bad, we'll make do without a genome */ }
	}

	if( g.get_base() ) 
	{
		std::string &r = *h.mutable_ref(), &q = *h.mutable_qry(), &c = *h.mutable_con() ;
		r.clear() ; q.clear() ; c.clear() ;
		DnaP ref = g.find_pos( h.sequence(), h.start_pos() ) ; 
		if( h.aln_length() < 0 ) ref = ref.reverse() + h.aln_length() + 1 ;

		for( size_t i = 0 ; i != h.cigar().size() ; ++i )
		{
			if( (uint8_t)h.cigar()[i] == 0 ) {
				r.push_back('~') ;
				q.push_back('~') ;
				c.push_back('~') ;
			}
			else if( (uint8_t)h.cigar()[i] < 128 ) for( size_t j = 0 ; j != (uint8_t)h.cigar()[i] ; ++j ) {
				r.push_back( from_ambicode( *ref ) ) ;
				q.push_back( *qry ) ;
				c.push_back( from_ambicode( *ref ) == *qry ? '*' : ' ' ) ;
				++ref; ++qry ;
			}
			else if( (uint8_t)h.cigar()[i] < 192 ) for( size_t j = 128 ; j != (uint8_t)h.cigar()[i] ; ++j, ++qry ) {
				r.push_back( '-' ) ;
				q.push_back( *qry ) ;
				c.push_back( ' ' ) ;
			}
			else for( size_t j = 192 ; j != (uint8_t)h.cigar()[i] ; ++j, ++ref ) {
				r.push_back( from_ambicode( *ref ) ) ;
				q.push_back( '-' ) ;
				c.push_back( ' ' ) ;
			}
		}
	}
	else
	{
		std::string &q = *h.mutable_qry() ;
		q.clear() ;

		for( size_t i = 0 ; i != h.cigar().size() ; ++i )
		{
			if( (uint8_t)h.cigar()[i] == 0 ) {
				q.push_back('~') ;
			}
			else if( (uint8_t)h.cigar()[i] < 128 ) for( size_t j = 0 ; j != (uint8_t)h.cigar()[i] ; ++j, ++qry ) {
				q.push_back( *qry ) ;
			}
			else if( (uint8_t)h.cigar()[i] < 192 ) for( size_t j = 128 ; j != (uint8_t)h.cigar()[i] ; ++j, ++qry ) {
				q.push_back( *qry ) ;
			}
			else for( size_t j = 192 ; j != (uint8_t)h.cigar()[i] ; ++j ) {
				q.push_back( '-' ) ;
			}
		}
	}
}

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
//! About the hit to convert: ANFO is happy representing more than one
//! hit per sequence, but SAM is not (or is it?).  We will convert the
//! \c best_to_genome entry, ignoring everything else.  This is both the
//! most useful alignment and the one best in line with the intended use
//! of SAM.
//!
//! About len: since we're doing semi-global alignments, it is actually
//! an error if the sequence length (taking trim points into account)
//! and the effective len from the cigar do not match.  However, files
//! with such errors exist, so we need to filter out those faulty
//! alignments.
//!
//! About the alignment score: This optional field of SAM is filled with
//! the raw ANFO alignment score, which is a bit difficult to interpret
//! if you're used to other aligners.  A lower score means a better
//! alignment, and a score of (length-20) * 18 is about the quality that
//! would give a zero score in BLAST.  Since SAM doesn't specify
//! anything about the alignment, we'll leave it that way.
//!
//! \todo calculate other SAM/BAM flags

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

SamWriter::bad_stuff SamWriter::protoHit_2_bam_Hit( const output::Result &result )
{

    if (!result.has_best_to_genome()) return no_hit ;
    if (!result.has_seqid()) return no_seqid ;
    if (!result.has_sequence()) return no_seq ;

    output::Hit hit = result.best_to_genome() ;

	if (streams::len_from_bin_cigar(hit.cigar()) != result.sequence().length()) return bad_cigar ;

	int mapq = !result.has_diff_to_next() ? 254 : std::min( 254,
			(int)( 0.5 + result.diff_to_next() / std::log(10.0) ) ) ;

	out_ << /*QNAME*/   result.seqid() << '\t'
         << /*FLAG */ ( hit.aln_length() < 0 ? bam_freverse : 0 ) << '\t'
         << /*RNAME*/   hit.sequence() << '\t'
         << /*POS*/     1 + hit.start_pos() << '\t'
		 << /*MAPQ*/    mapq << '\t' ;

	/*CIGAR*/ 
	if( hit.aln_length() >= 0 ) decode_binCigar( out_, hit.cigar().begin(),  hit.cigar().end() ) ;
	else                        decode_binCigar( out_, hit.cigar().rbegin(), hit.cigar().rend() ) ;
	
	out_ << '\t'
		 // We don't have paired end reads (or don't deal with them)
		 << /*MRNM*/    "*" << '\t'
		 << /*MPOS*/    "0" << '\t'
		 << /*ISIZE*/   "0" << '\t'

         << /*SEQ*/     result.sequence() << '\t' ;

	if( result.has_quality() ) /*QUAL*/   
		for (size_t i = 0; i < result.quality().size(); i++)
			out_ << char((uint8_t)result.quality()[i] + 33);
	else
		out_ << '*' ;

	/*[TAGS]*/ /*SCORE*/
	out_ << "\tAS:i:" << hit.score() << '\n' ;

    return goodness ;
}

const char *SamWriter::descr[] = { 0, "had no hit", "had multiple hits", "missed the sequence id"
	                             , "missed the sequence", "had a bad CIGAR" } ;

void SamWriter::put_footer( const Footer& f )
{
	for( int b = 1 ; b != bad_stuff_max ; ++b )
		if (discarded[b])
			std::clog << "SamWriter: " << discarded[b] << " reads " << descr[b] << std::endl;
	state_ = end_of_stream ;
}

} // namespace

