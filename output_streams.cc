#include "output_streams.h"

#include "conffile.h"
#include "output.pb.h"

#include <google/protobuf/text_format.h>
#include <google/protobuf/io/printer.h>
#include <cmath>
#include <iterator>
#include <sstream>

namespace {

void show_alignment(
		std::string::const_iterator qry,
		const output::Hit &h,
		const config::Config &conf,
		Genomes &genomes,
		bool for_fastq,
		std::string& r,
		std::string& q,
		std::string& c ) 
{
	CompactGenome &g = genomes[ h.genome_file() ] ;
	if( !g.get_base() ) 
	{
		// if creating FASTQ, we need the genome, else we can make do
		// without it
		try { CompactGenome( h.genome_file(), conf, MADV_RANDOM ).swap( g ) ; }
		catch(...) { if( for_fastq ) throw ; }
	}

	if( g.get_base() ) 
	{
		r.clear() ; q.clear() ; c.clear() ;
		DnaP ref = g.find_pos( h.sequence(), h.start_pos() ) ; 
		if( h.aln_length() < 0 ) ref = ref.reverse() + h.aln_length() + 1 ;

		for( int i = 0 ; i != h.cigar_size() ; ++i )
		{
			unsigned l = h.cigar(i) >> 3 ;
			switch( h.cigar(i) & 7 )
			{
				case output::Hit::Match:
					if( !l && !for_fastq ) {
						r.push_back('~') ;
						q.push_back('~') ;
						c.push_back('~') ;
					}

				case output::Hit::Mismatch:
					for( size_t j = 0 ; j != l ; ++j, ++ref, ++qry ) {
						r.push_back( from_ambicode( *ref ) ) ;
						q.push_back( *qry ) ;
						c.push_back( from_ambicode( *ref ) == *qry ? '*' : ' ' ) ;
					}
					break ;

				case output::Hit::Insert:
					for( size_t j = 0 ; j != l ; ++j, ++qry ) {
						r.push_back( '-' ) ;
						q.push_back( *qry ) ;
						c.push_back( ' ' ) ;
					}
					break ;

				case output::Hit::Delete:
					for( size_t j = 0 ; j != l ; ++j, ++ref ) {
						r.push_back( from_ambicode( *ref ) ) ;
						q.push_back( '-' ) ;
						c.push_back( ' ' ) ;
					}
					break ;
			}
		}
	}
	else
	{
		r.clear() ; q.clear() ; c.clear() ;

		for( int i = 0 ; i != h.cigar().size() ; ++i )
		{
			unsigned l = h.cigar(i) >> 3 ;
			switch( h.cigar(i) & 7 )
			{
				case output::Hit::Match:
					if( !l && !for_fastq ) q.push_back('~') ;

				case output::Hit::Mismatch:
				case output::Hit::Insert:
					for( size_t j = 0 ; j != l ; ++j, ++qry ) q.push_back( *qry ) ;
					break ;

				case output::Hit::Delete:
					for( size_t j = 0 ; j != l ; ++j ) q.push_back( '-' ) ;
					break ;
			}
		}
	}
}

} // namespace

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

void TextWriter::put_result( const Result& r )
{
	google::protobuf::TextFormat::Print( r, &fos_ ) ;
	google::protobuf::io::Printer p( &fos_, '`' ) ;
	for( int i = 0 ; i != r.hit_size() ; ++i )
	{
		std::map< std::string, std::string > vars ;
		show_alignment( r.read().sequence().begin(), r.hit(i), hdr_.config(), genomes_,
				false, vars["ref"], vars["qry"], vars["con"] ) ;
		p.Print( vars, "\nREF: `ref`\nQRY: `qry`\nCON: `con`\n" ) ;
	}
	p.Print( "\n\n" ) ;
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
std::ostream& decode_binCigar(std::ostream& s, Iter begin, Iter end )
{
	unsigned len = 0, op = Hit::Match ;

	for( ; begin != end ; ++begin )
	{
		unsigned o = *begin & 7, l = *begin >> 3  ;
		if( o == Hit::Mismatch ) o = Hit::Match ;
		if( o == op ) len += l ;
		else {
			if( len ) s << len << "MID"[op] ;
			op = o ;
			len = l ;
		}
	}
	if( len ) s << len << "MID"[op] ;
    return s ;
}

template <typename Iter> std::reverse_iterator<Iter> mk_rev_iter( Iter i )
{ return std::reverse_iterator<Iter>( i ) ; }

SamWriter::bad_stuff SamWriter::protoHit_2_bam_Hit( const output::Result &result )
{
    if (!has_hit_to(result,g_)) return no_hit ;
    if (!result.read().has_seqid()) return no_seqid ;
    if (!result.read().has_sequence()) return no_seq ;

    output::Hit hit = hit_to(result,g_) ;

	if (len_from_bin_cigar(hit.cigar()) != result.read().sequence().length()) return bad_cigar ;

	int mapq = !hit.has_diff_to_next() ? 254 : std::min( 254, hit.diff_to_next() ) ;

	out_ << /*QNAME*/   result.read().seqid() << '\t'
         << /*FLAG */ ( hit.aln_length() < 0 ? bam_freverse : 0 ) << '\t'
         << /*RNAME*/   hit.sequence() << '\t'
         << /*POS*/     1 + hit.start_pos() << '\t'
		 << /*MAPQ*/    mapq << '\t' ;

	/*CIGAR*/ 
	if( hit.aln_length() >= 0 ) decode_binCigar( out_, hit.cigar().begin(),  hit.cigar().end() ) ;
	else                        decode_binCigar( out_, mk_rev_iter( hit.cigar().end() ), mk_rev_iter( hit.cigar().begin() ) ) ;
	
	out_ << '\t'
		 // We don't have paired end reads (or don't deal with them)
		 << /*MRNM*/    "*" << '\t'
		 << /*MPOS*/    "0" << '\t'
		 << /*ISIZE*/   "0" << '\t'

         << /*SEQ*/     result.read().sequence() << '\t' ;

	if( result.read().has_quality() ) /*QUAL*/   
		for (size_t i = 0; i < result.read().quality().size(); i++)
			out_ << char((uint8_t)result.read().quality()[i] + 33);
	else
		out_ << '*' ;

	/*[TAGS]*/ /*SCORE*/
	out_ << "\tAS:i:" << hit.score() << '\n' ;

    return goodness ;
}

const char *SamWriter::descr[] = { "were converted", "had no hit", "had multiple hits", "missed the sequence id"
	                             , "missed the sequence", "had a bad CIGAR" } ;

void SamWriter::put_footer( const Footer& f )
{
	for( int b = 0 ; b != bad_stuff_max ; ++b )
		if (discarded[b]) {
			std::stringstream s ;
			s << "SamWriter: " << discarded[b] << " reads " << descr[b] ;
			console.output( Console::notice, s.str() ) ;
		}
	state_ = end_of_stream ;
}

void FastaWriter::put_result( const Result& r ) 
{
	if( has_hit_to( r, g_ ) )
	{
		const Hit &h = hit_to( r, g_ ) ;
		std::string ref, qry, con ;
		show_alignment( r.read().sequence().begin(), h, hdr_.config(), genomes_, true, ref, qry, con ) ;
		out_ << '>' << h.sequence() << ' '
			<< h.start_pos()
			<< "-+"[ h.aln_length() > 0 ]
			<< h.start_pos() + abs(h.aln_length()) - 1
			<< '\n' << ref << '\n' 
			<< '>' << r.read().seqid() 
			<< ( r.read().has_trim_right() ? " adapter cut off\n" : "\n" ) 
			<< qry << std::endl ;
	}
}

void TableWriter::put_result( const Result& r )
{
	if( !has_hit_to( r, g_ ) ) return ;
	int e = r.read().has_trim_right() ? r.read().trim_right() : r.read().sequence().size() ;
	int b = r.read().has_trim_left() ? r.read().trim_left() : 0 ;
	int diff = hit_to( r, g_ ).has_diff_to_next() ? hit_to( r, g_ ).diff_to_next() : 9999 ;

	out_ << e-b << '\t' << r.hit(0).score() << '\t' << diff << '\n' ;
}

} // namespace

