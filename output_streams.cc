#include "config.h"
#include "output_streams.h"

#include "conffile.h"
#include "output.pb.h"

#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/printer.h>
#include <google/protobuf/text_format.h>
#include <cmath>
#include <iterator>
#include <sstream>

namespace streams {
namespace {

void show_alignment(
		std::string::const_iterator qry,
		const output::Hit &h,
		bool for_fasta,
		std::string& r,
		std::string& q,
		std::string& c,
		int context = 0 ) 
{
	CompactGenome *g = 0 ;
	DnaP ref( 0 ) ;

	try {
		g = &Metagenome::find_sequence( h.genome_name(), h.sequence(), Metagenome::ephemeral ) ;
		if( g ) ref = g->find_pos( h.sequence(), h.start_pos() ) ; 
		if( ref && h.aln_length() < 0 ) ref = ref.reverse() + h.aln_length() + 1 ;
		if( !ref ) throw "sequence not found: " + h.sequence() ;
	}
	// if creating FASTA, we need the genome, else we can make do without it
	catch(...) { if( for_fasta ) throw ; }

	if( ref )
	{
		r.clear() ; q.clear() ; c.clear() ;
		if( context ) 
		{
			DnaP ref1 = ref ;
			for( int i = 0 ; i != context ; ++i )
			{
				r = from_ambicode( ref1[-1] ) + r ;
				if( ref1[-1] ) --ref1 ;
			}
			q = std::string( context, '-' ) ;
			c = std::string( context, ' ' ) ;
		}

		for( int i = 0 ; i != h.cigar_size() ; ++i )
		{
			unsigned l = cigar_len( h.cigar(i) ) ;
			switch( cigar_op( h.cigar(i) ) )
			{
				case output::Hit::Match:
					if( !l && !for_fasta ) {
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

				case output::Hit::SoftClip: break ; // ???
				case output::Hit::Skip: break ; // ???
				case output::Hit::HardClip: break ; // ???
				case output::Hit::Pad: break ; // ???
			}
		}

		for( int i = 0 ; i != context ; ++i )
		{
			r.push_back( from_ambicode( *ref ) ) ;
			if( *ref ) ++ref ;
		}
		q += std::string( context, '-' ) ;
		c += std::string( context, ' ' ) ;
	}
	else
	{
		r.clear() ; q.clear() ; c.clear() ;

		for( int i = 0 ; i != h.cigar().size() ; ++i )
		{
			unsigned l = cigar_len( h.cigar(i) ) ;
			switch( cigar_op( h.cigar(i) ) )
			{
				case output::Hit::Match:
					if( !l && !for_fasta ) {
						r.push_back('~') ;
						q.push_back('~') ;
						c.push_back('~') ;
					}
					else for( size_t j = 0 ; j != l ; ++j, ++qry ) {
						r.push_back('N') ;
						q.push_back( *qry ) ;
						c.push_back('*') ;
					}
					break ;

				case output::Hit::Mismatch:
					for( size_t j = 0 ; j != l ; ++j, ++qry ) {
						r.push_back('N') ;
						q.push_back( *qry ) ;
						c.push_back(' ') ;
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
					for( size_t j = 0 ; j != l ; ++j ) {
						r.push_back( 'N' ) ;
						q.push_back( '-' ) ;
						c.push_back( ' ' ) ;
					}
					break ;

				case output::Hit::SoftClip: break ; // ???
				case output::Hit::Skip: break ; // ???
				case output::Hit::HardClip: break ; // ???
				case output::Hit::Pad: break ; // ???
			}
		}
	}
}

} // namespace

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
		show_alignment( r.read().sequence().begin(), r.hit(i),
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

	for( Iter cur = begin ; cur != end ; )
	{
		unsigned o = cigar_op( *begin ), l = cigar_len( *begin ) ;
		if( o == Hit::Mismatch ) o = Hit::Match ;
		if( o == Hit::Insert && begin == cur ) o = Hit::SoftClip ;
		++cur ;
		if( o == Hit::Insert && end == cur ) o = Hit::SoftClip ;

		if( o == op ) len += l ;
		else {
			if( len ) s << len << "MIDNSHP"[op] ;
			op = o ;
			len = l ;
		}
	}
	if( len ) s << len << "MIDNSHP"[op] ;
    return s ;
}

template <typename Iter> std::reverse_iterator<Iter> mk_rev_iter( Iter i )
{ return std::reverse_iterator<Iter>( i ) ; }

inline char qual_to_sam( uint8_t q ) { return 33 + std::min( q, (uint8_t)93 ) ; }

SamWriter::bad_stuff SamWriter::protoHit_2_bam_Hit( const output::Result &result )
{
	if (!has_hit_to(result,g_)) return no_hit ;
	if (!result.read().has_seqid()) return no_seqid ;
	if (!result.read().has_sequence()) return no_seq ;

	output::Hit hit = hit_to(result,g_) ;

	if( len_from_bin_cigar( hit.cigar() ) != result.read().sequence().length() ) return bad_cigar ;

	int mapq = !hit.has_diff_to_next() ? 254 : std::min( 254, hit.diff_to_next() ) ;

	out_ << /*QNAME*/  result.read().seqid() << '\t'
		<< /*FLAG */ ( hit.aln_length() < 0 ? bam_freverse : 0 ) << '\t'
		<< /*RNAME*/   hit.sequence() << '\t'
		<< /*POS*/     1 + hit.start_pos() << '\t'
		<< /*MAPQ*/    mapq << '\t' ;

	if( hit.aln_length() >= 0 )
	{
		decode_binCigar( out_, hit.cigar().begin(),  hit.cigar().end() ) ; /*CIGAR*/ 
		// We don't have paired end reads (or don't deal with them)
		out_ << "\t*\t0\t0\t" // MRNM, MPOS, ISIZE
			<< /*SEQ*/ result.read().sequence() << '\t' ;

		if( result.read().has_quality() ) /*QUAL*/   
			for (size_t i = 0 ; i != result.read().quality().size() ; ++i )
				out_ << qual_to_sam( result.read().quality()[i] ) ;
		else
			out_ << '*' ;
	}
	else
	{
		// need to revcom sequence, reverse qual and cigar
		decode_binCigar( out_, hit.cigar().begin(),  hit.cigar().end() ) ; /*CIGAR*/ 
		// We don't have paired end reads (or don't deal with them)
		out_ << "\t*\t0\t0\t" ; // MRNM, MPOS, ISIZE
		const std::string& s = result.read().sequence() ;
		for( size_t i = s.size() ; i != 0 ; --i )
			switch( s[i-1] )
			{
				case 'A': case 'a': out_ << 'T' ; break ;
				case 'C': case 'c': out_ << 'G' ; break ;
				case 'G': case 'g': out_ << 'C' ; break ;
				case 'T': case 't':
				case 'U': case 'u': out_ << 'A' ; break ;
				default: out_ << s[i-1] ;
			}

		out_ << '\t' ;
		if( result.read().has_quality() ) /*QUAL*/   
			for (size_t i = result.read().quality().size() ; i != 0 ; --i )
				out_ << qual_to_sam( result.read().quality()[i-1] ) ;
		else
			out_ << '*' ;
	}

	/*[TAGS]*/ /*SCORE*/
	out_ << "\tAS:i:" << hit.score() << '\n' ;

	return goodness ;
}

const char *SamWriter::descr[] = { "were converted", "had no hit", "had multiple hits", "missed the sequence id"
	                             , "missed the sequence", "had a bad CIGAR" } ;

void SamWriter::put_footer( const Footer& f )
{
	Stream::put_footer( f ) ;
	for( int b = 0 ; b != bad_stuff_max ; ++b )
		if (discarded[b]) {
			std::stringstream s ;
			s << "SamWriter: " << discarded[b] << " reads " << descr[b] ;
			console.output( Console::notice, s.str() ) ;
		}
}

void FastaAlnWriter::put_result( const Result& r ) 
{
	if( has_hit_to( r, g_ ) )
	{
		const Hit &h = hit_to( r, g_ ) ;
		std::string ref, qry, con ;
		show_alignment( r.read().sequence().begin(), h, true, ref, qry, con, c_ ) ;
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

void FastqWriter::put_result( const Result& rr ) 
{
    const output::Read& r = rr.read() ;
	if( r.has_quality() ) 
	{
		out_ << '@' << r.seqid() ;
		if( r.has_description() ) out_ << ' ' << r.description() ;
		for( size_t i = 0 ; i < r.sequence().size() ; i += 50 )
			out_ << '\n' << r.sequence().substr( i, 50 ) ;
		out_ << "\n+\n" ;
        const std::string& q = r.quality() ;
		for( size_t i = 0 ; i < q.size() ; i += 50 )
		{
			for( size_t j = i ; j != i+50 && j != q.size() ; ++j )
				out_ << (char)std::min(126, 33 + q[j]) ;
			out_ << '\n' ;
		}
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

