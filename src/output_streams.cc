//    Copyright 2009 Udo Stenzel
//    This file is part of ANFO
//
//    ANFO is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    Anfo is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with Anfo.  If not, see <http://www.gnu.org/licenses/>.

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include "output_streams.h"

#include "conffile.h"
#include "output.pb.h"

#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/printer.h>
#include <google/protobuf/text_format.h>
#include <cmath>
#include <iterator>
#include <sstream>

#include <dlfcn.h>

namespace streams {

using namespace google::protobuf ;

bool GenTextAlignment::xform( Result& res ) 
{
	for( int i = 0 ; i != res.hit_size() ; ++i )
	{
		output::Hit& h = *res.mutable_hit(i) ;
		std::string::const_iterator qry =
			res.read().sequence().begin() + res.read().trim_left() ;

		std::string& r = *h.mutable_aln_ref() ;
		std::string& q = *h.mutable_aln_qry() ;

		GenomeHolder g = Metagenome::find_sequence( h.genome_name(), h.sequence() ) ;
		DnaP ref = g->find_pos( h.sequence(), h.start_pos() ) ; 
		if( h.aln_length() < 0 ) ref = ref.reverse() + h.aln_length() + 1 ;

		r.clear() ; q.clear() ;
		if( context_ ) 
		{
			DnaP ref1 = ref ;
			for( int i = 0 ; i != context_ ; ++i )
			{
				r = from_ambicode( ref1[-1] ) + r ;
				if( ref1[-1] ) --ref1 ;
			}
			q = std::string( context_, '-' ) ;
		}

		for( int i = 0 ; i != h.cigar_size() ; ++i )
		{
			unsigned l = cigar_len( h.cigar(i) ) ;
			switch( cigar_op( h.cigar(i) ) )
			{
				case output::Hit::Match:
				case output::Hit::Mismatch:
					for( size_t j = 0 ; j != l ; ++j, ++ref, ++qry ) {
						r.push_back( from_ambicode( *ref ) ) ;
						q.push_back( *qry ) ;
					}
					break ;

				case output::Hit::SoftClip:
				case output::Hit::Insert:
					for( size_t j = 0 ; j != l ; ++j, ++qry ) {
						r.push_back( '-' ) ;
						q.push_back( *qry ) ;
					}
					break ;

				case output::Hit::Delete:
					for( size_t j = 0 ; j != l ; ++j, ++ref ) {
						r.push_back( from_ambicode( *ref ) ) ;
						q.push_back( '-' ) ;
					}
					break ;

				case output::Hit::Skip: break ; // ???
				case output::Hit::HardClip: break ; // ???
				case output::Hit::Pad: break ; // ???
			}
		}

		for( int i = 0 ; i != context_ ; ++i )
		{
			r.push_back( from_ambicode( *ref ) ) ;
			if( *ref ) ++ref ;
		}
		q += std::string( context_, '-' ) ;
	}
	return true ;
}

typedef void (*PrintMethod)( const Message& m, const Reflection*, const FieldDescriptor*, ostream&, int ) ;

void generic_show_message( const Message&, std::ostream&, int ) ;
void generic_show_known_fields( const Message&, std::ostream&, int ) ;
void generic_show_unknown_fields( const UnknownFieldSet&, std::ostream&, int ) ;

// hides a field by printing nothing
extern "C" void show_noop( const Message& m, const Reflection*, const FieldDescriptor*, ostream& s, int ) {}

// CIGAR line: decode to string, same way SAM would do it
extern "C" void show_cigar( const Message& m, const Reflection* r, const FieldDescriptor* f, ostream& s, int indent ) 
{
	s << string( indent, ' ' ) << f->name() << ": " ;
	for( int i = 0 ; i != r->FieldSize( m, f ) ; ++i )
	{
		uint32_t c = r->GetRepeatedUInt32( m, f, i ) ;
		if( c == 0 ) s << '/' ;
		else s << cigar_len(c) << "MIDNSHP........M"[ cigar_op(c) ] ;
	}
	s << '\n' ;
}

// array of ints: formatted into comma-separated line
extern "C" void show_uint32_array( const Message& m, const Reflection* r, const FieldDescriptor* f, ostream& s, int indent ) 
{
	s << string( indent, ' ' ) << f->name() << ": " ;
	if( int l = r->FieldSize( m, f ) ) 
	{
		s << r->GetRepeatedUInt32( m, f, 0 ) ;
		for( int i = 1 ; i != l ; ++i )
			s << ',' << r->GetRepeatedUInt32( m, f, i ) ;
	}
	s << '\n' ;
}
extern "C" void show_int32_array( const Message& m, const Reflection* r, const FieldDescriptor* f, ostream& s, int indent ) 
{
	s << string( indent, ' ' ) << f->name() << ": " ;
	if( int l = r->FieldSize( m, f ) ) 
	{
		s << r->GetRepeatedInt32( m, f, 0 ) ;
		for( int i = 1 ; i != l ; ++i )
			s << ',' << r->GetRepeatedInt32( m, f, i ) ;
	}
	s << '\n' ;
}

// seen bases: same as array of ints, but we have four interleaved
// arrays
extern "C" void show_seen_bases( const Message& m, const Reflection* r, const FieldDescriptor* f, ostream& s, int indent ) 
{
	for( int i = 0 ; i != 4 ; ++i )
	{
		s << string( indent, ' ' ) << f->name() << '_' << "ACGT"[i] << ": " ;
		int l = r->FieldSize( m, f ) ;
		if( i <= l )
		{
			s << r->GetRepeatedUInt32( m, f, i ) ;
			for( int j = i+4 ; j < l ; j+=4 )
				s << ',' << r->GetRepeatedUInt32( m, f, j ) ;
		}
		s << '\n' ;
	}
}

// Hit: basically a generic message, but we add the textual,
// line-wrapped alignment (if present)
extern "C" void show_hit( const Message& m, const Reflection* r, const FieldDescriptor* f, ostream& s, int i ) 
{
	for( int j = 0 ; j != r->FieldSize( m, f ) ; ++j )
	{
		const output::Hit& h = static_cast<const output::Hit&>( r->GetRepeatedMessage( m, f, j ) ) ;
		s << string( i, ' ' ) << f->name() << " {\n" ;
		generic_show_known_fields( h, s, i+2 ) ;
		
		// prints a hit: REF & QRY are hidden anyway; here we add them 
		if( h.has_aln_ref() && h.has_aln_qry() && h.aln_ref().size() == h.aln_qry().size() ) {
			const std::string& r = h.aln_ref(), q = h.aln_qry() ;
			for( size_t k = 0 ; k < r.size() ; k += 50 )
			{
				s << '\n' ;
				s << string( i+2, ' ' ) << "REF: " << r.substr( k, 50 ) << '\n' ;
				s << string( i+2, ' ' ) << "QRY: " << q.substr( k, 50 ) << '\n' ;
				s << string( i+2, ' ' ) << "CNS: " ;
				for( size_t l = k ; l != k+50 && l != r.size() ; ++l )
					s << " *"[ r[l] == q[l] ] ;
				s << '\n' ;
			}
		}

		generic_show_unknown_fields( r->GetUnknownFields(m), s, i+2 ) ;
		s << string( i, ' ' ) << "}\n" ;
	}
}

void generic_show_string( const string& s, std::ostream& o, int off = 0 )
{
	o << '"' ;
	for( size_t i = 0 ; i != s.size() ; ++i )
	{
		unsigned c = (unsigned char)s[i] + off ;
		switch( c )
		{
			case '\r': o << "\\r" ; break ;
			case '\n': o << "\\n" ; break ;
			case '"':  o << "\\\"" ; break ;
			case '\\': o << "\\\\" ; break ;
			case '\t': o << "\\t" ; break ;
			default:
					   if( c >= ' ' && c <= '~' && c != '"' ) o << (char)c ;
					   else o << '\\' << (char)( '0' + ((c >> 6) & 7) )
                                      << (char)( '0' + ((c >> 3) & 7) )
                                      << (char)( '0' + ((c >> 0) & 7) ) ;
		}
	}
	o << '"' ;
}

extern "C" void show_quality( const Message& m, const Reflection* r, const FieldDescriptor* f, ostream& s, int i ) 
{
	s << string( i, ' ' ) << f->name() << ' ' ;
	generic_show_string( r->GetString( m, f ), s, 33 ) ;
	s << '\n' ;
}

void universal_print_field( const Message& m, const Reflection* r, const FieldDescriptor* f, ostream& s, int indent )
{
	if( f->options().HasExtension( output::show ) )
	{
		const std::string& method_name = f->options().GetExtension( output::show ) ;
		if( void* p = dlsym( RTLD_DEFAULT, method_name.c_str() ) ) return ((PrintMethod)p)( m, r, f, s, indent ) ;
		std::stringstream ss ;
		ss << "TextWriter: " << method_name << " declared for "
			<< m.GetDescriptor()->name() << "::" << f->name() << ", but not found." ;
		console.output( Console::warning, ss.str() ) ;
	}

	if( f->label() == FieldDescriptor::LABEL_REPEATED ) {
		for( int i = 0 ; i != r->FieldSize( m, f ) ; ++i ) 
		{
			s << string( indent, ' ' ) << f->name() ;
			if( f->cpp_type() != FieldDescriptor::CPPTYPE_MESSAGE ) s << ':' ;
			s << ' ' ;

			switch( f->cpp_type() )
			{
				case FieldDescriptor::CPPTYPE_INT32:   s << r->GetRepeatedInt32( m, f, i ) ; break ;
				case FieldDescriptor::CPPTYPE_INT64:   s << r->GetRepeatedInt64( m, f, i ) ; break ;
				case FieldDescriptor::CPPTYPE_UINT32:  s << r->GetRepeatedUInt32( m, f, i ) ; break ;
				case FieldDescriptor::CPPTYPE_UINT64:  s << r->GetRepeatedUInt64( m, f, i ) ; break ;
				case FieldDescriptor::CPPTYPE_DOUBLE:  s << r->GetRepeatedDouble( m, f, i ) ; break ;
				case FieldDescriptor::CPPTYPE_FLOAT:   s << r->GetRepeatedFloat( m, f, i ) ; break ;
				case FieldDescriptor::CPPTYPE_BOOL:    s << ( r->GetRepeatedBool( m, f, i ) ? "true" : "false" ) ; break ;
				case FieldDescriptor::CPPTYPE_ENUM:    s << r->GetRepeatedEnum( m, f, i )->name() ; break ;
				case FieldDescriptor::CPPTYPE_STRING:  generic_show_string( r->GetRepeatedString( m, f, i ), s ) ; break ;
				case FieldDescriptor::CPPTYPE_MESSAGE: generic_show_message( r->GetRepeatedMessage( m, f, i ), s, indent ) ; break ;
				default: std::cerr << "don't know how to print cpp_type " << f->cpp_type() << "\n" ;
			}
			s << '\n' ;
		}
	}
	else 
	{
		s << string( indent, ' ' ) << f->name() ;
		if( f->cpp_type() != FieldDescriptor::CPPTYPE_MESSAGE ) s << ':' ;
		s << ' ' ;

		switch( f->cpp_type() )
		{
			case FieldDescriptor::CPPTYPE_INT32:   s << r->GetInt32( m, f ) ; break ;
			case FieldDescriptor::CPPTYPE_INT64:   s << r->GetInt64( m, f ) ; break ;
			case FieldDescriptor::CPPTYPE_UINT32:  s << r->GetUInt32( m, f ) ; break ;
			case FieldDescriptor::CPPTYPE_UINT64:  s << r->GetUInt64( m, f ) ; break ;
			case FieldDescriptor::CPPTYPE_DOUBLE:  s << r->GetDouble( m, f ) ; break ;
			case FieldDescriptor::CPPTYPE_FLOAT:   s << r->GetFloat( m, f ) ; break ;
			case FieldDescriptor::CPPTYPE_BOOL:    s << ( r->GetBool( m, f ) ? "true" : "false" ) ; break ;
			case FieldDescriptor::CPPTYPE_ENUM:    s << r->GetEnum( m, f )->name() ; break ;
			case FieldDescriptor::CPPTYPE_STRING:  generic_show_string( r->GetString( m, f ), s ) ; break ;
			case FieldDescriptor::CPPTYPE_MESSAGE: generic_show_message( r->GetMessage( m, f ), s, indent ) ; break ;
			default: std::cerr << "don't know how to print cpp_type " << f->cpp_type() << "\n" ;
		}
		s << '\n' ;
	}
}

void generic_show_known_fields( const Message& m, std::ostream& s, int indent )
{
	const Reflection* reflection = m.GetReflection() ;
	std::vector< const FieldDescriptor* > fields;
	reflection->ListFields( m, &fields ) ;
	for( size_t i = 0 ; i < fields.size() ; i++ )
	{
		universal_print_field( m, reflection, fields[i], s, indent ) ;
	}
}

void generic_show_unknown_fields( const UnknownFieldSet& fields, std::ostream& s, int indent )
{
	for( int i = 0 ; i != fields.field_count() ; ++i )
	{
		UnknownField f = fields.field( i ) ;
		s << string( indent, ' ' ) << f.number() ;
		if( f.type() != UnknownField::TYPE_GROUP ) s << ':' ;
		s << ' ' ;
		switch( f.type() )
		{
			case UnknownField::TYPE_VARINT: 			s << f.varint() ; break ;
			case UnknownField::TYPE_FIXED32: 			s << f.fixed32() ; break ;
			case UnknownField::TYPE_FIXED64: 			s << f.fixed64() ; break ;
			case UnknownField::TYPE_LENGTH_DELIMITED:	generic_show_string( f.length_delimited(), s ) ; break ;
			case UnknownField::TYPE_GROUP: 				s << "{\n" ;
														generic_show_unknown_fields( f.group(), s, indent+2 ) ;
														s << string( indent, ' ' ) << "}" ;
														break ;
			default: 									s << "???" ; break ;
		}
		s << '\n' ;
	}
}

void generic_show_message( const Message& m, std::ostream& s, int indent )
{
	s << "{\n" ;
	generic_show_known_fields( m, s, indent+2 ) ;
	generic_show_unknown_fields( m.GetReflection()->GetUnknownFields(m), s, indent+2 ) ;
	s << string( indent, ' ' ) << "}" ;
}

void show_message( const std::string& name, const Message& m, std::ostream& s )
{
	s << name << ' ' ;
	generic_show_message( m, s, 0 ) ;
	s << "\n\n\n" ;
}

void TextWriter::put_header( const Header& h )
{
	Stream::put_header( h ) ;
	show_message( "Header", h, *out_ ) ;
}

void TextWriter::put_result( const Result& r )
{
	show_message( "Result", r, *out_ ) ;
}

void TextWriter::put_footer( const Footer& f )
{
	Stream::put_footer( f ) ;
	show_message( "Footer", f, *out_ ) ;
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
//! hit per sequence, but SAM is not (or is it?).  We will convert
//! either the best hit globallly or the best one to a given genome.
//! This is both the most useful alignment and the one best in line with
//! the intended use of SAM.
//!
//! About len: since we're doing semi-global alignments, it is actually
//! an error if the sequence length (taking trim points into account)
//! and the effective len from the cigar do not match.  However, files
//! with such errors exist, so we need to filter out those faulty
//! alignments.
//!
//! About the CIGAR line: parts of the samtools suite fail if an
//! alignment declares an insertion that overhangs the reference.  To
//! deal with that, we declare inserts at the ends of alignments as soft
//! clips instead.
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
	std::vector< std::pair< unsigned, unsigned > > new_cigar ;

	unsigned len = 0, op = Hit::Match ;
	for( Iter cur = begin ; cur != end ; ++cur )
	{
		unsigned o = cigar_op( *cur ), l = cigar_len( *cur ) ;
		if( o == Hit::Mismatch ) o = Hit::Match ;

		if( o == op ) len += l ;
		else {
			if( len ) new_cigar.push_back( std::make_pair( len, op ) ) ; 
			op = o ;
			len = l ;
		}
	}
	if( len ) new_cigar.push_back( std::make_pair( len, op ) ) ;

	if( !new_cigar.empty() )
	{
		if( new_cigar.front().second == Hit::Insert )
			new_cigar.front().second = Hit::SoftClip ;

		if( new_cigar.back().second == Hit::Insert )
			new_cigar.back().second = Hit::SoftClip ;
	}
	for( size_t i = 0 ; i != new_cigar.size() ; ++i )
		s << new_cigar[i].first << "MIDNSHP"[ new_cigar[i].second ] ;

    return s ;
}

template <typename Iter> std::reverse_iterator<Iter> mk_rev_iter( Iter i )
{ return std::reverse_iterator<Iter>( i ) ; }

inline char qual_to_sam( uint8_t q ) { return 33 + std::min( q, (uint8_t)93 ) ; }

SamWriter::bad_stuff SamWriter::protoHit_2_bam_Hit( const output::Result &result )
{
	if (!result.read().has_seqid()) return no_seqid ;
	if (!result.read().has_sequence()) return no_seq ;

	for( int i = 0 ; i != result.hit_size() ; ++i )
	{
		const output::Hit &hit = result.hit(0) ;

		if( len_from_bin_cigar( hit.cigar() ) != result.read().sequence().length() ) return bad_cigar ;

		int mapq = !hit.has_diff_to_next() ? 254 : std::min( 254, hit.diff_to_next() ) ;

		*out_ << /*QNAME*/  result.read().seqid() << '\t'
			<< /*FLAG */ ( hit.aln_length() < 0 ? bam_freverse : 0 ) << '\t'
			<< /*RNAME*/   hit.sequence() << '\t'
			<< /*POS*/     1 + hit.start_pos() << '\t'
			<< /*MAPQ*/    mapq << '\t' ;

		if( hit.aln_length() >= 0 )
		{
			decode_binCigar( *out_, hit.cigar().begin(),  hit.cigar().end() ) ; /*CIGAR*/ 
			// We don't have paired end reads (or don't deal with them)
			*out_ << "\t*\t0\t0\t" // MRNM, MPOS, ISIZE
				<< /*SEQ*/ result.read().sequence() << '\t' ;

			if( result.read().has_quality() ) /*QUAL*/   
				for (size_t i = 0 ; i != result.read().quality().size() ; ++i )
					*out_ << qual_to_sam( result.read().quality()[i] ) ;
			else
				*out_ << '*' ;
		}
		else
		{
			// need to revcom sequence, reverse qual and cigar
			decode_binCigar( *out_, mk_rev_iter( hit.cigar().end() ), mk_rev_iter( hit.cigar().begin() ) ) ; /*CIGAR*/ 
			// We don't have paired end reads (or don't deal with them)
			*out_ << "\t*\t0\t0\t" ; // MRNM, MPOS, ISIZE
			const std::string& s = result.read().sequence() ;
			for( size_t i = s.size() ; i != 0 ; --i )
				switch( s[i-1] )
				{
					case 'A': case 'a': *out_ << 'T' ; break ;
					case 'C': case 'c': *out_ << 'G' ; break ;
					case 'G': case 'g': *out_ << 'C' ; break ;
					case 'T': case 't':
					case 'U': case 'u': *out_ << 'A' ; break ;
					default: *out_ << s[i-1] ;
				}

			*out_ << '\t' ;
			if( result.read().has_quality() ) /*QUAL*/   
				for (size_t i = result.read().quality().size() ; i != 0 ; --i )
					*out_ << qual_to_sam( result.read().quality()[i-1] ) ;
			else
				*out_ << '*' ;
		}

		/*[TAGS]*/ /* score is actually phred likelihood */
		if( hit.has_score() ) *out_ << "\tUQ:i:" << hit.score() ;
		*out_ << '\n' ;
	}

	if( result.hit_size() == 0 )
	{
		// special treatment of unaligned sequence
		*out_ << /*QNAME*/  result.read().seqid() << '\t'
			<< /*FLAG */ bam_funmap
			<< "\t*\t0\t0\t*\t*\t0\t0\t" // RNAME, POS, MAPQ, CIGAR, MRNM, MPOS, ISIZE
			<< /*SEQ*/ result.read().sequence() << '\t' ;

		if( result.read().has_quality() ) /*QUAL*/   
			for (size_t i = 0 ; i != result.read().quality().size() ; ++i )
				*out_ << qual_to_sam( result.read().quality()[i] ) ;
		else
			*out_ << '*' ;

		/*[TAGS]*/
		*out_ << "\n" ;
		return no_hit ;
	}
	else return goodness ;
}

const char *SamWriter::descr[] = { "were converted", "had no hit", "had multiple hits", "missed the sequence id"
	                             , "missed the sequence", "had a bad CIGAR" } ;

void SamWriter::put_footer( const Footer& f )
{
	Stream::put_footer( f ) ;
	for( int b = 0 ; b != bad_stuff_max ; ++b )
		if (discarded[b]) {
			std::stringstream s ;
			s << "SamWriter: " << nm_ << ": " << discarded[b] << " reads " << descr[b] ;
			console.output( Console::notice, s.str() ) ;
		}
}

//! \todo Coordinates in here are quite probably wrong (good thing
//! nobody relies on them anyway).
void FastaAlnWriter::put_result( const Result& r ) 
{
	if( const Hit* h = hit_to( r ) )
	{
		if( h->has_aln_ref() && h->has_aln_qry() )
		{
			*out_ << '>' << h->sequence() << ' '
				<< h->start_pos()
				<< "-+"[ h->aln_length() > 0 ]
				<< h->start_pos() + abs(h->aln_length()) - 1
				<< '\n' << h->aln_ref() << '\n' 
				<< '>' << r.read().seqid() 
				<< ( r.read().has_trim_right() ? " adapter cut off\n" : "\n" ) 
				<< h->aln_qry() << std::endl ;
		}
	}
}

void FastqWriter::put_result( const Result& rr ) 
{
    const output::Read& r = rr.read() ;
	if( r.has_quality() ) 
	{
		*out_ << ( with_qual_ && r.has_quality() ? '@' : '>' ) << r.seqid() ;
		if( r.has_description() ) *out_ << ' ' << r.description() ;
		for( size_t i = 0 ; i < r.sequence().size() ; i += 50 )
			*out_ << '\n' << r.sequence().substr( i, 50 ) ;
		*out_ << '\n' ;
		if( with_qual_ && r.has_quality() ) 
		{
			*out_ << "+\n" ;
			const std::string& q = r.quality() ;
			for( size_t i = 0 ; i < q.size() ; i += 50 )
			{
				for( size_t j = i ; j != i+50 && j != q.size() ; ++j )
					*out_ << (char)std::min(126, 33 + (uint8_t)q[j]) ;
				*out_ << '\n' ;
			}
		}
	}
}

void TableWriter::put_result( const Result& r )
{
	if( const Hit *h = hit_to( r ) ) {
		int e = r.read().has_trim_right() ? r.read().trim_right() : r.read().sequence().size() ;
		int b = r.read().trim_left() ;
		int diff = h->has_diff_to_next() ? h->diff_to_next() : 9999 ;

		*out_ << e-b << '\t' << r.hit(0).score() << '\t' << diff << '\n' ;
	}
}

void WigCoverageWriter::put_result( const Result& r )
{
	// depth must be known at all
	if( !r.read().depth_size() ) return ;

	// assume one hit, but if there are more, take the best one
	const Hit* h = hit_to( r ) ;

	// must be on forward strand (by construction, in fact)
	if( !h || h->aln_length() <= 0 ) return ;

	*out_ << "fixedStep\tchrom=" << h->sequence() << "\tstart=" << h->start_pos() << "\tstep=1\n" ;

	int alnpos = 0 ;
	int cigar_maj = 0 ;
	size_t cigar_min = 0 ;
	while( alnpos != r.read().depth_size() && cigar_maj != h->cigar_size() ) {
		switch( cigar_op( h->cigar(cigar_maj) ) )
		{
			case Hit::Match:
			case Hit::Mismatch:
			case Hit::Delete:
				*out_ << r.read().depth(alnpos) << '\n' ;
				++alnpos ;
				break ;
			default:
				break ;
		} 
		++cigar_min ;
		while( cigar_maj != h->cigar_size() && cigar_min == cigar_len( h->cigar(cigar_maj) ) ) {
			++cigar_maj ;
			cigar_min = 0 ;
		}
	}
	*out_ << '\n' ;
}

} // namespace

