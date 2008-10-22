/*
 * FASTA to DNA file converter
 *
 * Reads a FASTA file with one ore more DNA or RNA sequences possibly
 * containing ambiguity codes from stdin, writes it in a compact
 * representation under the filename given as first argument, and writes
 * a human readable index to stdout, which should be stored.
 *
 * We chop the genome into contigs.  Every sequence start also starts a
 * new contig, stretches of more than 2 Ns separate contigs and are not
 * stored.  Ambiguity codes are recognized.
 *
 * A "compact genome" is a file encoding DNA as 4 bits per nucleotide.
 * The bits encode for A,C,T,G from LSB to MSB.  Combinations of bits
 * form ambiguity codes (0 is a gap, used as terminator, 15 is an N),
 * two left rotations reverse-complement a single nucleotide or a pair
 * of nucleotides.  In a byte, the first nucleotide occupies the 4 LSBs.
 *
 * The "compact genome" starts with the ACSII-signature "DNA0", followed
 * by the raw DNA.  Each contig is prepended with a single gap, then the
 * contigs are concatenated, then a final single gap is appended.  The
 * DNA is written out with two nucleotides per byte.  The final byte is
 * padded with a gap if necessary.  Note that every contig is always
 * preceded and succeeded by a gap.
 *
 * An index of the contigs is written to stdout as follows:
 *
 * index := "genome" name "{" scaffold* "}" ";"
 * scaffold := "sequence" name "{" contig* "}" ";"
 * contig := "contig" "{" offset range "}" ";"
 * offset := "offset" uint32 ";"
 * range := "range" uint32 "-" uint32 ";"
 *
 * Where offset is the number of half-bytes to seek into the file to get
 * to the start of the contig, and the range is the place of the contig
 * in the original sequence.
 */

#include <iostream>
#include <fstream>
#include <cctype>

static const int max_num_n = 2 ;

uint8_t decode_fna( char c ) {
	switch( c ) {
		case 'a': case 'A': return 1 ;
		case 'b': case 'B': return 1 ^ 0xf ;
		case 'c': case 'C': return 2 ;
		case 'd': case 'D': return 2 ^ 0xf ;
		case 'g': case 'G': return 8 ;
		case 'h': case 'H': return 8 ^ 0xf ;
		case 't': case 'T': return 4 ;
		case 'u': case 'U': return 4 ;
		case 'v': case 'V': return 4 ^ 0xf ;

		case 'n': case 'N': return 15 ;
		case 'm': case 'M': return  3 ;
		case 'k': case 'K': return 12 ;
		case 'w': case 'W': return  5 ;
		case 's': case 'S': return 10 ;
		case 'r': case 'R': return  9 ;
		case 'y': case 'Y': return  6 ;
		default:            return  0 ;
	}
}

inline bool encodes_nuc( char c ) { return decode_fna(c) != 0 ; }
inline uint8_t to_mask( char c ) { return decode_fna(c) & 0xf ; }

class FastaDecoder
{
	private:
		void (FastaDecoder::*state)( char ) ;

		unsigned position ;
		unsigned source_position ;
		unsigned num_n  ;
		uint8_t one_nt ;
		std::ostream& out ;
		std::ostream& dna ;
		
		
	public:
		FastaDecoder( unsigned p, std::ostream& dna_, std::ostream& out_ ) 
			: state( &FastaDecoder::s_idle ), position( p ), num_n( 0 ), one_nt( 0 )
			, dna( dna_ ), out( out_ ) {}

		void step( char c ) { (this->*state)( c ) ; }

		void begin_new_sequence()
		{
			source_position = 0 ;
			out << "\tsequence \"" ; 
		}
		void finish_sequence() { out << "\t} ;\n" ; }
		void begin_new_contig() 
		{ 
			put_nt( '-' ) ;
			num_n = 0 ;
			out << "\t\tcontig {\n\t\t\toffset "
				<< position << " ;\n\t\t\trange " 
				<< source_position << " - " ; 
		}
		void finish_contig()
		{ 
			out << source_position << " ;\n\t\t} ;\n" ;
		}
		void finish_genome()
		{
			put_nt( '-' ) ;
			if( position & 1 ) put_nt( '-' ) ;
			out << "} ;\n" ;
		}

		void put_char_lit( char c ) { out << c ; /* XXX */ }

		void put_nt( char c )
		{
			if( c == 'n' || c == 'N' ) 
			{
				++num_n ;
			}
			else
			{
				for( ; num_n ; --num_n ) put_nt_( '-' ) ;
				put_nt_( c ) ;
			}
		}

		void put_nt_( char c )
		{
			if( position & 1 )
			{
				dna.put( (to_mask(c) << 4) | one_nt ) ;
			}
			else
			{
				one_nt = to_mask(c) ;
			}
			++position ;
		}

		void s_idle( char c )
		{
			if( c == '>' ) {
				begin_new_sequence() ;
				state = &FastaDecoder::s_name ;
			}
		}

		void s_name( char c )
		{
			if( !c ) finish_genome() ;
			else if( std::isspace(c) || c == '\n' )
			{
				out << "\" {\n" ;
				if( c == '\n' ) state = &FastaDecoder::s_seq_begin ;
				else state = &FastaDecoder::s_info ;
			}
			else put_char_lit( c ) ;
		}

		void s_info( char c )
		{
			if( !c ) finish_genome() ; 
			else if( c == '\n' ) state = &FastaDecoder::s_seq_begin ;
		}

		void s_seq_begin( char c )
		{
			if( !c ) finish_genome() ;
			else if( c == ';' ) state = &FastaDecoder::s_info ;
			else if( c == '>' ) {
				begin_new_sequence() ;
				state = &FastaDecoder::s_name ;
			}
			else if( c == 'n' || c == 'N' ) {
				if( num_n == max_num_n ) {
					source_position += num_n+1 ;
					num_n = 0 ;
					state = &FastaDecoder::s_contig_gap ;
				}
				else
				{
					++source_position ;
					put_nt( c ) ;
				}
			}
			else if( encodes_nuc( c ) ) {
				begin_new_contig() ;
				++source_position ;
				put_nt( c ) ;
				state = &FastaDecoder::s_sequence ;
			}
		}

		void s_sequence( char c )
		{
			if( !c ) {
				finish_contig() ;
				finish_sequence() ;
				finish_genome() ;
			}
			else if( c == 'N' || c == 'n' && num_n == max_num_n ) {
				finish_contig() ;
				source_position += num_n+1 ;
				num_n = 0 ;
				state = &FastaDecoder::s_contig_gap ;
			}
			else if( encodes_nuc( c ) ) {
				++source_position ;
				put_nt( c ) ;
			}
			else if( c == '>' ) {
				finish_contig() ;
				finish_sequence() ;
				begin_new_sequence() ;
				state = &FastaDecoder::s_name ;
			}
		}

		void s_contig_gap( char c ) 
		{
			if( !c ) {
				finish_sequence() ;
				finish_genome() ;
			}
			else if( c == '>' ) {
				finish_sequence() ;
				begin_new_sequence() ;
				state = &FastaDecoder::s_name ;
			}
			else if( encodes_nuc( c ) && c != 'n' && c != 'N' ) {
				begin_new_contig() ;
				++source_position ;
				put_nt( c ) ;
				state = &FastaDecoder::s_sequence ;
			}
			else if( c == 'n' || c == 'N' ) {
				++source_position ;
			}
		}
} ;




int main( int argc, const char *argv[] )
{
	if( argc != 2 ) return 1 ;
	std::ofstream dna_file( argv[1] ) ;
	dna_file.write( "DNA0", 4 ) ;
	std::cout << "genome \"" << argv[1] << "\" {\n" ;
	FastaDecoder fd( 8, dna_file, std::cout ) ;
	for( int c ; (c = std::cin.get()) != std::istream::traits_type::eof() ; ) fd.step( c ) ;
	fd.step( 0 ) ;
}


