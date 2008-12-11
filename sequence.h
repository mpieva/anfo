#ifndef INCLUDED_SEQUENCE_H
#define INCLUDED_SEQUENCE_H

#include <climits>
#include <iostream>
#include <vector>

#include <stdint.h>

//! \brief Test macro whether this is a "small" system
// A system is considered "small" for our purposes, if long has no more
// than 32 bits.  This affects mainly the organization of Judy arrays.
#define SMALL_SYS (ULONG_MAX < 0x100000000)

/*! \defgroup typedefs Useful typedefs
 * These typedefs are used mostly for their documentation value;  C++
 * unfortunately won't be able to check their differences
 * most of the time.
 *
 * @{ */

//! \brief Type used to store short DNA sequences.
// Oligos are stored with two bits per bases where [0,1,2,3] mean
// [A,C,T,G].  The two LSBs contain the first base, up to 32 bases can
// be stored.
typedef uint64_t Oligo ;

//! \brief Type used to store a single nucleotide base.
// We encode [A,C,T,G] as [0,1,2,3], same as in \c Oligo.
typedef uint8_t  Nucleotide ;	// 0,1,2,3 == A,C,T,G

//! \brief Type used to store ambiguity codes.
// We encode [A,C,T,G] as [1,2,4,8].  A one shifted by the \c Nucleotide
// code gives the approriate ambiguity code, other codes can be created
// by combining bases with a logical OR.  Zero encodes a gap, 15 is an
// N.
typedef uint8_t  Ambicode ;		// 1,2,4,8 == A,C,T,G

//! @}

//! \brief Decodes a character to a nucleotide.
// Takes an arbitrary character and determines the IUPAC ambiguity code
// it stands for.  Small and capital letters are understood, everything
// unrecognized becomes a gap.

inline Ambicode to_ambicode( char c ) {
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

inline char from_ambicode( Ambicode a ) { return "-ACMTWYHGRSVKDBN"[a] ; }

//! \brief Complement an ambiguity code.
inline Ambicode complement( Ambicode w )
{ return 0xf & (w << 2 | w >> 2) ; }

//! \brief Reverse-complements a pair of ambiguity codes.
// \internal
inline uint8_t reverse_complement( uint8_t xy )
{
	return (xy & 0x03) << 6 |
		   (xy & 0x0c) << 2 |
		   (xy & 0x30) >> 2 |
		   (xy & 0xc0) >> 6 ;
}

//! \brief Checks whether a character codes for a nucleotide.
// This is euivalent to decoding the character and checking that it
// doesn't encode a gap.
inline bool encodes_nuc( char c ) { return to_ambicode(c) != 0 ; }

//! Dumb pointer to DNA
// A pointer with sub-byte precision is needed for our mmaped genomes.
// The assumption here is that this fits into 64 bits, which is true on
// any ix86_64 bit that doesn't implement full 64bit addresses, and that
// will be all of them for quite some time to come.  It's also true on a
// 32 bit system, obviously.  We also steal another bit to encode the
// strand we're pointing to.
//
// To make arithmetic easier, the encoding is as follows:  The main
// pointer is converted to an int64_t, it is then shifted left and the
// sub-byte index (just one bit) is shifted in.  Then the sign in
// inverted if and only if this is a pointer to the RC strand.  This way,
// arithmetic works on both strands the same way.  To dereference a
// pointer, the absolute value of the number is taken, interpreted as
// signed 63 bit value, sign extended and reinterpreted as pointer.
class DnaP
{
	private:
		int64_t p_ ;

	public:
		//! constructs a DNA pointer
		// The pointer is initialized from a pointer to bytes,
		// understood to point to the ambiguity code in the 4 LSBs.  It
		// is then converted to a reverse pointer and finally the offset
		// is added, taking directionality into account.
		// \param p pointer to ambiguity-encoded DNA
		// \param complement set to true to make a pointer to the
		//                   reverse-complemented strand
		// \param off offset to add to the final pointer
		explicit DnaP( uint8_t *p = 0, bool complement = false, int off = 0 ) 
		{
			p_ = (reinterpret_cast<int64_t>(p) << 1) & std::numeric_limits<int64_t>::max() ;
			if( complement ) p_ = -p_ ;
			p_ += off ;
		}

		void assign( uint8_t *p = 0, bool complement = false, int off = 0 )
		{
			p_ = (reinterpret_cast<int64_t>(p) << 1) & std::numeric_limits<int64_t>::max() ;
			if( complement ) p_ = -p_ ;
			p_ += off ;
		}

		Ambicode operator * () const { return (*this)[0] ; }
		Ambicode operator [] ( int64_t ix ) const {
			int64_t p = p_ + ix ;
			uint8_t w = *( reinterpret_cast<uint8_t*>( (std::abs(p) << 1) >> 2 ) ) ;
			w = 0xf & ( p & 1 ? w >> 4 : w ) ;
			return p < 0 ? complement(w) : w ;  
		}

		operator bool() const { return p_ ; }

		uint8_t *unsafe_ptr() const
		{ return reinterpret_cast<uint8_t*>( (std::abs(p_) << 1) >> 2 ) ; }

		DnaP &operator ++ () { p_++ ; return *this ; }
		DnaP &operator -- () { p_-- ; return *this ; }

		DnaP &operator += ( int64_t  o ) { p_ += o ; return *this ; }
		DnaP &operator += ( uint32_t o ) { p_ += o ; return *this ; }
		DnaP &operator += ( int32_t  o ) { p_ += o ; return *this ; }
		DnaP &operator -= ( int64_t  o ) { p_ -= o ; return *this ; }
		DnaP &operator -= ( uint32_t o ) { p_ -= o ; return *this ; }
		DnaP &operator -= ( int32_t  o ) { p_ -= o ; return *this ; }

		// hack to make this compatible with Judy arrays on both 32 and
		// 64 bit machines
#if SMALL_SYS
		unsigned long get() const { return p_ & ULONG_MAX ; }
		unsigned long high() const { return p_ >> sizeof( unsigned long ) * CHAR_BIT ; }
#else
		unsigned long get() const { return p_ ; }
#endif

		friend inline int64_t operator -  ( const DnaP &a, const DnaP &b ) { return a.p_  - b.p_ ; }
		friend inline bool    operator == ( const DnaP &a, const DnaP &b ) { return a.p_ == b.p_ ; }
		friend inline bool    operator != ( const DnaP &a, const DnaP &b ) { return a.p_ != b.p_ ; }
		
		friend inline std::ostream& operator << ( std::ostream& s, const DnaP &p )
		{ return s << std::hex << p.p_ ; }
} ;

inline DnaP operator + ( const DnaP& a, int64_t  o ) { DnaP b = a ; return b += o ; }
inline DnaP operator + ( const DnaP& a, int32_t  o ) { DnaP b = a ; return b += o ; }
inline DnaP operator + ( const DnaP& a, uint32_t o ) { DnaP b = a ; return b += o ; }
inline DnaP operator - ( const DnaP& a, int64_t  o ) { DnaP b = a ; return b -= o ; }
inline DnaP operator - ( const DnaP& a, int32_t  o ) { DnaP b = a ; return b -= o ; }
inline DnaP operator - ( const DnaP& a, uint32_t o ) { DnaP b = a ; return b -= o ; }

//! \brief Wrapper for easier output of gap-terminated sequences.
// \internal
struct Sequ
{
	Sequ( const DnaP p, int length = std::numeric_limits<int>::max() )
		: p_(p), l_(length) {}

	DnaP p_ ;
	int l_ ;
} ;

inline std::ostream& operator << ( std::ostream& s, const Sequ& d )
{
	for( DnaP p = d.p_, q = d.p_ + d.l_ ; *p && p != q ; ++p )
		s << from_ambicode( *p ) ;
	return s ;
}

//! \brief Sequence transformed into same representation used for genomes.
// We store the sequence in both forward and reverse direction so we can
// eaily refer to either strand by a DnaP.  This doesn't save space, but
// it does make the representation uniform.  Again, both sequences are
// gap-terminated on either end.
class PreparedSequence
{
	private:
		std::vector< uint8_t > seq ;

	public:
		PreparedSequence( const char* p ) 
		{
			// start out with a gap symbol
			int n = 0 ;
			for( Ambicode last = 0 ;; )
			{
				// this is intentional: we actually encode the final NUL
				while( *p && !encodes_nuc( *p ) ) ++p ;
				seq.push_back( last | to_ambicode( *p ) << 4 ) ;
				if( !*p ) break ;
				++n ;
				++p ;

				// also intentional: if we hit the final NUL here, we
				// encode it twice to achieve padding
				while( *p && !encodes_nuc( *p ) ) ++p ;
				last = to_ambicode( *p ) ;
				if( *p ) ++p, ++n ;
			}
		}

		DnaP forward() { return DnaP( &seq[0], false, 1 ) ; }
		DnaP reverse() { return DnaP( &seq[0], true, -1 ) ; }
} ;

#endif

