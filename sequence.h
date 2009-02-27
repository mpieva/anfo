#ifndef INCLUDED_SEQUENCE_H
#define INCLUDED_SEQUENCE_H

#include <cassert>
#include <climits>
#include <istream>
#include <vector>

#include <stdint.h>

//! \brief Test macro whether this is a "small" system
//!
//! A system is considered "small" for our purposes, if long has no more
//! than 32 bits.  This affects mainly the organization of Judy arrays.
#define SMALL_SYS (ULONG_MAX < 0x100000000)

/*! \defgroup typedefs Useful typedefs
 * These typedefs are used mostly for their documentation value;  C++
 * unfortunately won't be able to check their differences
 * most of the time.
 *
 * @{ */

//! \brief Type used to store short DNA sequences.
//!
//! Oligos are stored with two bits per bases where [0,1,2,3] mean
//! [A,C,T,G].  The two LSBs contain the first base, up to 32 bases can
//! be stored.
typedef uint64_t Oligo ;

//! \brief Type used to store a single nucleotide base.
//!
//! We encode [A,C,T,G] as [0,1,2,3], same as in \c Oligo.
typedef uint8_t  Nucleotide ;	// 0,1,2,3 == A,C,T,G

//! \brief Type used to store ambiguity codes.
//!
//! We encode [A,C,T,G] as [1,2,4,8].  A one shifted by the \c Nucleotide
//! code gives the approriate ambiguity code, other codes can be created
//! by combining bases with a logical OR.  Zero encodes a gap, 15 is an
//! N.
typedef uint8_t  Ambicode ;		// 1,2,4,8 == A,C,T,G

//! @}

//! \brief Decodes a character to a nucleotide.
//!
//! Takes an arbitrary character and determines the IUPAC ambiguity code
//! it stands for.  Small and capital letters are understood, everything
//! unrecognized becomes a gap.

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
//! \internal
inline uint8_t reverse_complement( uint8_t xy )
{
	return (xy & 0x03) << 6 |
		   (xy & 0x0c) << 2 |
		   (xy & 0x30) >> 2 |
		   (xy & 0xc0) >> 6 ;
}

//! \brief Checks whether a character codes for a nucleotide.
//! This is euivalent to decoding the character and checking that it
//! doesn't encode a gap.
inline bool encodes_nuc( char c ) { return to_ambicode(c) != 0 ; }

//! Dumb pointer to DNA
//! A pointer with sub-byte precision is needed for our mmaped genomes.
//! The assumption here is that this fits into 64 bits, which is true on
//! any ix86_64 bit that doesn't implement full 64bit addresses, and that
//! will be all of them for quite some time to come.  It's also true on a
//! 32 bit system, obviously.  We also steal another bit to encode the
//! strand we're pointing to.
//!
//! To make arithmetic easier, the encoding is as follows:  The main
//! pointer is converted to an int64_t, it is then shifted left and the
//! sub-byte index (just one bit) is shifted in.  Then the sign in
//! inverted if and only if this is a pointer to the RC strand.  This way,
//! arithmetic works on both strands the same way.  To dereference a
//! pointer, the absolute value of the number is taken, interpreted as
//! signed 63 bit value, sign extended and reinterpreted as pointer.
class DnaP
{
	private:
		int64_t p_ ;

	public:
		//! constructs a DNA pointer
		//! The pointer is initialized from a pointer to bytes,
		//! understood to point to the ambiguity code in the 4 LSBs.  It
		//! is then converted to a reverse pointer and finally the offset
		//! is added, taking directionality into account.
		//! \param p pointer to ambiguity-encoded DNA
		//! \param complement set to true to make a pointer to the
		//!                   reverse-complemented strand
		//! \param off offset to add to the final pointer
		explicit DnaP( const uint8_t *p = 0, bool complement = false, int off = 0 ) { assign( p, complement, off ) ; }

		void assign( const uint8_t *p = 0, bool complement = false, int off = 0 )
		{
			p_ = (reinterpret_cast<int64_t>(p) << 1) & std::numeric_limits<int64_t>::max() ;
			if( complement ) p_ = -p_ ;
			p_ += off ;
			assert( unsafe_ptr() == p ) ;
		}

		Ambicode operator * () const { return (*this)[0] ; }
		Ambicode operator [] ( int64_t ix ) const {
			int64_t p = p_ + ix ;
			int64_t p2 = std::abs(p) ;
			if( (p2 << 1) < 0 ) p2 |= std::numeric_limits<int64_t>::min() ;
			uint8_t w = *( reinterpret_cast<uint8_t*>( p2 >> 1 ) ) ;
			w = 0xf & ( p & 1 ? w >> 4 : w ) ;
			return p < 0 ? complement(w) : w ;  
		}

		operator const void*() const { return (const void*)p_ ; }

		bool is_reversed() const { return p_ < 0 ; }
		DnaP abs() const { DnaP q ; q.p_ = std::abs(p_) ; return q ; }
		DnaP reverse() const { DnaP q ; q.p_ = -p_ ; return q ; }
		DnaP reverse_if( bool p ) const { return p ? reverse() : *this ; }

		//! \brief returns a pointer to the underlying storage.
		//! \internal
		//! If you call this function for anything, you're on your own.
		const uint8_t *unsafe_ptr() const
		{ 
			int64_t p2 = std::abs(p_) ;
			if( (p2 << 1) < 0 ) p2 |= std::numeric_limits<int64_t>::min() ;
			return reinterpret_cast<uint8_t*>( p2 >> 1 ) ;
		}

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

class QDnaP
{
	private:
		const uint16_t *p_ ;

	public:
		explicit QDnaP( const uint16_t *p = 0 ) : p_(p) {}

		QDnaP &operator = ( const QDnaP& p ) { p_ = p.p_ ; return *this ; }

		Ambicode operator * () const { return *p_ & 0xf ; }
		uint8_t qual( size_t ix = 0 ) const { return p_[ix] >> 8 ; }

		Ambicode operator [] ( size_t ix ) const { return p_[ix] & 0xf ; }

		operator const void*() const { return p_ ; }

		QDnaP &operator ++ () { p_++ ; return *this ; }
		QDnaP &operator -- () { p_-- ; return *this ; }

		QDnaP &operator += ( int64_t  o ) { p_ += o ; return *this ; }
		QDnaP &operator += ( uint32_t o ) { p_ += o ; return *this ; }
		QDnaP &operator += ( int32_t  o ) { p_ += o ; return *this ; }
		QDnaP &operator -= ( int64_t  o ) { p_ -= o ; return *this ; }
		QDnaP &operator -= ( uint32_t o ) { p_ -= o ; return *this ; }
		QDnaP &operator -= ( int32_t  o ) { p_ -= o ; return *this ; }

#if SMALL_SYS
		unsigned long high() const { return 0 ; }
#endif
		unsigned long get() const { return (unsigned long)p_ ; }

		friend inline size_t operator -  ( const QDnaP &a, const QDnaP &b ) { return a.p_  - b.p_ ; }
		friend inline bool   operator == ( const QDnaP &a, const QDnaP &b ) { return a.p_ == b.p_ ; }
		friend inline bool   operator != ( const QDnaP &a, const QDnaP &b ) { return a.p_ != b.p_ ; }
		
		friend inline std::ostream& operator << ( std::ostream& s, const QDnaP &p )
		{ return s << std::hex << p.p_ ; }
} ;

inline QDnaP operator + ( const QDnaP& a, int64_t  o ) { QDnaP b = a ; return b += o ; }
inline QDnaP operator + ( const QDnaP& a, int32_t  o ) { QDnaP b = a ; return b += o ; }
inline QDnaP operator + ( const QDnaP& a, uint32_t o ) { QDnaP b = a ; return b += o ; }
inline QDnaP operator - ( const QDnaP& a, int64_t  o ) { QDnaP b = a ; return b -= o ; }
inline QDnaP operator - ( const QDnaP& a, int32_t  o ) { QDnaP b = a ; return b -= o ; }
inline QDnaP operator - ( const QDnaP& a, uint32_t o ) { QDnaP b = a ; return b -= o ; }

//! \brief Wrapper for easier output of gap-terminated sequences.
//! \internal
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

//! \brief sequence with quality scores
//! We store four qualities for the four possible bases.  If the
//! sequence was read from a 4Q file, that's okay.  If it came from
//! FASTQ or even FASTA, they are just an extrapolation.  That's why we
//! keep a flag saying which data set is actually valid.  Qualities are
//! actual probabilities (they should sum to one), because calculation
//! would involve logarithms anyway and we needn't fear loss of
//! precision here.
//! 
//! In theory, we don't need to store the bases themselves.  However, to
//! support weird base callers and to keep the legacy indexing code
//! working, we still store them.
class QSequence
{
	public:
		//! \brief An ambiguity code and four qualities together.
		//! The qualities are invalid if the ambiguity code encodes a
		//! gap.
		struct Base {
			uint8_t ambicode ;
			float   qualities[4] ;

			Base() : ambicode(0) { qualities[0] = qualities[1] = qualities[2] = qualities[3] = 0 ; }
			Base( uint8_t a, int q_score ) ;
		} ;

		enum Validity { bases_only, bases_with_quality, qualities_only, bases_with_qualities } ;
	private:
		std::vector< Base > seq_ ;
		std::string name_ ;
		std::string description_ ;
		Validity validity_ ;

	public:
		typedef std::vector< Base >::const_iterator const_iterator ;
		typedef std::vector< Base >::iterator iterator ;
		typedef Base const *const_pointer ;
		typedef Base const &const_reference ;
		typedef Base pointer ;
		typedef Base &reference ;
		typedef Base value_type ;

		//! \brief creates an empty sequence
		//! The sequence is terminated correctly by two gaps and has no
		//! bases.  It is basically considered invalid and used as a
		//! placeholder where necessary.
		QSequence() : seq_(), name_(), description_(), validity_( bases_only )
		{
			seq_.push_back( Base() ) ; 
			seq_.push_back( Base() ) ; 
		}

		//! \brief creates a sequence from ACSII encoded bases
		//! All qualities are set to a reasonable constant, which is
		//! equivalent to a Q score of 30 by default.  Useful mostly for
		//! debugging.
		QSequence( const char* p, int q_score = 30 ) ;

		const_pointer start() const { return &seq_[1] ; }
		unsigned length() const { return seq_.size() - 2 ; }

		const std::string &get_name() const { return name_ ; }
		const std::string &get_descr() const { return description_ ; } 

		Validity get_validity() const { return validity_ ; }

		std::string as_string() const {
			const_iterator a = seq_.begin(), b = seq_.end() ;
			std::string r ;
			for( ++a, --b ; a != b ; ++a ) r.push_back( from_ambicode( a->ambicode ) ) ;
			return r ;
		}

		reference operator [] ( size_t ix ) { return seq_[1+ix] ; }
		const_reference operator [] ( size_t ix ) const { return seq_[1+ix] ; }

		// Ambicode operator [] ( size_t ix ) const { return seq_[1+ix]->ambicode ; }

		//! \todo How do you calculate the quality value of an ambiguity
		//! code given the qualities of the individual bases?  Just
		//! adding them doesn't seem right somehow.
		
		// uint8_t qual( size_t ix ) const { return seq_[1+ix] >> 8 ; }
		// void qual( size_t ix, uint8_t q ) { seq_[1+ix] = seq_[1+ix] & 0xff | (uint16_t)q << 8 ; }

		friend std::istream& read_fastq( std::istream& s, QSequence& qs, bool solexa_scores ) ;
} ;

//! \brief reads a sequence from a FASTA, FASTQ or 4Q file
//! To unify all formats (well-defined or not...) this function accepts
//! both '@' and '>' as header delimiter.  Any lines starting ';' and
//! following the header are treated as extended description.  Quality
//! scores are optional and introduced by '+'; they can be either ASCII
//! (U Rockefeller idiocy), or 33-based raw Phred scale values (Sanger
//! standard), or 64-based Solexa scale values (Solexa idiocy).  The tag
//! "SQ " (4Q definition) is ignored in a sequence line.  Quality scores
//! can also be introduced by '*', in which case the quality line tag is
//! read and the remainder is parsed as quality scores in 4Q format.  If
//! no bases are found, but four Q scores are available, bases are
//! "called" internally.
//!
//!	\todo Add sensible error checking (e.g. consistent line lengths, no
//!	      additional junk, no wrong characters).
//!	\todo Going from log domain to probabilities and back to log domain
//!	      may cause loss of precision.  Maybe unaltered Q scores should
//!	      be stored, too.
//!
//! \param s stream to read from
//! \param qs sequence-with-quality to read into
//! \param solexa_scores 
//!     If set, quality scores are converted from Solexa to Phred
//!     conventions and the origin for raw scores is 64.  Else the
//!     origin for raw scores is 33 and no conversion is done.
//! \return the stream \c s
std::istream& read_fastq( std::istream& s, QSequence& qs, bool solexa_scores = false ) ;

#endif

