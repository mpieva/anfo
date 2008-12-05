#ifndef INCLUDED_SEQUENCE_H
#define INCLUDED_SEQUENCE_H

#include <vector>

#include <stdint.h>

// Useful typedefs.  These are used mostly for their documentation
// value, C++ unfortunately won't be able to check their differences
// most of the time.

typedef uint64_t Oligo ;
typedef uint8_t  Nucleotide ;	// 0,1,2,3 == A,C,T,G
typedef uint8_t  Ambicode ;		// 1,2,4,8 == A,C,T,G

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

inline bool encodes_nuc( char c ) { return to_ambicode(c) != 0 ; }

// dumb pointer to DNA
// A pointer with sub-byte precision is needed.  The assumption is that
// this fits into 64 bits, which is true on any ix86_64 bit that doesn't
// implement full 64bit addresses, and that will be all of them for
// quite some time to come.

class DnaP
{
	private:
		int64_t p_ ;

	public:
		explicit DnaP( void const *p = 0, int off = 0 ) 
			: p_( (reinterpret_cast<int64_t>(p) << 1) + off ) {}

		Ambicode operator[]( uint32_t ix ) const {
			int64_t p = p_ + ix ;
			uint8_t w = *( reinterpret_cast<uint8_t*>(p >> 1) ) ;
			if( p & 1 ) return (w >> 4) & 0xf ; else return w & 0xf ;
		}

		operator void const * () const { return reinterpret_cast<void*>( p_ ) ; }

		void assign( void const *p = 0, int off = 0 )
			{ p_ = (reinterpret_cast<int64_t>(p) << 1) + off ; }

		void *unsafe_ptr() const { return reinterpret_cast<void*>( p_ >> 1 ) ; }
} ;

class PreparedSequence
{
	private:
		std::vector< uint8_t > forward ;
		std::vector< uint8_t > reverse ;

	public:
		PreparedSequence( const char* ) ;
} ;

#endif

