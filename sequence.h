#ifndef INCLUDED_SEQUENCE_H
#define INCLUDED_SEQUENCE_H

#include <climits>
#include <iostream>
#include <vector>

#include <stdint.h>

#if ULONG_MAX < 0x100000000
#define SMALL_SYS 1
#else
#undef SMALL_SYS 
#endif

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

inline uint8_t reverse_complement( uint8_t xy )
{
	return (xy & 0x03) << 6 |
		   (xy & 0x0c) << 2 |
		   (xy & 0x30) >> 2 |
		   (xy & 0xc0) >> 6 ;
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
		explicit DnaP( uint8_t *p = 0, int off = 0 ) 
			: p_( (reinterpret_cast<int64_t>(p) << 1) + off ) {}

		Ambicode operator[]( int64_t ix ) const {
			int64_t p = p_ + ix ;
			uint8_t w = *( reinterpret_cast<uint8_t*>(p >> 1) ) ;
			if( p & 1 ) return (w >> 4) & 0xf ; else return w & 0xf ;
		}

		operator bool() const { return p_ ; }

		void assign( uint8_t *p = 0, int off = 0 )
			{ p_ = (reinterpret_cast<int64_t>(p) << 1) + off ; }

		uint8_t *unsafe_ptr() const { return reinterpret_cast<uint8_t*>( p_ >> 1 ) ; }

		DnaP &operator += ( int64_t o ) { p_ += o ; return *this ; }
		DnaP &operator -= ( int64_t o ) { p_ += o ; return *this ; }

		// hack to make this compatible with Judy arrays on both 32 and
		// 64 bit machines
#if SMALL_SYS
		unsigned long get() const { return p_ & ULONG_MAX ; }
		unsigned long high() const { return p_ >> sizeof( unsigned long ) * CHAR_BIT ; }
#else
		unsigned long get() const { return p_ ; }
#endif

		friend inline int64_t operator - ( DnaP a, DnaP b ) { return a.p_ - b.p_ ; }
} ;

inline DnaP operator + ( DnaP a, int64_t o ) { DnaP b = a ; return b += o ; }
inline DnaP operator - ( DnaP a, int64_t o ) { DnaP b = a ; return b -= o ; }

class PreparedSequence
{
	private:
		std::vector< uint8_t > forward_seq ;
		std::vector< uint8_t > reverse_seq ;

		DnaP forward_ ;
		DnaP reverse_ ;
	public:
		PreparedSequence( const char* p ) 
		{
			// start out with a gap symbol
			int n = 0 ;
			for( Ambicode last = 0 ;; )
			{
				// this is intentional: we actually encode the final NUL
				while( *p && !encodes_nuc( *p ) ) ++p ;
				forward_seq.push_back( last | to_ambicode( *p ) << 4 ) ;
				if( !*p ) break ;
				++n ;
				++p ;

				// also intentional: if we hit the final NUL here, we
				// encode it twice to achieve padding
				while( *p && !encodes_nuc( *p ) ) ++p ;
				last = to_ambicode( *p ) ;
				if( *p ) ++p, ++n ;
			}
			std::clog << "recognized " << n << " bases" << std::endl ;

			std::transform( forward_seq.rbegin(), forward_seq.rend(),
					std::back_inserter( reverse_seq ), reverse_complement ) ;

			forward_.assign( &forward_seq[0], 1 ) ; 
			reverse_.assign( &reverse_seq[0], (n % 2) + 1 ) ;
		}

		DnaP forward() { return forward_ ; }
		DnaP reverse() { return reverse_ ; }
} ;

#endif

