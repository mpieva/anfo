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

#ifndef INCLUDED_LOGDOM_H
#define INCLUDED_LOGDOM_H

#include <math.h>
#include <stdint.h>

//! \brief representation of floating point number in log domain
//! This class is supposed to behave like ordinary floating point
//! numbers, but the representation is in the log domain, so we don't
//! lose precision due to multiplication of many small values.  For ease
//! of conversion, we'll use the phred scale.  That means more work in
//! additions, but less on conversion.  Assuming \f$ q_2 \geq q_1 \f$
//! the formula for additions is:
//!
//! \f[ q_s = q_1 - \frac{10}{\ln 10} \ln \left( 1 + \exp \left( \frac{\ln 10}{10} (q_1-q_2) \right) \right) \f]
//!
//! It turns out that IEEE floating point handles corner cases well:
//! from_float(0) gives an internal representation of \f$ -\inf \f$,
//! which works fine in calculations, except that adding two such
//! numbers gives NaN (so don't do that).

class Logdom {
	private:
		double v_ ;

		explicit Logdom( double v ) : v_(v) {}

	public:
		//! \brief constructs a value of one(!)
		//! Why one and not zero?  Because we tend to multiply things,
		//! and one is the neutral element in that case.
		Logdom() : v_(0) {}

		bool is_finite() const { return isfinite( v_ ) ; }

		//! \brief calculates log( 1 + exp x )
		//! Calculation is done in a way that preserves accuracy, the
		//! base of the logarithm is chosen to coincide with the Phred
		//! scale.
		//! \todo This function need to be accurate only for x < 0; it
		//!       could probably be approximated using a Taylor series.
		static double ld1pexp10( double x ) { return -10.0 / log(10.0) * log1p(  exp( -log(10.0)/10.0*x ) ) ; }

		//! \brief calculates log( 1 - exp x )
		//! \see ld1pexp10()
		static double ld1mexp10( double x ) { return -10.0 / log(10.0) * log1p( -exp( -log(10.0)/10.0*x ) ) ; }

		static Logdom null() { return from_float( 0 ) ; }
		static Logdom one() { return from_float( 1 ) ; }

		//! \brief converts an ordinary number to log domain
		static Logdom from_float( double v ) { return Logdom( -10 * log(v) / log(10.0) ) ; }

		//! \brief converts a Phred-quality score to log domain
		//! This will normally result in a very small number; the Phred
		//! score represents the probability of something being wrong.
		static Logdom from_phred( int p ) { return Logdom( p ) ; }

		int to_phred() const { return int( v_ + 0.5 ) ; }
		uint8_t to_phred_byte() const { int p = to_phred() ; return p > 255 ? 255 : p ; }
		double to_float() const { return exp( -log(10.0) * v_ / 10.0 ) ; }

		// multiplication is addition, division is subtraction
		Logdom& operator *= ( Logdom b ) { v_ += b.v_ ; return *this ; }
		Logdom& operator /= ( Logdom b ) { v_ -= b.v_ ; return *this ; }

		Logdom operator * ( Logdom b ) const { return Logdom( v_ + b.v_ ) ; }
		Logdom operator / ( Logdom b ) const { return Logdom( v_ - b.v_ ) ; }

		// addition w/o sacrificing precision is a bit harder
		Logdom operator + ( Logdom b ) const
		{
			return Logdom( v_ <= b.v_
					?   v_ + ld1pexp10( b.v_ -   v_ ) 
					: b.v_ + ld1pexp10(   v_ - b.v_ ) ) ;
		}

		Logdom operator - ( Logdom b ) const
		{
			return Logdom( v_ <= b.v_
					?   v_ + ld1mexp10( b.v_ -   v_ ) 
					: b.v_ + ld1mexp10(   v_ - b.v_ ) ) ; 
		}

		Logdom& operator += ( Logdom b )
		{
			v_ = v_ <= b.v_
				?   v_ + ld1pexp10( b.v_ -   v_ ) 
				: b.v_ + ld1pexp10(   v_ - b.v_ ) ;
			return *this ;
		}

		Logdom& operator -= ( Logdom b )
		{
			v_ = v_ <= b.v_
				?   v_ + ld1mexp10( b.v_ -   v_ ) 
				: b.v_ + ld1mexp10(   v_ - b.v_ ) ; 
			return *this ;
		}

		bool operator >  ( Logdom b ) { return b.v_ >  v_ ; }
		bool operator >= ( Logdom b ) { return b.v_ >= v_ ; }
		bool operator <  ( Logdom b ) { return b.v_ <  v_ ; }
		bool operator <= ( Logdom b ) { return b.v_ <= v_ ; }
} ;

inline Logdom operator + ( double a, Logdom b ) { return Logdom::from_float( a ) + b ; }
inline Logdom operator + ( Logdom a, double b ) { return a + Logdom::from_float( b ) ; }
inline Logdom operator - ( double a, Logdom b ) { return Logdom::from_float( a ) - b ; }
inline Logdom operator - ( Logdom a, double b ) { return a - Logdom::from_float( b ) ; }
inline Logdom operator * ( double a, Logdom b ) { return Logdom::from_float( a ) * b ; }
inline Logdom operator * ( Logdom a, double b ) { return a * Logdom::from_float( b ) ; }
inline Logdom operator / ( double a, Logdom b ) { return Logdom::from_float( a ) / b ; }
inline Logdom operator / ( Logdom a, double b ) { return a / Logdom::from_float( b ) ; }

template <typename A, typename B> 
A lerp( B p, A a, A b )
{ return (1.0-p) * a + p * b ; }

#endif
