#ifndef INCLUDED_LOGDOM_H
#define INCLUDED_LOGDOM_H

#include <math.h>

//! \brief representation of floating point number in log domain
//! This class is supposed to behave like ordinary floating point
//! numbers, but the representation is in the log domain, so we don't
//! lose precision due to multiplication of small values.  For ease of
//! conversion, we'll use the phred scale.
class Logdom {
	private:
		double v_ ;

		explicit Logdom( double v ) : v_(v) {}

	public:
		//! \brief constructs a value of one(!)
		//! Why one and not zero?  Because we tend to multiply things,
		//! and one is the neutral element in that case.
		Logdom() : v_(0) {}

		//! \brief converts an ordinary number to log domain
		static Logdom from_float( double v ) { return Logdom( log(v) ) ; }

		//! \brief converts a Phred-quality score to log domain
		//! This will normally result in a very small number; the Phred
		//! score represents the probability of something being wrong.
		static Logdom from_phred( int p ) { return Logdom( -log(10.0)/10.0*p ) ; }

		int to_phred() const { return int( -10.0/log(10.0)*v_ + 0.5 ) ; }
		double to_float() const { return exp( v_ ) ; }

		// multiplication is addition, division is subtraction
		Logdom& operator *= ( Logdom b ) { v_ += b.v_ ; return *this ; }
		Logdom& operator /= ( Logdom b ) { v_ -= b.v_ ; return *this ; }

		Logdom operator * ( Logdom b ) const { return Logdom( v_ + b.v_ ) ; }
		Logdom operator / ( Logdom b ) const { return Logdom( v_ - b.v_ ) ; }

		// addition w/o sacrificing performance is a bit harder
		Logdom operator + ( Logdom b ) const
		{
			return Logdom( v_ >= b.v_
					?   v_ + log1p( exp( b.v_ -   v_ ) ) 
					: b.v_ + log1p( exp(   v_ - b.v_ ) ) ) ;
		}

		Logdom operator - ( Logdom b ) const
		{
			return Logdom( v_ >= b.v_
					?   v_ + log1p( -exp( b.v_ -   v_ ) ) 
					: b.v_ + log1p( -exp(   v_ - b.v_ ) ) ) ; 
		}

		Logdom& operator += ( Logdom b )
		{
			v_ = v_ >= b.v_
				?   v_ + log1p( exp( b.v_ -   v_ ) ) 
				: b.v_ + log1p( exp(   v_ - b.v_ ) ) ;
			return *this ;
		}

		Logdom& operator -= ( Logdom b )
		{
			v_ = v_ >= b.v_
				?   v_ + log1p( -exp( b.v_ -   v_ ) ) 
				: b.v_ + log1p( -exp(   v_ - b.v_ ) ) ; 
			return *this ;
		}

		bool operator >  ( Logdom b ) { return b.v_ <  v_ ; }
		bool operator >= ( Logdom b ) { return b.v_ <= v_ ; }
		bool operator <  ( Logdom b ) { return b.v_ >  v_ ; }
		bool operator <= ( Logdom b ) { return b.v_ >= v_ ; }
} ;

#endif
