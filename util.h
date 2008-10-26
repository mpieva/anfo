#ifndef INCLUDED_UTIL_H
#define INCLUDED_UTIL_H

#include <cerrno>
#include <string>

template< typename T >
void throw_errno_if_minus1( T x, const char* a, const char* b = 0 )
{ throw_errno_if_eq( x, (T)(-1), a, b ) ; }

template< typename T >
void throw_errno_if_null( T x, const char* a, const char* b = 0 )
{ throw_errno_if_eq( x, (T)0, a, b ) ; }

template< typename T >
void throw_errno_if_eq( T x, T y, const char* a, const char* b = 0 )
{
	if( x == y )
	{
		std::string msg ;
#ifdef _GNU_SOURCE
		msg.append( program_invocation_short_name ).append( ": " ) ;
#endif
		msg.append( strerror( errno ) ).append( " while " ).append( a ) ;
		if( b ) msg.append( " " ).append( b ) ;
		throw msg ;
	}
}


int main( int argc, const char * const argv[] ) ;

#endif
