#ifndef INCLUDED_UTIL_H
#define INCLUDED_UTIL_H

#include <cerrno>
#include <string>

#include <sys/mman.h>
#include <unistd.h>


template< typename T >
T throw_errno_if_minus1( T x, const char* a, const char* b = 0 )
{ return throw_errno_if_eq( x, (T)(-1), a, b ) ; }

template< typename T >
T throw_errno_if_null( T x, const char* a, const char* b = 0 )
{ return throw_errno_if_eq( x, (T)0, a, b ) ; }

template< typename T >
T throw_errno_if_eq( T x, T y, const char* a, const char* b = 0 )
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
	return x ;
}

//! \brief near drop-in for write(2)
//! Would you believe it, write(2) sometimes decides to write less than
//! what was requested.  No idea why, but a loop around it solves it.  
inline void mywrite( int fd, const void* buf, size_t count, const char* msg = 0 )
{
	while( count > 0 ) 
	{
		ssize_t w = write( fd, buf, count ) ;
		if( w == -1 && errno == EINTR ) w = 0 ;
		throw_errno_if_minus1( w, "writing", msg ) ;
		count -= w ;
		buf = (const char*)buf + w ;
	}
}

#ifndef _BSD_SOURCE
// fake for systems that don't provide madvise()
inline int madvise( void*, size_t, int ) { return 0 ; }
static const int MADV_NORMAL = 0 ;
static const int MADV_SEQUENTIAL = 0 ;
static const int MADV_WILLNEED = 0 ;
#endif

int main( int argc, const char * argv[] ) ;

#endif
