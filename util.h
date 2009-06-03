#ifndef INCLUDED_UTIL_H
#define INCLUDED_UTIL_H

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <cerrno>
#include <iosfwd>
#include <sstream>
#include <string>

#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>

#if HAVE_FCNTL_H
#include <fcntl.h>
#endif

#if HAVE_UNISTD_H
#include <unistd.h>
#endif

struct Exception { virtual void print_to( std::ostream& ) const = 0 ;
                   virtual ~Exception() {} } ;

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
		msg.append( strerror( errno ) ).append( " while " ).append( a ) ;
		if( b ) msg.append( " " ).append( b ) ;
		throw msg ;
	}
	return x ;
}

template< typename T >
T throw_if_negative( T x, const char* a, const char* b = 0 )
{
	if( x < 0 )
	{
		std::stringstream msg ;
		msg << x << " while " << a ;
		if( b ) msg << ' '<< b ;
		throw msg.str() ;
	}
	return x ;
}

template< typename T >
T throw_if_not_null( T x, const char* a, const char* b = 0 )
{
	if( x != 0 )
	{
		std::stringstream msg ;
		msg << x << " while " << a ;
		if( b ) msg << ' '<< b ;
		throw msg.str() ;
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


namespace {
	inline int open1( const std::string& s ) {
		int fd = open( s.c_str(), O_RDONLY ) ;
		if( errno != ENOENT ) throw_errno_if_minus1( fd, "opening", s.c_str() ) ;
		return fd ;
	}

	inline int open2( const std::string& s, const std::string& ext ) {
		int fd = open1( s ) ;
		if( fd == -1 ) fd = open1( s + "." + ext ) ;
		return fd ;
	}
} ;

//! \brief searches a file in the path, opens it for reading
//! This is equivalent to calling open( name, O_RDONLY ), but if the
//! file is not there, it additionally tried to find it with an
//! extension, then a supplied environmenmt variable is followed, then
//! an explicit list of paths.
//!
//! \param name name of the file
//! \param ext optional extension (without the dot)
//! \param env_var optional name of environment variable in same format
//!                as PATH
//! \param begin start iterator of list of directories
//! \param end end iterator of list of directories
//! \return opened file descriptor

template <typename Iter> int path_open(
		const std::string& name, const std::string& ext,
		const char *env_var, Iter begin, Iter end )
{
	int fd = open2( name, ext ) ;
	if( env_var ) {
		if( const char *p = getenv(env_var) ) {
			if( fd == -1 && *p ) {
				std::string path = p ;

				size_t end, start = 0 ;
				for( ; fd == -1 && (end = path.find( ':', start )) != std::string::npos ; start = end+1 )
					fd = open2( path.substr( start, end ) + "/" + name, ext ) ;
				if( fd == -1 ) fd = open2( path.substr( start ) + "/" + name, ext ) ;
			}
		}
	}
	for( ; fd == -1 && begin != end ; ++begin )
		fd = open2( *begin + ("/" + name), ext ) ;
	return throw_errno_if_minus1( fd, "opening", name.c_str() ) ;
}

template< typename T > struct delete_ptr { void operator()( T* p ) const { delete p ; } } ;

#ifndef _BSD_SOURCE
// fake for systems that don't provide madvise()
inline int madvise( void*, size_t, int ) { return 0 ; }
static const int MADV_NORMAL = 0 ;
static const int MADV_SEQUENTIAL = 0 ;
static const int MADV_WILLNEED = 0 ;
#endif

int main( int argc, const char * argv[] ) ;

/*! \brief tries to set proc title for ps
 *
 * This is borderline malpractice: since Linux is missing the
 * appropriate API, we directly overwrite argv[0].  This \e will trash
 * argv, so don't try to access that after setting a title, and this \e
 * may trash other things, should my assumption that argv is terminated
 * by "\0\0" turn out wrong.
 *
 * Ignoring the above, it is quite helpful, though...
 * \param title new program title to be displayed
 */
void set_proc_title( const char *title ) ;

//! \brief global flag to signal shutdown
//! Some long running algoritms will periodically check this variable
//! and abort cleanly if it is not null.  Use this to cause an early
//! shutdown, e.g. on reception of a signal.
extern volatile int exit_with ;

enum PacketTag { packet_config, packet_read, packet_result, packet_quit } ;


#endif
