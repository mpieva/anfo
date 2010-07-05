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

#ifndef INCLUDED_UTIL_H
#define INCLUDED_UTIL_H

#include <cerrno>
#include <cstring>
#include <cstdlib>
#include <deque>
#include <fcntl.h>
#include <iosfwd>
#include <sstream>
#include <string>
#include <unistd.h>

#include <sys/types.h>
#include <sys/stat.h>

struct Exception { virtual void print_to( std::ostream& ) const = 0 ;
                   virtual ~Exception() {} } ;

inline std::ostream& operator << ( std::ostream& s, const Exception& e ) { e.print_to( s ) ; return s ; }

inline void throw_errno( const char* a, const char* b = 0 )
{
	std::string msg ;
	msg.append( strerror( errno ) ).append( " while " ).append( a ) ;
	if( b ) msg.append( " " ).append( b ) ;
	throw msg ;
}

template< typename T >
T throw_errno_if_minus1( T x, const char* a, const char* b = 0 )
{ return throw_errno_if_eq( x, (T)(-1), a, b ) ; }

template< typename T >
T throw_errno_if_null( T x, const char* a, const char* b = 0 )
{ return throw_errno_if_eq( x, (T)0, a, b ) ; }

template< typename T >
T throw_errno_if_eq( T x, T y, const char* a, const char* b = 0 )
{
	if( x == y ) throw_errno( a, b ) ;
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

inline void throw_strerror_if_not_null( int x, const char* a, const char* b = 0 )
{
	if( x != 0 )
	{
		std::stringstream msg ;
		msg << strerror( x ) << " while " << a ;
		if( b ) msg << ' '<< b ;
		throw msg.str() ;
	}
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

//! \brief near drop-in for read(2)
inline void myread( int fd, void* buf, size_t count, const char* msg = 0 )
{
	while( count > 0 ) 
	{
		ssize_t w = read( fd, buf, count ) ;
		if( w == -1 && errno == EINTR ) w = 0 ;
		throw_errno_if_minus1( w, "reading", msg ) ;
		count -= w ;
		buf = (char*)buf + w ;
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
//! \return opened file descriptor

inline int path_open(
		const std::string& name, const std::string& ext,
		const char *env_var )
{
	int fd = open2( name, ext ) ;
	if( env_var ) {
		if( const char *p = getenv(env_var) ) {
			if( fd == -1 && *p ) {
				std::string path = p ;

				size_t end, start = 0 ;
				for( ; fd == -1 && (end = path.find( ':', start )) != std::string::npos ; start = end+1 )
					fd = open2( path.substr( start, end-start ) + "/" + name, ext ) ;
				if( fd == -1 ) fd = open2( path.substr( start ) + "/" + name, ext ) ;
			}
		}
	}
	return throw_errno_if_minus1( fd, "opening", name.c_str() ) ;
}

template< typename T > struct delete_ptr { void operator()( T* p ) const { delete p ; } } ;

int wrap_main( int argc, const char * argv[], int(*)(int,const char*[]) ) ;
#define WRAPPED_MAIN \
	int main_(int,const char*[]) ; \
	int main(int argc, const char*  argv[]) { return wrap_main( argc, argv, main_ ) ; } \
	int main_(int argc, const char* argv[])

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

//! \brief creates a temporary file 
//! The name is constructed from a directory taken from the environment
//! variables \c ANFO_TEMP, \c TMPDIR, \c TEMP, or \c TMP, whichever
//! first one is set, and a unique file name.  The file is opened and
//! then unlinked, and both the name and the file descriptor are
//! returned.  Note that the name can only serve informational purposes,
//! since it's already gone by the time this function returns.
//! \param name contains the file name on return
//! \return the file descriptor
int mktempfile( std::string* name = 0 ) ;

//! \brief here we send progress notification and log information
//! Output goes to the controlling terminal (/dev/tty) if possible,
//! progress notifications are displayed in the last line and everything
//! is updating sensibly.  Might grow into a logging system, too.
//! 
//! XXX this isn't thread safe.  should be changed...
class Console 
{
	public: 
		enum Loglevel { debug, info, notice, warning, error, critical } ;
		Loglevel loglevel ;

	private:
		int fd_ ;
		int next_ ;
		typedef std::deque< std::pair< int, std::string > > Chans ;
		Chans chans_ ;

	public:
		Console() : loglevel(warning), fd_( open( "/dev/tty", O_WRONLY ) ), next_(0) {}
		~Console() {
            if( fd_ >= 0 ) {
                if( !chans_.empty() ) {
                    chans_.clear() ;
                    update() ; 
                }
                close( fd_ ) ;
            }
        }

		int alloc_chan() { return ++next_ ; }
		void free_chan( int c ) ;
		void progress( int c, Loglevel l, const std::string& s ) ;
		void update() ; 

		void output( Loglevel, const std::string& ) ;

		void set_quiet() { loglevel = error ; }
		void more_verbose() { loglevel = loglevel > debug ? Loglevel(loglevel-1) : debug ; }
} ;

extern Console console ;
extern std::string program_name ;

template <typename T> void perr(const T& e)
{
	std::stringstream ss ;
	ss << program_name << "[" << getpid() << "]: " << e ;
	console.output( Console::error, ss.str() ) ; 
}

//! \brief a channel for progress notification
class Chan 
{
	private:
		int n_ ;

	public:
		Chan() : n_( console.alloc_chan() ) {}
		~Chan() { close() ; }
		void close() { console.free_chan( n_ ) ; }
		void operator()( Console::Loglevel l, const std::string& s ) { console.progress( n_, l, s ) ; }
} ;

template< typename T > class Holder
{
	private:
		T* s_ ;

	public:
		Holder<T>( const Holder<T>& rhs ) : s_( rhs.operator->() ) { if( s_ ) ++s_->refcount_ ; }
		template< typename U > Holder<T>( const Holder<U>& rhs ) : s_( rhs.operator->() ) { if( s_ ) ++s_->refcount_ ; }
		template< typename U > Holder<T>( U *s ) : s_( s ) { if( s_ ) ++s_->refcount_ ; }
		Holder<T>() : s_(0) {}
		~Holder<T>() { if( s_ && --s_->refcount_ == 0 ) T::cleanup( s_ ) ; }

		template< typename U > Holder<T>& operator = ( U *s ) 
		{ Holder<T>( s ).swap( *this ) ; return *this ; }

		Holder<T>& operator = ( const Holder<T>& s ) 
		{ Holder<T>( s ).swap( *this ) ; return *this ; }

		template< typename U > Holder<T>& operator = ( const Holder<U>& s ) 
		{ Holder<T>( s ).swap( *this ) ; return *this ; }

		void swap( Holder<T>& s ) { std::swap( s_, s.s_ ) ; }

		T* operator -> () const { return s_ ; }
		T& operator * () const { return *s_ ; }

		operator const void* () const { return s_ ; }
} ;

namespace std { template < typename T > void swap( Holder<T>& a, Holder<T>& b ) { a.swap( b ) ; } 

} ;

#endif
