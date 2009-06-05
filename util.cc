#include "util.h"

#include <cstdlib>
#include <cstring>
#include <iostream>

volatile int exit_with = 0 ;

namespace {
	std::string program_name ;
	template <typename T> void perr(const T& e) {
		std::cerr << "\r" << program_name << "[" << getpid() << "]: "
			      << e << "\33[K" << std::endl ;
	}
	std::ostream& operator << ( std::ostream& s, const Exception& e ) {
		e.print_to( s ) ; return s ;
	}
} ;

extern int main_( int argc, const char * argv[] ) ;
int main( int argc, const char * argv[] )
{
	program_name = argv[0] ;
	try { return main_( argc, argv ) ; }
	catch( const std::string& e ) { perr( e ) ; }
	catch( const char *e ) { perr( e ) ; }
	catch( char *e ) { perr( e ) ; }
	catch( const Exception& e ) { perr( e ) ; }
	catch( const std::exception& e ) { perr( e.what() ) ; }
	catch( ... ) { perr( "Oh noes!" ) ; }
	return 1 ;
}

void set_proc_title( const char *title ) 
{
	extern char* __progname_full ;
	extern char* __progname ;
	static char* pe = 0 ;
	static char* pa = 0 ;

	if( !pe ) {
		char* p = __progname ;
		pa = __progname_full ;
		for( pe = pa ; pe[0] || pe[1] ; ++pe ) ;
		while( *p && pa != pe ) *pa++ = *p++ ;
		if( pa != pe ) *pa++ = ':' ;
		if( pa != pe ) *pa++ = ' ' ;
	}

	char* pf = pa ;
	if( pf != pe ) {
		while( *title && pf != pe ) *pf++ = *title++ ;
		while( pf != pe ) *pf++ = 0 ;
	}
}

int mktempfile( std::string* name )
{
	const char *suffix = "/anfo_sort_XXXXXX" ;
	const char *base = getenv("ANFO_TEMP") ;
	if( !base ) base = getenv("TMPDIR") ;
	if( !base ) base = getenv("TEMP") ;
	if( !base ) base = getenv("TMP") ;
	if( !base ) base = "." ;

	char  n1[ strlen(base) + strlen(suffix) + 1 ] ;
	char *n2 = n1 ;
	while( *base ) *n2++ = *base++ ;
	while( *suffix ) *n2++ = *suffix++ ;
	*n2 = 0 ;
    int fd = throw_errno_if_minus1( mkstemp( n1 ), "making temp file" ) ;
	throw_errno_if_minus1( unlink( n1 ), "unlinking temp name" ) ;
	if( name ) *name = n1 ;
	return fd ;
}
