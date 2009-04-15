#include "util.h"
#include <iostream>

volatile int exit_with = 0 ;

extern int main_( int argc, const char * argv[] ) ;
int main( int argc, const char * argv[] )
{
	try { return main_( argc, argv ) ; }
	catch( const std::string& e ) { std::cerr << "\r\33[K" << e << std::endl ; }
	catch( const char *e ) { std::cerr << "\r\33[K" << e << std::endl ; }
	catch( char *e ) { std::cerr << "\r\33[K" << e << std::endl ; }
	catch( const Exception& e ) { std::cerr << "\r\33[K" ; e.print_to( std::cerr ) ; std::cerr << std::endl ; }
	catch( const std::exception& e ) { std::cerr << "\r\33[K" << e.what() << std::endl ; }
	catch( ... ) { std::cerr << "\r\33[KOh noes!" << std::endl ; }
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

