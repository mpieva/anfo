#include "util.h"
#include <iostream>

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

