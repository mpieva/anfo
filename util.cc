#include <util.h>
#include <iostream>

int main_( int argc, const char * const argv[] ) ;
int main( int argc, const char * const argv[] )
{
	try {
		return main_( argc, argv ) ;
	}
	catch( const std::string& e ) { std::cerr << e << std::endl ; }
	catch( const char *e ) { std::cerr << e << std::endl ; }
	catch( const std::exception& e ) { std::cerr << e.what() << std::endl ; }
	catch( ... ) { std::cerr << "Oh noes!" << std::endl ; }
	return 1 ;
}

