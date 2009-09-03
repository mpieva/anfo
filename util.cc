#include "config.h"
#include "util.h"

#include <cstdlib>
#include <cstring>
#include <iostream>

#include <sys/ioctl.h>

volatile int exit_with = 0 ;
std::string program_name ;
Console console ;

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

namespace {
	const char *describe[] = { "[debug] ", "[info] ", "[notice] ",
		"[warning] ", "[error] ", "[critical] " } ;
} ;

void Console::output( Loglevel l, const std::string& s ) 
{
	if( l < loglevel ) return ;
	if( fd_ < 0 ) { std::clog << s << std::endl ; }
	if( s.empty() ) return ;
	mywrite( fd_, "\r\e[K", 4 ) ;
	mywrite( fd_, describe[l], strlen( describe[l] ) ) ;
	mywrite( fd_, s.data(), s.size() ) ;
	if( s[ s.size()-1 ] != '\n' ) mywrite( fd_, "\n", 1 ) ;
	update() ;
}

void Console::update()
{
	if( fd_ < 0 ) return ;
	int width = 79 ;
	struct winsize ws;
    if( 0 == ioctl( fd_, TIOCGWINSZ, &ws ) ) width = ws.ws_col-1 ;

	std::string line = "\r\e[K" ;
	for( std::map<int, std::string>::const_reverse_iterator ch = chans_.rbegin() ;
		width >= 3 && ch != chans_.rend() ; ++ch )
	{
		line.push_back( '[' ) ;
		line.append( ch->second.substr( 0, width-2 ) ) ;
		line.push_back( ']' ) ;
		line.push_back( ' ' ) ;
		width -= ch->second.size()+3 ;
	}
	mywrite( fd_, line.data(), line.size() - (line[line.size()-1] == ' ') ) ;
}


