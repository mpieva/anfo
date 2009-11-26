#include "../config.h"

#include "index.h"
#include "stream.h"
#include "util.h"

#include <elk/scheme.h>

#include <iostream>

using namespace std ;
using namespace streams ;

// generic interfacing to ELK: primitive procedures
typedef Object (*P)() ;

// ugly hack: this entry point is neither declared nor documented, but I
// don't see how else to enter the REPL
extern "C" Object General_Load( Object, Object ) ;


// the wrapped stream type: construction, destruction, predicates
static int t_stream ;

struct StreamWrapper {
	Object o_ ;
	StreamHolder h_ ;

	StreamWrapper( StreamHolder h ) : h_(h) {}
} ;

extern "C" Object terminate_stream( Object o )
{
	((StreamWrapper*)POINTER(o))->~StreamWrapper() ;
	Deregister_Object( o ) ;
	return Void ;
}

Object wrap_stream( StreamHolder h )
{
	Object o = Alloc_Object( sizeof( StreamWrapper ), t_stream, 0 ) ;
	new( POINTER(o) ) StreamWrapper( h ) ;
	Register_Object( o, const_cast<char*>("anfostreams") , terminate_stream, 0 ) ;
	return o ;
}

StreamHolder obj_to_stream( Object o ) { return TYPE( o ) != t_stream ? StreamHolder() : ((StreamWrapper*)POINTER(o))->h_ ; }

extern "C" {
Object p_is_stream( Object d ) { return TYPE(d) == t_stream ? True : False; }
int compare_stream_wrappers( Object a, Object b ) { return obj_to_stream(a) == obj_to_stream(b) ; }

int print_stream_wrapper( Object p, Object o, int, int, int )
{ Printf( p, "#[anfo-stream %p]", (const void*)((StreamWrapper*)POINTER(o))->h_ ) ; return 0 ; }

} // extern C

// misc. ANFO primitves
extern "C" Object p_set_verbosity( Object v ) { console.loglevel = (Console::Loglevel)Get_Integer( v ) ; return Void ; }
extern "C" Object p_use_mmap( Object v ) { Metagenome::nommap = !Truep( v ) ; return Void ; }

// ANFO stream constructors

int main_( int argc, const char**argv )
{
	Elk_Init( argc, (char**)argv, 0, 0 ) ;

	t_stream = Define_Type( 0, "anfo-stream", 
        0, sizeof( StreamWrapper ),
		compare_stream_wrappers, 
		compare_stream_wrappers, 
		print_stream_wrapper, 0 ) ;

	Define_Primitive( (P)p_is_stream, "anfo-stream?", 1, 1, EVAL ) ;
	Define_Primitive( (P)p_set_verbosity, "set-verbosity", 1, 1, EVAL ) ;
	Define_Primitive( (P)p_use_mmap, "use-mmap", 1, 1, EVAL ) ;

	// XXX this is nonsense: REPL is entered by loading toplevel.scm;
	// there has to be a clean way to do that (without accessing
	// General_Load)
	Object file = Make_String( "toplevel.scm", 12 ) ;
	General_Load( file, The_Environment ) ;				// enters REPL
	return 0 ;
}

