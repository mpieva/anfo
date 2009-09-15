#include "config.h"

//! \page anfo-shell ANFO-tool with scripting language.
//! The idea: use the Guile interpreter as a scripting engine for
//! golfing ANFO files.  Expose the stream infrastruture as a Scheme
//! primitive along with suitable functions.  More scripting might be
//! possible, it could make sense to access the protobuf messages
//! through the built-in reflection interface---think 'customizable
//! filters'.

//! \todo Check what happens to string we pull out of SCM values.  They
//! should be copied and properly freed.

//! \todo Deal with errors: C++ exceptions need to be caught, formatted
//!       and reflected back into Scheme.

//! \todo Include many, many more streams...

#include <libguile.h>

#include "stream.h"

using namespace streams ;

static scm_t_bits stream_tag ;
     
extern "C" size_t free_stream( SCM stream_smob )
{
	delete (Stream*) SCM_SMOB_DATA( stream_smob ) ;
	return 0 ;
}

extern "C" SCM s_make_input_stream( SCM name )
{
	char *n = scm_to_locale_string( name ) ;
	Stream *s = make_input_stream( n, false, 33 ) ;
	SCM smob ;
	SCM_NEWSMOB( smob, stream_tag, s ) ;
	return smob ;
}

extern "C" SCM s_make_output_stream( SCM name )
{
	char *n = scm_to_locale_string( name ) ;
	Stream *s = new AnfoWriter( n ) ;
	SCM smob ;
	SCM_NEWSMOB( smob, stream_tag, s ) ;
	return smob ;
}

extern "C" SCM s_transfer( SCM input, SCM output ) 
{
	scm_assert_smob_type( stream_tag, input ) ;
	scm_assert_smob_type( stream_tag, output ) ;

	Stream* in  = (Stream*)SCM_SMOB_DATA(  input ) ;
	Stream *out = (Stream*)SCM_SMOB_DATA( output ) ;

	out->put_header( in->fetch_header() ) ;
	SCM_TICK;
	while( in->get_state() == Stream::have_output && out->get_state() == Stream::need_input )
	{
		out->put_result( in->fetch_result() ) ;
		SCM_TICK;
	}
	out->put_footer( in->fetch_footer() ) ;
	return SCM_BOOL_T ;
}

typedef SCM (*FCN)() ;

extern "C" void inner_main( void *closure, int argc, char **argv )
{
	console.more_verbose() ;
	console.more_verbose() ;
	stream_tag = scm_make_smob_type ("stream", sizeof (Stream*) ) ;
	scm_c_define_gsubr( "read-file", 1, 0, 0, (FCN)s_make_input_stream ) ;
	scm_c_define_gsubr( "write-file", 1, 0, 0, (FCN)s_make_output_stream ) ;
	scm_c_define_gsubr( "transfer", 2, 0, 0, (FCN)s_transfer ) ;

	scm_set_smob_free( stream_tag, free_stream ) ;
	scm_shell( argc, argv ) ;
}

int main_( int argc, const char** argv )
{
	scm_boot_guile( argc, (char**)argv, inner_main, 0 ) ;
	return 0 ;
}
