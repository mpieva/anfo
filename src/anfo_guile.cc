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

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

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

#include "misc_streams.h"
#include "output_streams.h"
#include "stream.h"
#include "ducttape.h"

using namespace streams ;

namespace {

scm_t_bits stream_tag ;
     
extern "C" size_t free_stream( SCM stream_smob )
{
	delete (Stream*) SCM_SMOB_DATA( stream_smob ) ;
	return 0 ;
}

class scm_to_str
{
	private:
		char *s_ ;

	public:
		scm_to_str( SCM s ) : s_( scm_is_string( s ) ? scm_to_locale_string( s ) : 0 ) {}
		~scm_to_str() { free( s_ ) ; }
		operator const char* () const { return s_ ; }
		const char* alt( const char* a ) const { return s_ ? s_ : a ; }
} ;

SCM mk_str( Stream* s )
{
	SCM smob ;
	SCM_NEWSMOB( smob, stream_tag, s ) ;
	return smob ;
}

extern "C" SCM scm_make_input_stream( SCM name, SCM sol_scale, SCM origin )
{
	return mk_str( make_input_stream(
				scm_to_str( name ), scm_is_true( sol_scale ), scm_to_int( origin ) ) ) ;
}

extern "C" SCM scm_make_output_stream( SCM name, SCM level )
{
	bool exp = scm_to_int( level ) > 50 ;
	return mk_str( scm_is_true( name )
		? new AnfoWriter( scm_to_str( name ), exp ) : new AnfoWriter( 1, "<stdout>", exp ) ) ;
}

extern "C" SCM scm_transfer( SCM input, SCM output ) 
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

extern "C" SCM scm_duct_tape( SCM genome, SCM name )
{ return mk_str( new DuctTaper( scm_to_str( genome ), scm_to_str( name ) ) ) ; }

extern "C" SCM scm_write_stats( SCM fn, SCM genome )
{ return mk_str( new StatStream( scm_to_str( fn ), scm_to_str( genome ) ) ) ; }

typedef SCM (*FCN)() ;

} ; // namespace

extern "C" void init_anfo_guile()
{
	console.more_verbose() ; // optional... should be an option?
	console.more_verbose() ;

	stream_tag = scm_make_smob_type ("stream", sizeof (Stream*) ) ;
	scm_set_smob_free( stream_tag, free_stream ) ;
	scm_c_define_gsubr( "prim-read-file", 3, 0, 0, (FCN)scm_make_input_stream ) ;
	scm_c_define_gsubr( "prim-write-anfo-file", 3, 0, 0, (FCN)scm_make_output_stream ) ;
	scm_c_define_gsubr( "prim-duct-tape", 3, 0, 0, (FCN)scm_duct_tape ) ;
	scm_c_define_gsubr( "prim-write-stats", 2, 0, 0, (FCN)scm_write_stats ) ;
	scm_c_define_gsubr( "transfer", 2, 0, 0, (FCN)scm_transfer ) ;
}

