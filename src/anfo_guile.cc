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

//! \page anfo-guile Guile extension to control ANFO.   
//! The idea: use the Guile interpreter as a scripting engine for
//! golfing ANFO files.  To keep the interface manageable, we expose a
//! single function (run-anfo) that gets a description(!) of the stream
//! tree as argument and interprets it accordingly.  Scheme can be used
//! to compute the tree.

//! \todo Check what happens to string we pull out of SCM values.  They
//!       should be copied and freed properly.

//! \todo Deal with errors: C++ exceptions need to be caught, formatted
//!       and reflected back into Scheme.

//! \todo Include many, many more streams...
//!
//! \todo Deal with file descriptors/ports in a regular way; e.g.
//!       strings are file names, ints are file descriptors, #f is
//!       stdout or stdin, ports have their fd extracted, other scheme
//!       ports are wrapped (and garbage collected or something?).

#include <libguile.h>

#include "misc_streams.h"
#include "output_streams.h"
#include "stream.h"
#include "ducttape.h"

#include <strings.h>

using namespace streams ;

namespace {

typedef SCM (*FCN)() ;

const char *mk_str( SCM s, deque<char*> &stab )
{
	stab.push_back( scm_to_locale_string( s ) ) ;
	return stab.back() ;
}

extern "C" SCM scm_verbosity( SCM v )
{
	if( scm_is_integer( v ) ) console.loglevel = (Console::Loglevel)scm_to_int( v ) ;
	else if( !scm_is_symbol( v ) ) console.set_quiet() ;
	else {
		if( scm_is_eq( scm_from_locale_symbol(    "debug" ), v ) ) console.loglevel = Console::debug ;
		if( scm_is_eq( scm_from_locale_symbol(     "info" ), v ) ) console.loglevel = Console::info ;
		if( scm_is_eq( scm_from_locale_symbol(   "notice" ), v ) ) console.loglevel = Console::notice ;
		if( scm_is_eq( scm_from_locale_symbol(  "warning" ), v ) ) console.loglevel = Console::warning ;
		if( scm_is_eq( scm_from_locale_symbol(    "error" ), v ) ) console.loglevel = Console::error ;
		if( scm_is_eq( scm_from_locale_symbol( "critical" ), v ) ) console.loglevel = Console::critical ;
	}
	return SCM_BOOL_T ;
}

/*
scm_t_bits stream_tag ;
     
extern "C" size_t free_stream( SCM stream_smob )
{
	// XXX delete (Stream*) SCM_SMOB_DATA( stream_smob ) ;
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
		operator string () const { return s_ ; }
		const char* alt( const char* a ) const { return s_ ? s_ : a ; }
} ;

SCM mk_str( Stream* s )
{
	SCM smob ;
	SCM_NEWSMOB( smob, stream_tag, s ) ;
	return smob ;
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

extern "C" SCM scm_write_text( SCM fn )
{
	return mk_str( scm_is_string( fn )
			? new TextWriter( scm_to_str( fn ) )
			: new TextWriter( 1 ) ) ;
}

extern "C" SCM scm_write_sam( SCM fn, SCM genome )
{
	return mk_str( scm_is_string( fn )
			? new SamWriter( scm_to_str( fn ), scm_to_str( genome ) )
			: new SamWriter( cout.rdbuf(), scm_to_str( genome ) ) ) ;
}

extern "C" SCM scm_write_glz( SCM fn )
{
	return mk_str( scm_is_string( fn )
			? new GlzWriter( scm_to_str( fn ) )
			: new GlzWriter( 1 ) ) ;
}

extern "C" SCM scm_write_3aln( SCM fn ) 
{
	return mk_str( scm_is_string( fn )
			? new ThreeAlnWriter( scm_to_str( fn ) )
			: new ThreeAlnWriter( cout.rdbuf() ) ) ;
}

extern "C" SCM scm_write_fasta( SCM fn, SCM genome, SCM context )
{
	scm_to_str g( genome ) ;
	int c = scm_is_integer( context ) ? scm_to_int( context ) : 0 ;
	return mk_str( scm_is_string( fn ) 
			? new FastaAlnWriter( scm_to_str( fn ), g, c )
			: new FastaAlnWriter( cout.rdbuf(), g, c ) ) ;
} 

extern "C" SCM scm_write_fastq( SCM fn )
{
	return mk_str( scm_is_string( fn )
			? new FastqWriter( scm_to_str( fn ) )
			: new FastqWriter( cout.rdbuf() ) ) ;
}

extern "C" SCM scm_write_table( SCM fn, SCM genome )
{ 
	return mk_str( scm_is_string( fn )
			? new TableWriter( scm_to_str( fn ), scm_to_str( genome ) )
			: new TableWriter( cout.rdbuf(), scm_to_str( genome ) ) ) ;
}

extern "C" SCM scm_duct_tape( SCM genome, SCM name )
{ return mk_str( new DuctTaper( scm_to_str( genome ), scm_to_str( name ) ) ) ; }

extern "C" SCM scm_write_stats( SCM fn, SCM genome )
{ return mk_str( new StatStream( scm_to_str( fn ), scm_to_str( genome ) ) ) ; }

SCM scm_make_merger( SCM args, StreamBundle* b ) 
{
	for( ; scm_is_pair( args ) ; args = scm_cdr( args ) )
	{
		scm_assert_smob_type( stream_tag, scm_car( args ) ) ;
		Stream* s  = (Stream*)SCM_SMOB_DATA( scm_car( args ) ) ;
		b->add_stream( s ) ;
	}
	assert( scm_is_null( args ) ) ;
	return mk_str( b ) ;
}

extern "C" SCM scm_merge( SCM args )
{ return scm_make_merger( args, new MergeStream() ) ; }

extern "C" SCM scm_join( SCM args )
{ return scm_make_merger( args, new BestHitStream() ) ; }

extern "C" SCM scm_concat( SCM args )
{ return scm_make_merger( args, new ConcatStream() ) ; }

extern "C" SCM scm_mega_merge( SCM args )
{ return scm_make_merger( args, new MegaMergeStream() ) ; }

extern "C" SCM scm_sort_by_pos( SCM mem, SCM handles, SCM genome )
{
 return mk_str( new SortingStream<by_genome_coordinate>( 
			 scm_to_uint( mem ) * 1024U * 1024U, scm_to_uint( handles ),
			 by_genome_coordinate( scm_to_str( genome ) ) ) ) ;
}

extern "C" SCM scm_sort_by_name( SCM mem, SCM handles )
{
 return mk_str( new SortingStream<by_seqid>( 
			 scm_to_uint( mem ) * 1024U * 1024U, scm_to_uint( handles ) ) ) ;
}

extern "C" SCM scm_filter_by_length( SCM l )
{ return mk_str( new LengthFilter( scm_to_int( l ) ) ) ; }

extern "C" SCM scm_filter_by_score( SCM s, SCM l, SCM g )
{ return mk_str( new ScoreFilter( scm_to_double( s ), scm_to_double( l ), scm_to_str( g ) ) ) ; }

extern "C" SCM scm_filter_by_mapq( SCM q, SCM g )
{ return mk_str( new MapqFilter( scm_to_str( g ), scm_to_int( q ) ) ) ; }

extern "C" SCM scm_ensure_hit( SCM g, SCM s )
{ return mk_str( new HitFilter( scm_to_str( g ), scm_to_str( s ) ) ) ; }

extern "C" SCM scm_delete_hit( SCM g, SCM s )
{ return mk_str( new IgnoreHit( scm_to_str( g ), scm_to_str( s ) ) ) ; } // XXX

extern "C" SCM scm_ensure_multi( SCM m )
{ return mk_str( new MultiFilter( scm_is_integer( m ) ? scm_to_int( m ) : 2 ) ) ; }

extern "C" SCM scm_subsample( SCM r )
{ return mk_str( new Subsample( scm_to_double( r ) ) ) ; }

extern "C" SCM scm_edit_header( SCM e )
{ return mk_str( new RepairHeaderStream( scm_to_str( e ) ) ) ; } // XXX

extern "C" SCM scm_rmdup( SCM s, SCM i, SCM q )
{ return mk_str( new RmdupStream( scm_to_double( s ), scm_to_double( i ), scm_to_int( q ) ) );}

extern "C" SCM scm_regions_only( SCM fn )
{ return mk_str( new InsideRegion( scm_to_str( fn ) ) ) ; }

extern "C" SCM scm_not_regions( SCM fn )
{ return mk_str( new OutsideRegion( scm_to_str( fn ) ) ) ; }
*/


//!	- a string: read the file
//!	- an integer: read from the fd
//!	- a scheme port: read from the port (?) 

Stream* stream_from_literal( SCM lit, deque<char*> &stab, bool solexa = false, int origin = 33 )
{
	if( scm_is_string( lit ) )
	{
		return make_input_stream( mk_str( lit, stab ), solexa, origin ) ;
	}
	else if( scm_is_integer( lit ) )
	{
		return make_input_stream( scm_to_int( lit ), "<pipe>", solexa, origin ) ;
	}
	// handle port?  How?
	else
		throw "cannot handle literal" ;
}

Stream* stream_from_sym( SCM sym, SCM args, deque<char*> &stab )
{
	SCM arg[5] ;
	for( int i = 0 ; i != 5 ; ++i )
	{
		if( scm_is_pair( args ) )
		{
			arg[i] = scm_car( args ) ;
			args = scm_cdr( args ) ;
		}
		else arg[i] = SCM_BOOL_F ;
	}

	if( scm_is_eq( sym, scm_from_locale_symbol( "read-file" ) ) )
	{
		cerr << "reading anything: " << scm_to_locale_string( scm_symbol_to_string( sym ) ) << ' ' << scm_to_locale_string( arg[0] ) << endl ;
		return stream_from_literal( arg[0], stab, scm_is_true( arg[1] ),
				scm_is_integer( arg[2] ) ? scm_to_int( arg[2] ) : 33 ) ;
	}
	else if( scm_is_eq( sym, scm_from_locale_symbol( "write-anfo" ) ) )
	{
		cerr << "writing native" << endl ;
		bool exp = scm_is_integer( arg[1] ) ? scm_to_int( arg[1] ) > 50 : false ;
		SCM name = arg[0] ;
		return scm_is_false( name )   ? new AnfoWriter( 1, "<stdout>", exp ) :
			   scm_is_string( name )  ? new AnfoWriter( mk_str( name, stab ), exp ) :
			   scm_is_integer( name ) ? new AnfoWriter( scm_to_int( name ), "<pipe>", exp ) :
			   throw "cannot handle weird file name" ;
	}
	else if( scm_is_eq( sym, scm_from_locale_symbol( "filter-qual" ) ) )
	{
		SCM qual = arg[0] ;
		return new QualFilter( scm_is_integer( qual ) ? scm_to_int( qual ) : 30 ) ;
	}
	else throw "I don't understand" ;
}

Stream* stream_from_list( SCM list, deque<char*> &stab ) ;

//!	Possible arguments:
//!	- a keyword: stream without arguments
//!	- a list starting with a keyword: stream with arguments
//!	- any other list: compose the elements
//!	- just about everything else: read a file
Stream* stream_from_args( SCM args, deque<char*> &stab )
{
	if( scm_is_symbol( args ) )
	{
		cerr << "symbol found" << endl ;
		return stream_from_sym( args, SCM_EOL, stab ) ;
	}
	else if( scm_is_pair( args ) && scm_is_symbol( scm_car( args ) ) )
	{
		cerr << "list w/ symbol found" << endl ;
		return stream_from_sym( scm_car( args ), scm_cdr( args ), stab ) ;
	}
	else if( scm_is_pair( args ) )
	{
		cerr << "list found" << endl ;
		return stream_from_list( args, stab ) ;
	}
	else 
	{
		cerr << "literal found" << endl ;
		return stream_from_literal( args, stab ) ;
	}
}

Stream* stream_from_list( SCM list, deque<char*> &stab )
{
	Compose *s = new Compose ;
	for( ; scm_is_pair( list ) ; list = scm_cdr( list ) )
	{
		cerr << "recurse for stream" << endl ;
		s->add_stream( stream_from_args( scm_car( list ), stab ) ) ;
	}
	return s ; 
}

//! \brief interprets a description of a streaming operation
//! What to do?  We effectively construct a single stream as composition
//! of whatever, then determine its status.  Something should result
//! from that...  we might even return a tree of status values!
//!
//!	
//!	If a file argument is needed:
//!	- #f: stdin or stdout
//!	- string: file name
//!	- int: file descriptor
//!	- scheme port: read or write to port (?)
//!
//! We need some memory management: strings are malloc()ed, we pass them
//! around, but finally need to free them.  We pass around a deque of
//! pointers.

extern "C" SCM scm_anfo_run( SCM args )
{
	try
	{
		// construct whatever stream is needed
		deque<char*> stringtab ;
		Stream *s = stream_from_args( args, stringtab ) ;

		// top level is probably of type 'Compose'; get it to calculate
		if( Compose *c = dynamic_cast<Compose*>( s ) ) c->update_status() ;

		// cleanup
		for_each( stringtab.begin(), stringtab.end(), free ) ;

		// extract some sort of result and return it?
		return SCM_BOOL_T ;
	} 
	catch( const Exception& e ) { 
		stringstream s ;
		s << e ;
		scm_throw( scm_from_locale_symbol( "anfo-error" ), scm_from_locale_string( s.str().c_str() ) ) ;
	}
	catch( const string& s ) { scm_throw( scm_from_locale_symbol( "anfo-error" ), scm_from_locale_string( s.c_str() ) ) ; }
	catch( const char* s ) { scm_throw( scm_from_locale_symbol( "anfo-error" ), scm_from_locale_string( s ) ) ; }
	catch( const exception& e ) { scm_throw( scm_from_locale_symbol( "anfo-error" ), scm_from_locale_string( e.what() ) ) ; }
	catch( ... ) { scm_throw( scm_from_locale_symbol( "anfo-error" ), SCM_BOOL_F ) ; }
	return SCM_BOOL_F ;
}


} ; // namespace

extern "C" void init_anfo_guile()
{
	// stream_tag = scm_make_smob_type ("stream", sizeof (Stream*) ) ;
	// scm_set_smob_free( stream_tag, free_stream ) ;

	scm_c_define_gsubr( "verbosity",       1, 0, 0, (FCN)scm_verbosity ) ;
	scm_c_define_gsubr( "anfo-run",        0, 0, 1, (FCN)scm_anfo_run ) ;

	/*
	scm_c_define_gsubr( "prim-sort-by-name",    2, 0, 0, (FCN)scm_sort_by_name ) ;
	scm_c_define_gsubr( "prim-sort-by-pos",     3, 0, 0, (FCN)scm_sort_by_pos ) ;
	scm_c_define_gsubr(     "filter-by-length", 1, 0, 0, (FCN)scm_filter_by_length ) ;
	scm_c_define_gsubr( "prim-filter-by-score", 3, 0, 0, (FCN)scm_filter_by_score ) ;
	scm_c_define_gsubr( "prim-filter-by-mapq",  2, 0, 0, (FCN)scm_filter_by_mapq ) ;
	scm_c_define_gsubr( "prim-ensure-hit",      2, 0, 0, (FCN)scm_ensure_hit ) ;
	scm_c_define_gsubr( "prim-delete-hit",      2, 0, 0, (FCN)scm_delete_hit ) ;
	scm_c_define_gsubr(      "filter-by-qual",  1, 0, 0, (FCN)scm_filter_by_qual ) ;
	scm_c_define_gsubr(  "ensure-multiplicity", 1, 0, 0, (FCN)scm_ensure_multi ) ;
	scm_c_define_gsubr(      "subsample",       1, 0, 0, (FCN)scm_subsample ) ;
	scm_c_define_gsubr( "prim-edit-header",     1, 0, 0, (FCN)scm_edit_header ) ;
	scm_c_define_gsubr( "prim-rmdup",           2, 0, 0, (FCN)scm_rmdup ) ;
	scm_c_define_gsubr(      "regions-only",    1, 0, 0, (FCN)scm_regions_only ) ;
	scm_c_define_gsubr(      "not-regions",     1, 0, 0, (FCN)scm_not_regions ) ;

	scm_c_define_gsubr( "prim-duct-tape",       3, 0, 0, (FCN)scm_duct_tape ) ;

	scm_c_define_gsubr(      "merge-streams",   0, 0, 1, (FCN)scm_merge ) ;
	scm_c_define_gsubr(      "join-streams",    0, 0, 1, (FCN)scm_join ) ;
	scm_c_define_gsubr(      "concat-streams",  0, 0, 1, (FCN)scm_concat ) ;
	scm_c_define_gsubr(   "mega-merge-streams", 0, 0, 1, (FCN)scm_mega_merge ) ;

	scm_c_define_gsubr( "prim-read-file",       3, 0, 0, (FCN)scm_make_input_stream ) ;

	scm_c_define_gsubr( "prim-write-anfo-file", 3, 0, 0, (FCN)scm_make_output_stream ) ;
	scm_c_define_gsubr(      "write-text",      1, 0, 0, (FCN)scm_write_text ) ;
	scm_c_define_gsubr( "prim-write-sam",       2, 0, 0, (FCN)scm_write_sam ) ;
	scm_c_define_gsubr(      "write-glz",       1, 0, 0, (FCN)scm_write_glz ) ;
	scm_c_define_gsubr(      "write-3aln",      1, 0, 0, (FCN)scm_write_3aln ) ;
	scm_c_define_gsubr( "prim-write-fasta",     3, 0, 0, (FCN)scm_write_fasta ) ;
	scm_c_define_gsubr(      "write-fastq",     1, 0, 0, (FCN)scm_write_fastq ) ;
	scm_c_define_gsubr( "prim-write-table",     2, 0, 0, (FCN)scm_write_table ) ;
	scm_c_define_gsubr( "prim-write-stats",     2, 0, 0, (FCN)scm_write_stats ) ;

	scm_c_define_gsubr(      "transfer",        2, 0, 0, (FCN)scm_transfer ) ;
	scm_c_define_gsubr(      "chain",           0, 0, 1, (FCN)scm_chain ) ;
	*/
}

