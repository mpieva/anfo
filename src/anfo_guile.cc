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

extern "C" SCM scm_chain( SCM args )
{
	Compose* c = new Compose ;
	for( ; scm_is_pair( args ) ; args = scm_cdr( args ) )
	{
		scm_assert_smob_type( stream_tag, scm_car( args ) ) ;
		Stream* s  = (Stream*)SCM_SMOB_DATA( scm_car( args ) ) ;
		c->add_stream( s ) ;
	}
	assert( scm_is_null( args ) ) ;
	return mk_str( c ) ;
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

extern "C" SCM scm_verbosity( SCM v )
{
	if( scm_is_integer( v ) ) console.loglevel = (Console::Loglevel)scm_to_int( v ) ;
	else if( !scm_is_symbol( v ) ) console.set_quiet() ;
	else {
		char *s = scm_to_locale_string( scm_symbol_to_string( v ) ) ;
		if( !strcasecmp( s, "debug" ) ) console.loglevel = Console::debug ;
		if( !strcasecmp( s, "info" ) ) console.loglevel = Console::info ;
		if( !strcasecmp( s, "notice" ) ) console.loglevel = Console::notice ;
		if( !strcasecmp( s, "warning" ) ) console.loglevel = Console::warning ;
		if( !strcasecmp( s, "error" ) ) console.loglevel = Console::error ;
		if( !strcasecmp( s, "critical" ) ) console.loglevel = Console::critical ;
		free( s ) ;
	}
	return SCM_BOOL_T ;
}

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

extern "C" SCM scm_filter_by_qual( SCM q )
{ return mk_str( new QualFilter( scm_to_int( q ) ) ) ; }

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


typedef SCM (*FCN)() ;

} ; // namespace

extern "C" void init_anfo_guile()
{
	stream_tag = scm_make_smob_type ("stream", sizeof (Stream*) ) ;
	scm_set_smob_free( stream_tag, free_stream ) ;

	scm_c_define_gsubr( "verbosity",       1, 0, 0, (FCN)scm_verbosity ) ;
	// scm_c_define_gsubr( "anfo-run

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
}

