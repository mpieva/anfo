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

//! \page anfo-elk ELK extension to control ANFO.   
//! The idea: use the ELK interpreter as a scripting engine for
//! golfing ANFO files.
//!
//! \todo Deal with errors: C++ exceptions need to be caught, formatted
//!       and reflected back into Scheme.
//! \todo Deal with scheme ports and raw file descriptors somehow.

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include "ducttape.h"
#include "index.h"
#include "misc_streams.h"
#include "output_streams.h"
#include "stream.h"
#include "util.h"

#include <elk/scheme.h>
#include <google/protobuf/io/zero_copy_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>

#include <iostream>

using namespace std ;
using namespace streams ;
using namespace google::protobuf::io ;

// generic interfacing to ELK: primitive procedures
typedef Object (*P)() ;

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
	StreamWrapper *w = new( POINTER(o) ) StreamWrapper( h ) ;
	Register_Object( o, const_cast<char*>("anfostreams") , terminate_stream, 0 ) ;
	w->h_->get_state() ;
	return o ;
}

string object_to_string( Object o, const string& def = "" ) 
{ 
	if( TYPE(o) == T_Symbol ) o = SYMBOL(o)->name ;
	return TYPE(o) == T_String ? string( STRING(o)->data, STRING(o)->size ) : def ; 
}

pair< ZeroCopyOutputStream*, string > open_any_output_zc( Object o )
{
	switch( TYPE(o) )
	{
		case T_Symbol:
		case T_String:
			{
				string nm = object_to_string(o) ;
				FileOutputStream *s = new FileOutputStream( 
						throw_errno_if_minus1(
							open( nm.c_str(), O_WRONLY | O_CREAT ),
							"opening file" ) ) ;
				s->SetCloseOnDelete( true ) ;
				return make_pair( s, nm ) ;
			}
		case T_Fixnum:
			return make_pair( new FileOutputStream( Get_Exact_Integer( o ) ), "<pipe>" ) ;

		case T_Boolean:
			if( !Truep(o) ) return make_pair( new FileOutputStream( 1 ), "<stdout>" ) ;

		case T_Port: // needs support code
		default:
			throw "can't handle file argument" ;
	}
}

pair< std::ostream*, string > open_any_output_std( Object o )
{
	switch( TYPE(o) )
	{
		case T_Symbol:
		case T_String:
			return make_pair( new ofstream( object_to_string(o).c_str() ), object_to_string(o) ) ;

		case T_Boolean:
			if( !Truep(o) ) return make_pair( new ostream( cout.rdbuf() ), "<stdout>" ) ;

		case T_Fixnum: // needs support code
		case T_Port: // needs support code
		default:
			break ;
	}
	throw "can't handle file argument" ;
}

pair< std::istream*, string > open_any_input_std( Object o )
{
	switch( TYPE(o) )
	{
		case T_Symbol:
		case T_String:
			return make_pair( new ifstream( object_to_string(o).c_str() ), object_to_string(o) ) ;

		case T_Boolean:
			if( !Truep(o) ) return make_pair( new istream( cin.rdbuf() ), "<stdin>" ) ;

		case T_Fixnum: // needs support code
		case T_Port: // needs support code
		default:
			break ;
	}
	throw "can't handle file argument" ;
}

StreamHolder obj_to_stream( Object o, bool sol = false , int ori = 33 )
{
	if( TYPE(o) == t_stream ) return ((StreamWrapper*)POINTER(o))->h_ ;
	switch( TYPE(o) )
	{
		case T_Primitive:
		case T_Compound:
			return obj_to_stream( Funcall( o, Null, 0 ), sol, ori ) ;

		case T_Symbol:
		case T_String:
			return make_input_stream( object_to_string(o).c_str(), sol, ori ) ;

		case T_Boolean:
			if( !Truep(o) ) return make_input_stream( 0, "<stdin>", sol, ori ) ;

		case T_Fixnum: // needs support code
		case T_Port: // needs support code
			break ;
	}
	throw "can't handle file argument" ;
}

vector<string> obj_to_genomes( Object o )
{
	vector<string> r ;
	if( TYPE(o) == T_String || TYPE(o) == T_Symbol )
		r.push_back( object_to_string( o ) ) ;
	else
		for( ; TYPE(o) == T_Pair ; o = Cdr(o) ) 
			r.push_back( object_to_string( Car(o) ) ) ;
	return r ;
}

Object wrap_streams( StreamBundle *m_, int argc, Object *argv )
{
	Holder< StreamBundle > m( m_ ) ;
	for( Object *o = argv ; o != argv + argc ; ++o )
		m->add_stream( obj_to_stream( *o ) ) ;
	return wrap_stream( m ) ;
}


#define WRAP( fun, proto, args ) \
	Object wrap_##fun proto ; \
	Object fun proto { \
		try { return wrap_##fun args ; } \
		catch( const std::string& e ) { Primitive_Error( e.c_str() ) ; } \
		catch( const char *e ) { Primitive_Error( e ) ; } \
		catch( const Exception& e ) { stringstream ss ; ss << e ; Primitive_Error( ss.str().c_str() ) ; } \
		catch( const std::exception& e ) { Primitive_Error( e.what() ) ; } \
		catch( ... ) { Primitive_Error( "unhandled C++ exception" ) ; } \
		return Void ; \
	} \
	Object wrap_##fun proto 

extern "C" {

Object p_is_stream( Object d ) { return TYPE(d) == t_stream ? True : False; }
int compare_stream_wrappers( Object a, Object b ) { return obj_to_stream(a) == obj_to_stream(b) ; }

int print_stream_wrapper( Object p, Object o, int, int, int )
{ Printf( p, "#[anfo-stream %p]", (const void*)((StreamWrapper*)POINTER(o))->h_ ) ; return 0 ; }

// misc. ANFO primitves
SYMDESCR verbosity_syms[] = {
	{ (char*) "debug",    Console::debug },
	{ (char*) "info",     Console::info },
	{ (char*) "notice",   Console::notice },
	{ (char*) "warning",  Console::warning },
	{ (char*) "error",    Console::error },
	{ (char*) "critical", Console::critical },
	{ 0, 0 } } ;

WRAP( p_set_verbosity, ( Object v ), (v) )
{
	if( TYPE(v) == T_Symbol )
		console.loglevel = (Console::Loglevel) Symbols_To_Bits( v, 0, verbosity_syms ) ;
	else 
		console.loglevel = (Console::Loglevel) Get_Integer( v ) ;
	return Void ;
}

WRAP( p_use_mmap, ( Object v ), (v) ) { Metagenome::nommap = !Truep( v ) ; return Void ; }

// ANFO stream constructors

// Output
WRAP( p_write_native, ( Object f, Object c ), (f,c) ) { return wrap_stream( new ChunkedWriter( open_any_output_zc( f ), Get_Integer( c ) ) ) ; }
WRAP( p_write_text,   ( Object f ), (f) ) { return wrap_stream( new TextWriter( open_any_output_zc( f ) ) ) ; }
WRAP( p_write_sam,    ( Object f ), (f) ) { return wrap_stream( new SamWriter( open_any_output_std( f ) ) ) ; } 
WRAP( p_write_glz,    ( Object f ), (f) ) { return wrap_stream( new GlzWriter( open_any_output_zc( f ) ) ) ; }
WRAP( p_write_3aln,   ( Object f ), (f) ) { return wrap_stream( new ThreeAlnWriter( open_any_output_std( f ) ) ) ; }
WRAP( p_write_fastq,  ( Object f ), (f) ) { return wrap_stream( new FastqWriter( open_any_output_std( f ) ) ) ; }
WRAP( p_write_table,  ( Object f ), (f) ) { return wrap_stream( new TableWriter( open_any_output_std( f ) ) ) ; }
WRAP( p_write_fasta,  ( Object f ), (f) ) { return wrap_stream( new FastaAlnWriter( open_any_output_std( f ) ) ) ; }

// Processors
WRAP( p_duct_tape, ( Object n ), (n) ) { return wrap_stream( new DuctTaper( object_to_string( n, "contig" ) ) ) ; }
WRAP( p_add_alns, ( Object c ), (c) ) { return wrap_stream( new GenTextAlignment( Get_Integer( c ) ) ) ; }
WRAP( p_rmdup, ( Object s, Object i, Object q ), (s,i,q) ) { return wrap_stream( new RmdupStream( Get_Double(s), Get_Double(i), Get_Integer(q) ) ) ; }
WRAP( p_write_stats, ( Object f ), (f) ) { return wrap_stream( new StatStream( object_to_string( f ) ) ) ; }

// Filters
WRAP( p_sort_by_pos, ( Object mem, Object handles, Object genomes ), (mem,handles,genomes) )
{ return wrap_stream( new SortingStream<by_genome_coordinate>(
			Get_Integer( mem ) * 1024U * 1024U, Get_Integer( handles ),
			by_genome_coordinate( obj_to_genomes( genomes ) ) ) ) ; }

WRAP( p_sort_by_name, ( Object mem, Object handles ), (mem,handles) )
{ return wrap_stream( new SortingStream<by_seqid>(
			Get_Integer( mem ) * 1024U * 1024U, Get_Integer( handles ) ) ) ; }

WRAP( p_filter_by_length, ( Object l ), (l) ) { return wrap_stream( new LengthFilter( Get_Integer( l ) ) ) ; }

WRAP( p_filter_by_score, ( Object slope, Object len, Object genomes ), (slope,len,genomes) )
{ return wrap_stream( new ScoreFilter( Get_Double( slope ), Get_Double( len ), obj_to_genomes( genomes ) ) ) ; }

WRAP( p_filter_by_mapq, ( Object mapq, Object genomes ), (mapq,genomes) )
{ return wrap_stream( new MapqFilter( obj_to_genomes( genomes ), Get_Integer( mapq ) ) ) ; }

WRAP( p_require_best_hit, ( Object genomes, Object sequences ), (genomes,sequences) )
{ return wrap_stream( new RequireBestHit( obj_to_genomes( genomes ), obj_to_genomes( sequences ) ) ) ; }

WRAP( p_require_hit, ( Object genomes, Object sequences ), (genomes,sequences) )
{ return wrap_stream( new RequireHit( obj_to_genomes( genomes ), obj_to_genomes( sequences ) ) ) ; }

WRAP( p_ignore_hit, ( Object genomes, Object sequences ), (genomes,sequences) )
{ return wrap_stream( new IgnoreHit( obj_to_genomes( genomes ), obj_to_genomes( sequences ) ) ) ; }

WRAP( p_only_genome, ( Object genomes ), (genomes) ) { return wrap_stream( new OnlyGenome( obj_to_genomes( genomes ) ) ) ; }
WRAP( p_filter_multi, ( Object m ), (m) ) { return wrap_stream( new MultiFilter( Get_Integer( m ) ) ) ; }
WRAP( p_subsample, ( Object r ), (r) ) { return wrap_stream( new Subsample( Get_Double( r ) ) ) ; }
WRAP( p_sanitize, (), () ) { return wrap_stream( new Sanitizer() ) ; }
WRAP( p_edit_header, ( Object e ), (e) ) { return wrap_stream( new RepairHeaderStream( object_to_string( e, "" ) ) ) ; }
WRAP( p_inside_region, ( Object f ), (f) ) { return wrap_stream( new InsideRegion( open_any_input_std( f ) ) ) ; }
WRAP( p_outside_region, ( Object f ), (f) ) { return wrap_stream( new OutsideRegion( open_any_input_std( f ) ) ) ; }


// Mergers  (I'll drop the unholy mega-merge here, that's far better
// scripted in Scheme)

WRAP( p_merge, (  int argc, Object *argv ), (argc,argv) ) { return wrap_streams( new MergeStream,   argc, argv ) ; }
WRAP( p_join, (   int argc, Object *argv ), (argc,argv) ) { return wrap_streams( new BestHitStream, argc, argv ) ; }
WRAP( p_concat, ( int argc, Object *argv ), (argc,argv) ) { return wrap_streams( new ConcatStream,  argc, argv ) ; }

// Composition

//! \brief top-level ELK call.
//! Gets one input stream and many output streams, copies between them.
//! Might one day extract a tree of results and return it, but for now
//! returns the final exit code (just like anfo-tool).
WRAP( p_anfo_run, ( int argc, Object *argv ), (argc,argv) )
{
	StreamHolder inp = obj_to_stream( argv[0] ) ;
	FanOut out ;
	for( Object *o = argv+1 ; o != argv+argc ; ++o )
		out.add_stream( obj_to_stream( *o ) ) ;

	out.put_header( inp->fetch_header() ) ;
	while( inp->get_state() == Stream::have_output && out.get_state() == Stream::need_input )
		out.put_result( inp->fetch_result() ) ;
	out.put_footer( inp->fetch_footer() ) ;
	return Make_Integer( out.fetch_footer().exit_code() ) ;
}

//! \brief stream composition.
//! Gets many streams as argument, ties them into a chain and returns a
//! new stream.
WRAP( p_chain, ( int argc, Object *argv ), (argc,argv) )
{
	Holder< Compose > c( new Compose ) ;
	for( Object *o = argv ; o != argv+argc ; ++o )
		c->add_stream( obj_to_stream( *o ) ) ;
	return wrap_stream( c ) ;
}

WRAP( p_read_file, ( Object fn, Object sol_scores, Object ori ), (fn,sol_scores,ori) )
{ return wrap_stream( obj_to_stream( fn, Truep( sol_scores ), Get_Integer( ori ) ) ) ; }

// init code

void elk_finit_libanfo() {}
void elk_init_libanfo() 
{
	t_stream = Define_Type( 0, "anfo-stream", 
        0, sizeof( StreamWrapper ),
		compare_stream_wrappers, 
		compare_stream_wrappers, 
		print_stream_wrapper, 0 ) ;

	Define_Primitive( (P)p_is_stream,        "anfo-stream?",        1, 1, EVAL ) ;
	Define_Primitive( (P)p_set_verbosity,    "set-verbosity",       1, 1, EVAL ) ;
	Define_Primitive( (P)p_use_mmap,         "use-mmap",            1, 1, EVAL ) ;

	Define_Primitive( (P)p_read_file,        "read-file",           3, 3, EVAL ) ;
	Define_Primitive( (P)p_write_native,     "write-native",        2, 2, EVAL ) ;
	Define_Primitive( (P)p_write_text,       "write-text",          1, 1, EVAL ) ;
	Define_Primitive( (P)p_write_sam,        "write-sam",           1, 1, EVAL ) ;
	Define_Primitive( (P)p_write_glz,        "write-glz",           1, 1, EVAL ) ;
	Define_Primitive( (P)p_write_3aln,       "write-three-aln",     1, 1, EVAL ) ;
	Define_Primitive( (P)p_write_fastq,      "write-fastq",         1, 1, EVAL ) ;
	Define_Primitive( (P)p_write_table,      "write-table",         1, 1, EVAL ) ;
	Define_Primitive( (P)p_write_fasta,      "write-fasta",         1, 1, EVAL ) ;

	Define_Primitive( (P)p_duct_tape,        "duct-tape",           1, 1, EVAL ) ;
	Define_Primitive( (P)p_add_alns,         "add-alns",            1, 1, EVAL ) ;
	Define_Primitive( (P)p_rmdup,            "rmdup",               3, 3, EVAL ) ;
	Define_Primitive( (P)p_write_stats,      "write-stats",         1, 1, EVAL ) ;

	Define_Primitive( (P)p_filter_by_length, "filter-length",       1, 1, EVAL ) ;
	Define_Primitive( (P)p_filter_by_score,  "filter-score",        3, 3, EVAL ) ;
	Define_Primitive( (P)p_filter_by_mapq,   "filter-mapq",         2, 2, EVAL ) ;
	Define_Primitive( (P)p_filter_multi,     "filter-multiplicity", 1, 1, EVAL ) ;
	Define_Primitive( (P)p_subsample,        "subsample",           1, 1, EVAL ) ;
	Define_Primitive( (P)p_sanitize,         "sanitize",            0, 0, EVAL ) ;
	Define_Primitive( (P)p_edit_header,      "edit-header",         1, 1, EVAL ) ;
	Define_Primitive( (P)p_inside_region,    "inside-region",       1, 1, EVAL ) ;
	Define_Primitive( (P)p_outside_region,   "outside-region",      1, 1, EVAL ) ;
	Define_Primitive( (P)p_require_best_hit, "require-best-hit",    2, 2, EVAL ) ;
	Define_Primitive( (P)p_require_hit,      "require-hit",         2, 2, EVAL ) ;
	Define_Primitive( (P)p_ignore_hit,       "ignore-hit",          2, 2, EVAL ) ;
	Define_Primitive( (P)p_only_genome,		 "only-genome",         1, 1, EVAL ) ;

	Define_Primitive( (P)p_sort_by_pos,      "sort-pos",            3, 3, EVAL ) ;
	Define_Primitive( (P)p_sort_by_name,     "sort-name",           3, 3, EVAL ) ;

	Define_Primitive( (P)p_merge,            "merge",               1, MANY, VARARGS ) ; 
	Define_Primitive( (P)p_join,             "join",                1, MANY, VARARGS ) ; 
	Define_Primitive( (P)p_concat,           "concat",              1, MANY, VARARGS ) ; 

	Define_Primitive( (P)p_anfo_run,         "anfo-run",            2, MANY, VARARGS ) ;
	Define_Primitive( (P)p_chain,            "chain", 	            1, MANY, VARARGS ) ;
}

} // extern C

