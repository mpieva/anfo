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

// don't do anything unless Elk is present
#if HAVE_ELK_SCHEME_H

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

#include <glob.h>
#include <sys/resource.h>

// Somewhat silly workaround for linking problems: protobuf needs
// pthread, but doesn't automtically link it; we don't need it, and Elk
// won't provide it.  So here's a dependency on libpthread, it gets
// linked and protobuf is happy again.
#include <pthread.h>
pthread_once_t anfo_elk_once_control = PTHREAD_ONCE_INIT;
extern "C" void anfo_elk_init_routine() {}

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
	StreamWrapper* s = (StreamWrapper*)POINTER(o) ;
    ((StreamWrapper*)POINTER(o))->~StreamWrapper() ;
	Deregister_Object( o ) ;
	return Void ;
}

extern "C" Object p_delete_stream( Object o ) 
{
	StreamWrapper* s = (StreamWrapper*)POINTER(o) ;
	if( !s->h_ ) Primitive_Error( "BUG: double deletion of stream" ) ;
	StreamHolder().swap( s->h_ ) ;
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
				if( !nm.empty() && nm[0] == '|' )
					return make_PipeOutputStream( nm.substr(1) ) ;
				else {
					FileOutputStream *s = new FileOutputStream( 
							throw_errno_if_minus1(
								open( nm.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0666 ),
								"opening file" ) ) ;
					s->SetCloseOnDelete( true ) ;
					return make_pair( s, nm ) ;
				}
			}
		case T_Fixnum:
			return make_pair( new FileOutputStream( Get_Exact_Integer( o ) ), "<pipe>" ) ;

		case T_Boolean:
			if( !Truep(o) ) return make_pair( new FileOutputStream( 1 ), "<stdout>" ) ;

		default:
			throw "can't handle file argument" ;
	}
}

pair< std::ostream*, string > open_any_output_std( Object o )
{
	pair< ZeroCopyOutputStream*, string > p = open_any_output_zc( o ) ;
	return make_pair( new zero_copy_ostream( p.first ), p.second ) ;
}

// note: this is only used to initialize the "region filters".  Use of a
// scheme data structure might be beneficial here.
pair< std::istream*, string > open_any_input_std( Object o )
{
	switch( TYPE(o) )
	{
		case T_Symbol:
		case T_String:
			return make_pair( new ifstream( object_to_string(o).c_str() ), object_to_string(o) ) ;

		case T_Boolean:
			if( !Truep(o) ) return make_pair( new istream( cin.rdbuf() ), "<stdin>" ) ;

		default:
			break ;
	}
	throw "can't handle file argument" ;
}

inline StreamHolder obj_to_stream( Object o ) // , bool sol = false , int ori = 33 )
{
	Check_Type( o, t_stream ) ;
	StreamHolder& h = ((StreamWrapper*)POINTER(o))->h_ ;
	if( !h ) Primitive_Error( "BUG: use of deleted stream" ) ;
	return h ;
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

int print_stream_wrapper( Object o, Object p, int, int, int )
{
	StreamWrapper *w = (StreamWrapper*)POINTER(o) ;
	Printf( p, "#[anfo-stream %s %p]",
			w && w->h_ ? w->h_->type_name().c_str() : "Null",
			w ? (const void*)(w->h_) : 0 ) ;
	return 0 ; 
}

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
WRAP( p_write_text,   ( Object f ), (f) ) { return wrap_stream( new TextWriter( open_any_output_std( f ) ) ) ; }
WRAP( p_write_sam,    ( Object f ), (f) ) { return wrap_stream( new SamWriter( open_any_output_std( f ) ) ) ; } 
WRAP( p_write_glz,    ( Object f ), (f) ) { return wrap_stream( new GlzWriter( open_any_output_zc( f ) ) ) ; }
WRAP( p_write_3aln,   ( Object f ), (f) ) { return wrap_stream( new ThreeAlnWriter( open_any_output_std( f ) ) ) ; }
WRAP( p_write_fastq,  ( Object f ), (f) ) { return wrap_stream( new FastqWriter( open_any_output_std( f ), true ) ) ; }
WRAP( p_write_table,  ( Object f ), (f) ) { return wrap_stream( new TableWriter( open_any_output_std( f ) ) ) ; }
WRAP( p_write_fasta,  ( Object f ), (f) ) { return wrap_stream( new FastaAlnWriter( open_any_output_std( f ) ) ) ; }
WRAP( p_write_wig,    ( Object f ), (f) ) { return wrap_stream( new WigCoverageWriter( open_any_output_std( f ) ) ) ; }

// Processors
WRAP( p_duct_tape, ( Object n ), (n) ) { return wrap_stream( new DuctTaper( object_to_string( n, "contig" ) ) ) ; }
WRAP( p_add_alns, ( Object c ), (c) ) { return wrap_stream( new GenTextAlignment( Get_Integer( c ) ) ) ; }
WRAP( p_rmdup, ( Object s, Object i, Object q ), (s,i,q) ) { return wrap_stream( new RmdupStream( Get_Double(s), Get_Double(i), Get_Integer(q) ) ) ; }

// Evaluators
WRAP( p_stats, (), () ) { return wrap_stream( new StatStream() ) ; }
WRAP( p_mismatches, (), () ) { return wrap_stream( new MismatchStats() ) ; }
WRAP( p_divergence, ( Object p, Object s, Object b ), (p,s,b) ) { return wrap_stream( new DivergenceStream( object_to_string(p), object_to_string(s), Get_Integer(b) ) ) ; }

// Filters
WRAP( p_sort_by_pos, ( Object mem, Object handles, Object genomes ), (mem,handles,genomes) )
{ return wrap_stream( new SortingStream<by_genome_coordinate>(
			Get_Integer( mem ), Get_Integer( handles ), by_genome_coordinate( obj_to_genomes( genomes ) ) ) ) ; }

WRAP( p_sort_by_name, ( Object mem, Object handles ), (mem,handles) )
{ return wrap_stream( new SortingStream<by_seqid>( Get_Integer( mem ), Get_Integer( handles ) ) ) ; }

WRAP( p_filter_by_score, ( Object slope, Object len, Object genomes ), (slope,len,genomes) )
{ return wrap_stream( new ScoreFilter( Get_Double( slope ), Get_Double( len ), obj_to_genomes( genomes ) ) ) ; }

WRAP( p_filter_by_qual, ( Object qual ), (qual) )
{ return wrap_stream( new QualFilter( Get_Double( qual ) ) ) ; }

WRAP( p_filter_chain, ( Object left, Object right, Object file ), (left,right,file) )
{ return wrap_stream( new AgreesWithChain( object_to_string(left), object_to_string(right), open_any_input_std( file ) ) ) ; }

WRAP( p_mask_by_qual, ( Object qual ), (qual) )
{ return wrap_stream( new QualMasker( Get_Integer( qual ) ) ) ; }

WRAP( p_filter_by_mapq, ( Object mapq, Object genomes ), (mapq,genomes) )
{ return wrap_stream( new MapqFilter( obj_to_genomes( genomes ), Get_Integer( mapq ) ) ) ; }

WRAP( p_require_best_hit, ( Object genomes, Object sequences ), (genomes,sequences) )
{ return wrap_stream( new RequireBestHit( obj_to_genomes( genomes ), obj_to_genomes( sequences ) ) ) ; }

WRAP( p_require_hit, ( Object genomes, Object sequences ), (genomes,sequences) )
{ return wrap_stream( new RequireHit( obj_to_genomes( genomes ), obj_to_genomes( sequences ) ) ) ; }

WRAP( p_ignore_hit, ( Object genomes, Object sequences ), (genomes,sequences) )
{ return wrap_stream( new IgnoreHit( obj_to_genomes( genomes ), obj_to_genomes( sequences ) ) ) ; }

WRAP( p_only_genome,    ( Object g ), (g) ) { return wrap_stream( new OnlyGenome( obj_to_genomes( g ) ) ) ; }
WRAP( p_filter_by_len,  ( Object l ), (l) ) { return wrap_stream( new LengthFilter( Get_Integer( l ) ) ) ; }
WRAP( p_filter_multi,   ( Object m ), (m) ) { return wrap_stream( new MultiFilter( Get_Integer( m ) ) ) ; }
WRAP( p_subsample,      ( Object r ), (r) ) { return wrap_stream( new Subsample( Get_Double( r ) ) ) ; }
WRAP( p_edit_header,    ( Object e ), (e) ) { return wrap_stream( new RepairHeaderStream( object_to_string( e, "" ) ) ) ; }
WRAP( p_inside_region,  ( Object f ), (f) ) { return wrap_stream( new InsideRegion( open_any_input_std( f ) ) ) ; }
WRAP( p_outside_region, ( Object f ), (f) ) { return wrap_stream( new OutsideRegion( open_any_input_std( f ) ) ) ; }
WRAP( p_sanitize,       (),           ( ) ) { return wrap_stream( new Sanitizer() ) ; }

WRAP( p_get_summary,    ( Object o ), (o) ) { return obj_to_stream( o )->get_summary() ; }


// Mergers  (I'll drop the unholy mega-merge here, that's far better
// scripted in Scheme)

WRAP( p_merge,  ( int argc, Object *argv ), (argc,argv) ) { return wrap_streams( new MergeStream,    argc, argv ) ; }
WRAP( p_join,   ( int argc, Object *argv ), (argc,argv) ) { return wrap_streams( new NearSortedJoin, argc, argv ) ; }
WRAP( p_concat, ( int argc, Object *argv ), (argc,argv) ) { return wrap_streams( new ConcatStream,   argc, argv ) ; }

WRAP( p_read_file, ( Object fn, Object sol_scores, Object origin ), (fn,sol_scores,origin) )
{
	bool sol = Truep( sol_scores ) ;
	int ori = Get_Integer( origin ) ;
	switch( TYPE(fn) )
	{
		// case T_Primitive:
		// case T_Compound:
			// return obj_to_stream( Funcall( o, Null, 0 ), sol, ori ) ;

		case T_Symbol:
		case T_String:
			return wrap_stream( new UniversalReader( object_to_string(fn), 0, sol, ori ) ) ;

		case T_Fixnum:
			return wrap_stream( new UniversalReader( "<pipe>", new FileInputStream( Get_Integer(fn) ), sol, ori ) ) ;

		case T_Boolean:
			if( !Truep(fn) ) return wrap_stream( new UniversalReader( "<stdin>", new FileInputStream(0), sol, ori ) ) ;

		default:
			Primitive_Error( "can't handle file argument ~s", fn ) ;
	}
}

// Composition

//! \brief top-level ELK call.
//! Gets exactly one input stream and one output stream, transfers
//! between them.  Returns results extracted from them.
WRAP( p_anfo_run, ( Object inpo, Object outo ), (inpo,outo) )
{
	StreamHolder inp = obj_to_stream( inpo ) ;
	StreamHolder out = obj_to_stream( outo ) ;

	out->put_header( inp->fetch_header() ) ;
	while( inp->get_state() == Stream::have_output && out->get_state() == Stream::need_input )
		out->put_result( inp->fetch_result() ) ;
	out->put_footer( inp->fetch_footer() ) ;
	return Null ;
}

//! \brief output stream replication
//! Gets many streams as argument, bundles them up so they all receive
//! the same input.
WRAP( p_tee, ( int argc, Object *argv ), (argc,argv) )
{
	Holder< FanOut > f( new FanOut ) ;
	for( Object *o = argv ; o != argv+argc ; ++o )
		f->add_stream( obj_to_stream( *o ) ) ;
	return wrap_stream( f ) ;
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

WRAP( p_glob, ( Object path ), (path) )
{
	glob_t the_glob ;
	throw_errno_if_minus1( glob( Get_String( path ), GLOB_MARK, 0, &the_glob ), "globbing" ) ;

	GC_Node ;
	Object r = Null ;
	GC_Link( r ) ;
	for( int i = the_glob.gl_pathc ; i != 0 ; --i )
		r = Cons( Make_String( the_glob.gl_pathv[i-1], strlen( the_glob.gl_pathv[i-1] ) ), r ) ;
	GC_Unlink ;
	globfree( &the_glob ) ;
	return r ;
}

Object p_version() { return Make_String( PACKAGE_VERSION, strlen(PACKAGE_VERSION) ) ; }
Object p_limit_core( Object amt ) { 
	struct rlimit lim ;
	getrlimit( RLIMIT_AS, &lim ) ;
	lim.rlim_cur = 1024*1024 * Get_Integer( amt ) ;
	setrlimit( RLIMIT_AS, &lim ) ;
	return Void ;
}

// init code

void elk_init_libanfo() 
{
	t_stream = Define_Type( 0, "anfo-stream", 
        0, sizeof( StreamWrapper ),
		compare_stream_wrappers, 
		compare_stream_wrappers, 
		print_stream_wrapper, 0 ) ;

	Define_Primitive( (P)p_is_stream,        "anfo-stream?",        1, 1, EVAL ) ;
	Define_Primitive( (P)p_set_verbosity,    "set-verbosity!",      1, 1, EVAL ) ;
	Define_Primitive( (P)p_use_mmap,         "use-mmap!",           1, 1, EVAL ) ;
	Define_Primitive( (P)p_limit_core,       "limit-core!",         1, 1, EVAL ) ;
	Define_Primitive( (P)p_get_summary,      "get-summary",         1, 1, EVAL ) ;

	Define_Primitive( (P)p_read_file,        "prim-read-file",      3, 3, EVAL ) ;
	Define_Primitive( (P)p_write_native,     "prim-write-native",   2, 2, EVAL ) ;
	Define_Primitive( (P)p_write_text,       "prim-write-text",     1, 1, EVAL ) ; 
	Define_Primitive( (P)p_write_sam,        "prim-write-sam",      1, 1, EVAL ) ;
	Define_Primitive( (P)p_write_glz,        "prim-write-glz",      1, 1, EVAL ) ;
	Define_Primitive( (P)p_write_3aln,       "prim-write-threealn", 1, 1, EVAL ) ;
	Define_Primitive( (P)p_write_fastq,      "prim-write-fastq",    1, 1, EVAL ) ;
	Define_Primitive( (P)p_write_table,      "prim-write-table",    1, 1, EVAL ) ;
	Define_Primitive( (P)p_write_fasta,      "prim-write-fasta",    1, 1, EVAL ) ;
	Define_Primitive( (P)p_write_wig,        "prim-write-wiggle",   1, 1, EVAL ) ;

	Define_Primitive( (P)p_duct_tape,        "prim-duct-tape",      1, 1, EVAL ) ;
	Define_Primitive( (P)p_add_alns,         "prim-add-alns",       1, 1, EVAL ) ;
	Define_Primitive( (P)p_rmdup,            "prim-rmdup",          3, 3, EVAL ) ;

	Define_Primitive( (P)p_filter_by_len,    "prim-filter-length",  1, 1, EVAL ) ;
	Define_Primitive( (P)p_filter_by_qual, 	 "prim-filter-qual",    1, 1, EVAL ) ;
	Define_Primitive( (P)p_mask_by_qual,     "prim-mask-qual",      1, 1, EVAL ) ;
	Define_Primitive( (P)p_filter_multi,     "prim-filter-multi", 	1, 1, EVAL ) ;
	Define_Primitive( (P)p_subsample,        "prim-subsample",      1, 1, EVAL ) ;
	Define_Primitive( (P)p_sanitize,         "prim-sanitize",       0, 0, EVAL ) ;
	Define_Primitive( (P)p_edit_header,      "prim-edit-header",    1, 1, EVAL ) ;
	Define_Primitive( (P)p_filter_by_score,  "prim-filter-score",   3, 3, EVAL ) ;
	Define_Primitive( (P)p_filter_by_mapq,   "prim-filter-mapq",    2, 2, EVAL ) ;
	Define_Primitive( (P)p_filter_chain, 	 "prim-filter-chain",   3, 3, EVAL ) ;  // redesign?
	Define_Primitive( (P)p_inside_region,    "prim-inside-region",  1, 1, EVAL ) ;  // redesign?
	Define_Primitive( (P)p_outside_region,   "prim-outside-region", 1, 1, EVAL ) ;  // redesign?
	Define_Primitive( (P)p_require_best_hit, "prim-require-bht",    2, 2, EVAL ) ;
	Define_Primitive( (P)p_require_hit,      "prim-require-hit",    2, 2, EVAL ) ;
	Define_Primitive( (P)p_ignore_hit,       "prim-ignore-hit",     2, 2, EVAL ) ;
	Define_Primitive( (P)p_only_genome,		 "prim-only-genomes",   1, 1, EVAL ) ;

	Define_Primitive( (P)p_sort_by_pos,      "prim-sort-pos",       3, 3, EVAL ) ;
	Define_Primitive( (P)p_sort_by_name,     "prim-sort-name",      2, 2, EVAL ) ;

	Define_Primitive( (P)p_merge,            "prim-merge",          1, MANY, VARARGS ) ; 
	Define_Primitive( (P)p_join,             "prim-join",           1, MANY, VARARGS ) ; 
	Define_Primitive( (P)p_concat,           "prim-concat",         1, MANY, VARARGS ) ; 

	Define_Primitive( (P)p_stats,            "prim-stats",          0, 0, EVAL ) ;
	Define_Primitive( (P)p_divergence,       "prim-divergence",     3, 3, EVAL ) ;
	Define_Primitive( (P)p_mismatches,       "prim-mismatches",     0, 0, EVAL ) ;

	Define_Primitive( (P)p_delete_stream,    "prim-delete-stream",  1, 1, EVAL ) ;
	Define_Primitive( (P)p_anfo_run,         "prim-anfo-run",       2, 2, EVAL ) ;
	Define_Primitive( (P)p_tee,              "prim-tee",            1, MANY, VARARGS ) ;
	Define_Primitive( (P)p_chain,            "prim-chain", 	        1, MANY, VARARGS ) ;

	// not exactly Anfo, but damn practical
	Define_Primitive( (P)p_glob,             "glob",                1, 1, EVAL ) ;
	Define_Primitive( (P)p_version,          "anfo-version",        0, 0, EVAL ) ;

	// part of the aforementioned workaround
	pthread_once( &anfo_elk_once_control, anfo_elk_init_routine ) ;
}

void elk_finit_libanfo()
{
    // force GC to run d'tors, makes sure nothing was forgotten
    Terminate_Type( t_stream ) ;
}

} // extern C

#endif
