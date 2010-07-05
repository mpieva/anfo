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

#include "compress_stream.h"
#include "ducttape.h"
#include "misc_streams.h"
#include "output_streams.h"
#include "stream.h"
#include "util.h"

#include <ostream>

#include <glob.h>
#include <popt.h>
#include <sys/resource.h>

using namespace std ;
using namespace output ;
using namespace streams ;

//! \page anfo_stream stream-like operations on ANFO files
//!
//! This program builds a chain out of a variety of filters for ANFO
//! streams, allowing various filters, merging and format conversions to
//! be applied.  The options on the command line are read from left to
//! right.  Those that name filters will build a chain of filters to be
//! applied to every input file until a merge operation (merge,
//! mega-merge, or concat) is reached, then they build up a chain to be
//! applied to the result of the merge.  Some filters will write a file
//! containing a subset of the alignments.  Afterwards, all the input
//! files are "sucked" through the filter pipeline, driven by the
//! calculation of a few statistics, which are printed at the end.

struct ParamBlock {
	float slope ;
	float intercept ;
	const char* genome ;
	const char* arg ;

	ParamBlock( float s, float i, const char* g, const char* a )
		: slope(s), intercept(i), genome(g), arg(a) {}
} ;

typedef void (*G)( ostream&, const ParamBlock& ) ;

template< typename S > struct FilterParams : public ParamBlock {
	typedef S* (*F)( const ParamBlock& ) ;

	F maker ;
	G describe ;

	FilterParams( const ParamBlock& p, const char* a, F m, G d )
		: ParamBlock(p.slope, p.intercept, p.genome, a), maker(m), describe(d) {}
} ;

typedef std::vector< FilterParams< Stream > > FilterStack ;

float parse_float( const char* a )
{ 
	char *e ;
	float f = strtod( a, &e ) ;
	if( a && *a && !*e ) return f ;
	throw "expected real number, but found \"" + string(a) + "\"" ;
}
int parse_int( const char* a )
{ 
	char *e ;
	int i = strtol( a, &e, 10 ) ;
	if( a && *a && !*e ) return i ;
	throw "expected integer, but found \"" + string(a) + "\"" ;
}
int parse_int( const char* a, int d )
{ 
	char *e ;
	if( !a || !*a ) return d ;
	int i = strtol( a, &e, 10 ) ;
	if( !*e ) return i ;
	throw "expected integer, but found \"" + string(a) + "\"" ;
}

vector<string> split_string( const string& s )
{
	if( s.empty() ) return vector<string>() ;
	vector<string> t ;

	for( string::size_type r, l = 0 ; ; l = r+1 )
	{
		r = s.find( ':', l ) ;
		t.push_back( s.substr( l, r-l ) ) ;
		if( r == string::npos ) break ;
	}
	return t ;
}
vector<string> split_string( const char* s )
{ return s ? split_string( string(s) ) : vector<string>() ; }

bool is_stdout( const char* a ) { return !a || 0 == strcmp( a, "-" ) ; }
const char* parse_fn( const char* a ) { return is_stdout( a ) ? "<stdout>" : a ; }

pair< ZeroCopyOutputStream*, string > make_output_stream_zc( const char* a )
{
	return is_stdout( a )
		? make_pair( new FileOutputStream( 1 ), "<stdout>" ) 
		: make_pair( new FileOutputStream( throw_errno_if_minus1( 
						creat( a, 0666 ), "opening", a ) ) , a ) ;
}

pair< ostream*, string > make_output_stream_std( const char* a )
{
	return is_stdout( a )
		? make_pair( new ostream( cout.rdbuf() ), "<stdout>" )
		: make_pair( static_cast<ostream*>( new ofstream( a ) ), a ) ;
}

pair< istream*, string > make_input_stream_std( const char* a )
{
	return is_stdout( a )
		? make_pair( new istream( cin.rdbuf() ), "<stdin>" ) 
		: make_pair( static_cast<istream*>( new ifstream( a ) ), a ) ;
}

Stream* mk_sort_by_pos( const ParamBlock& p )
{ return new SortingStream<by_genome_coordinate>( parse_int( p.arg, 1024 ), 256,
		by_genome_coordinate( split_string( p.genome ) ) ) ; }

void desc_sort_by_pos( ostream& ss, const ParamBlock& p )
{ 
	ss << "sort by position on genomes [" ;
	vector<string> gs = split_string( p.genome ) ;
	if( !gs.empty() ) { 
		copy( gs.begin(), gs.end()-1, ostream_iterator<string>( ss, ", " ) ) ;
		ss << gs.back() ;
	}
	ss << "], using " << + parse_int( p.arg, 1024 ) <<  " MB" ;
}

Stream* mk_sort_by_name( const ParamBlock& p )
{ return new SortingStream<by_seqid>( parse_int( p.arg, 1024 ) ) ; }

void desc_sort_by_name( ostream& ss, const ParamBlock& p )
{ ss << "sort by sequence id, using " << parse_int( p.arg, 1024 ) << " MB" ; }

Stream* mk_filter_by_length( const ParamBlock& p )
{ return new LengthFilter( parse_int( p.arg ) ) ; }

void desc_filter_by_length( ostream& ss, const ParamBlock& p )
{ ss << "remove alignments shorter than " << parse_int( p.arg ) ; }

Stream* mk_filter_by_score( const ParamBlock& p )
{ return new ScoreFilter( p.slope, p.intercept, split_string( p.genome ) ) ; }

void desc_filter_by_score( ostream& ss, const ParamBlock& p )
{
	ss << "remove alignments to " << (p.genome?p.genome:"any genome") 
		<< " scoring worse than ( " << p.slope << " * ( L - " << p.intercept << " ) )" ;
}

Stream* mk_filter_by_mapq( const ParamBlock& p )
{ return new MapqFilter( split_string( p.genome ), parse_int( p.arg ) ) ; }

void desc_filter_by_mapq( ostream& ss, const ParamBlock& p )
{
	ss << "remove alignments where MAPQ" ;
	if( p.genome ) ss << " on genome " << p.genome ;
	ss << " is below " << parse_int( p.arg ) ;
}

Stream* mk_filter_by_hit( const ParamBlock& p )
{ return new OnlyGenome( split_string( p.genome ) ) ; }

void desc_filter_by_hit( ostream& ss, const ParamBlock& p )
{
	ss << "remove hits" ;
	if( p.genome && *p.genome ) ss << " to genomes other than " << p.genome ; 
}

Stream* mk_delete_hit( const ParamBlock& p )
{ return new IgnoreHit( split_string( p.genome ), split_string( p.arg ) ) ; }

void desc_delete_hit( ostream& ss, const ParamBlock& p )
{
	ss << "delete hits" ;
	if( p.arg && *p.arg ) ss << " to sequence " << p.arg ;
	if( p.genome && *p.genome ) ss << " in genome " << p.genome ;
}

Stream* mk_require_hit( const ParamBlock& p )
{ return new RequireHit( split_string( p.genome ), split_string( p.arg ) ) ; }

void desc_require_hit( ostream& ss, const ParamBlock& p )
{
	ss << "delete records without hits" ;
	if( p.arg && *p.arg ) ss << " to sequence " << p.arg ;
	if( p.genome && *p.genome ) ss << " in genome " << p.genome ;
}

Stream* mk_filter_qual( const ParamBlock& p )
{ return new QualFilter( parse_float( p.arg ) ) ; }

void desc_filter_qual( ostream& ss, const ParamBlock& p )
{ ss << "remove reads with avg. quality below " << parse_float( p.arg ) ; }

Stream* mk_mask_qual( const ParamBlock& p )
{ return new QualMasker( parse_int( p.arg ) ) ; }

void desc_mask_qual( ostream& ss, const ParamBlock& p )
{ ss << "mask bases with quality below " << parse_int( p.arg ) ; }

Stream* mk_filter_multi( const ParamBlock& p )
{ return new MultiFilter( parse_int( p.arg, 2 ) ) ; }

void desc_filter_multi( ostream& ss, const ParamBlock& p )
{ ss << "retain only sequences that were seen at least " << parse_int( p.arg, 2 ) << " times" ; }

Stream* mk_subsample( const ParamBlock& p )
{ return new Subsample( parse_float( p.arg ) ) ; }

void desc_subsample( ostream& ss, const ParamBlock& p )
{ ss << "subsample a " << parse_float(p.arg) << " fraction of sequences" ; }

Stream* mk_edit_header( const ParamBlock& p )
{ return new RepairHeaderStream( p.arg ? p.arg : "" ) ; }

void desc_edit_header( ostream& ss, const ParamBlock& p )
{ ss << "invoke " << (p.arg?p.arg:" text editor ") << " on stream's header" ; }

Stream* mk_rmdup( const ParamBlock& p )
{ return new RmdupStream( p.slope, p.intercept, parse_int( p.arg, 127 ) ) ; }

void desc_rmdup( ostream& ss, const ParamBlock& p )
{
	ss << "coalesce duplicates as long as score is no worse than ( "
		<< p.slope << " * ( L - " << p.intercept << " ) ), "
		<< "limit Q score to " << parse_int( p.arg, 127 ) ;
}

Stream* mk_add_alns( const ParamBlock& p )
{ return new GenTextAlignment( parse_int( p.arg, 0 ), false ) ; }

void desc_add_alns( ostream& ss, const ParamBlock& p )
{ ss << "add textual alignments with " << parse_int( p.arg, 0 ) << "nt of context." ; }

Stream* mk_regions_only( const ParamBlock& p )
{ return new InsideRegion( make_input_stream_std( p.arg ) ) ; }

void desc_regions_only( ostream& ss, const ParamBlock& p )
{ ss << "keep results only in regions read from " << p.arg ; }

Stream* mk_not_regions( const ParamBlock& p )
{ return new OutsideRegion( make_input_stream_std( p.arg ) ) ; }

void desc_not_regions( ostream& ss, const ParamBlock& p )
{ ss << "keep results outside regions read from " << p.arg ; }

StreamBundle* mk_merge( const ParamBlock& )
{ return new MergeStream() ; }

void desc_merge( ostream& ss, const ParamBlock& )
{ ss << "merge sorted streams" ; }

StreamBundle* mk_join( const ParamBlock& )
{ return new NearSortedJoin() ; }

void desc_join( ostream& ss, const ParamBlock& )
{ ss << "join near-sorted streams and retain best hits to each genome" ; }

StreamBundle* mk_concat( const ParamBlock& )
{ return new ConcatStream() ; }

void desc_concat( ostream& ss, const ParamBlock& )
{ ss << "concatenate streams" ; }

Stream* mk_output( const ParamBlock& p )
{ return is_stdout( p.arg ) ? new ChunkedWriter( 1, 99, "<stdout>" ) : new ChunkedWriter( p.arg, 99 ) ; } 

void desc_output( ostream& ss, const ParamBlock& p )
{ ss << "write native output to " << parse_fn( p.arg ) ; }

Stream* mk_output_text( const ParamBlock& p )
{ return new TextWriter( make_output_stream_std( p.arg ) ) ; }

void desc_output_text( ostream& ss, const ParamBlock& p )
{ ss << "write in text format to " << parse_fn( p.arg ) ; }

Stream* mk_output_sam( const ParamBlock& p )
{ return new SamWriter( make_output_stream_std( p.arg ) ) ; }

void desc_output_sam( ostream& ss, const ParamBlock& p )
{ ss << "write alignments in SAM format to " << parse_fn( p.arg ) ; }

Stream* mk_output_glz( const ParamBlock& p )
{ return new GlzWriter( make_output_stream_zc( p.arg ) ) ; }

void desc_output_glz( ostream& ss, const ParamBlock& p )
{ ss << "write contigs in GLZ format to " << parse_fn( p.arg ) ; }

void desc_output_3aln( ostream& ss, const ParamBlock& p )
{ ss << "write contigs in 3ALN format to " << parse_fn( p.arg ) ; }

Stream* mk_output_fasta( const ParamBlock& p )
{ return new FastaAlnWriter( make_output_stream_std( p.arg ) ) ; }

void desc_output_fasta( ostream& ss, const ParamBlock& p )
{ ss << "write best alignments(!) in FASTA format to " << parse_fn( p.arg ) ; }

Stream* mk_output_fastq( const ParamBlock& p )
{ return new FastqWriter( make_output_stream_std( p.arg ), true ) ; }

void desc_output_fastq( ostream& ss, const ParamBlock& p )
{ ss << "write sequences(!) in FASTQ format to " << parse_fn( p.arg ) ; }

Stream* mk_output_table( const ParamBlock& p )
{ return new TableWriter( make_output_stream_std( p.arg ) ) ; }

void desc_output_table( ostream& ss, const ParamBlock& p )
{ ss << "write useless table to " << parse_fn( p.arg ) ; }

Stream* mk_output_wiggle( const ParamBlock& p )
{ return new WigCoverageWriter( make_output_stream_std( p.arg ) ) ; }

void desc_output_wiggle( ostream& ss, const ParamBlock& p )
{ ss << "write depth-of-coverage in WIG format to " << parse_fn( p.arg ) ; }

Stream* mk_duct_tape( const ParamBlock& p )
{ return new DuctTaper( p.arg ? p.arg : "contig" ) ; }

void desc_duct_tape( ostream& ss, const ParamBlock& p )
{ 
	ss << "mock-assemble hits" ;
	if( p.genome ) ss << " to genome " << p.genome ;
	ss << ", name contigs '" << (p.arg ? p.arg : "contig" ) << '\'' ;
}

Stream* mk_sanitize( const ParamBlock& )
{ return new Sanitizer() ; }

void desc_sanitize( ostream& ss, const ParamBlock& )
{ ss << "remove debugging information" ; }

const char *poptGetOptArg1( poptContext con )
{
	const char *p = poptGetOptArg( con ) ;
	if( !p || *p != '-' || !p[1] ) return p ;

	console.output( Console::warning, string("poptGetOptArg: ") + p + " treated as parameter" ) ;
	return p ;
}

struct optdef {
	FilterParams<Stream>::F filter_maker ;
	FilterParams<StreamBundle>::F merge_maker ;
	FilterParams<Stream>::F output_maker ;
	G description ;
} ;

WRAPPED_MAIN
{
	GOOGLE_PROTOBUF_VERIFY_VERSION ;
	enum {
		opt_none, opt_sort_pos, opt_sort_name, opt_filter_length,
		opt_filter_score, opt_filter_mapq, opt_filter_hit,
		opt_delete_hit, opt_require_hit, opt_filter_qual, opt_mask_qual,
		opt_subsample, opt_sanitize, opt_filter_multi, opt_regions_only,
		opt_not_regions, opt_edit_header, opt_add_alns, opt_rmdup,
		opt_duct_tape, opt_merge, opt_join, opt_concat,
		opt_output, opt_output_text, opt_output_sam, opt_output_glz,
		opt_output_3aln, opt_output_fasta, opt_output_fastq,
		opt_output_table, opt_output_wiggle, opt_version,
	} ;

	struct optdef optdefs[] = {
		{ 0,                   0, 0,                         0 },
		{ mk_sort_by_pos,      0, 0,          desc_sort_by_pos },
		{ mk_sort_by_name,     0, 0,         desc_sort_by_name },
		{ mk_filter_by_length, 0, 0,     desc_filter_by_length },
		{ mk_filter_by_score,  0, 0,      desc_filter_by_score },
		{ mk_filter_by_mapq,   0, 0,       desc_filter_by_mapq },
		{ mk_filter_by_hit,    0, 0,        desc_filter_by_hit },
		{ mk_delete_hit,       0, 0,           desc_delete_hit },
		{ mk_require_hit,      0, 0,          desc_require_hit },
		{ mk_filter_qual,      0, 0,          desc_filter_qual },
		{ mk_mask_qual,        0, 0,            desc_mask_qual },
		{ mk_subsample,        0, 0,            desc_subsample },
		{ mk_sanitize,         0, 0,             desc_sanitize },
		{ mk_filter_multi,     0, 0,         desc_filter_multi },
		{ mk_regions_only,     0, 0,         desc_regions_only },
		{ mk_not_regions,      0, 0,          desc_not_regions },
		{ mk_edit_header,      0, 0,          desc_edit_header },
		{ mk_add_alns,         0, 0,             desc_add_alns },
		{ mk_rmdup,            0, 0,                desc_rmdup },
		{ mk_duct_tape,        0, 0,            desc_duct_tape },
		{ 0, mk_merge,            0,                desc_merge },
		{ 0, mk_join,             0,                 desc_join },
		{ 0, mk_concat,           0,               desc_concat },
		{ 0, 0, mk_output,                         desc_output },
		{ 0, 0, mk_output_text,               desc_output_text },
		{ 0, 0, mk_output_sam,                 desc_output_sam },
		{ 0, 0, mk_output_glz,                 desc_output_glz },
		{ 0, 0, mk_output_fasta,             desc_output_fasta },
		{ 0, 0, mk_output_fastq,             desc_output_fastq },
		{ 0, 0, mk_output_table,             desc_output_table },
		{ 0, 0, mk_output_wiggle,           desc_output_wiggle },
		{ 0, 0, 0, 0 }
	} ;

	ParamBlock param( 7.5, 20.0, 0, 0 ) ;
	int POPT_ARG_DFLT = POPT_ARG_FLOAT | POPT_ARGFLAG_SHOW_DEFAULT ;
	int core_limit = 0, dry_run = 0 ;

	struct poptOption options[] = {
		{ "sort-pos",      's', POPT_ARG_INT,    0, opt_sort_pos,      "sort by alignment position [using <n MiB memory]", "n" },
		{ "sort-name",     'S', POPT_ARG_INT,    0, opt_sort_name,     "sort by read name [using <n MiB memory]", "n" },
		{ "filter-length", 'l', POPT_ARG_INT,    0, opt_filter_length, "filter for length of at least L", "L" },
		{ "filter-score",  'f', POPT_ARG_NONE,   0, opt_filter_score,  "filter for max score", 0 },
		{ "filter-mapq",    0 , POPT_ARG_INT,    0, opt_filter_mapq,   "remove alignments with MAPQ below Q", "Q" },
		{ "only-genome",    0 , POPT_ARG_NONE,   0, opt_filter_hit,    "keep only hits (in G/anywhere) to SEQ/anything", "SEQ" },
		{ "delete-hit",     0 , POPT_ARG_STRING, 0, opt_delete_hit,    "delete hits (in G/anywhere) to SEQ/anything", "SEQ" },
		{ "require-hit",    0 , POPT_ARG_STRING, 0, opt_require_hit,   "delete records without a hit (in G/anywhere) to SEQ/anything", "SEQ" },
		{ "filter-qual",    0 , POPT_ARG_FLOAT,  0, opt_filter_qual,   "delete reads with avg. quality below Q", "Q" },
		{ "mask-qual",      0 , POPT_ARG_INT,    0, opt_mask_qual,     "delete bases with quality below Q", "Q" },
		{ "subsample",      0,  POPT_ARG_FLOAT,  0, opt_subsample,     "subsample a fraction F of the results", "F" },
		{ "sanitize",       0 , POPT_ARG_NONE,   0, opt_sanitize,      "remove debugging information", 0 },
		{ "multiplicity",   0 , POPT_ARG_INT,    0, opt_filter_multi,  "keep reads with multiplicity above N", "N" },
		{ "edit-header",    0 , POPT_ARG_STRING, 0, opt_edit_header,   "invoke editor ED on the stream's header", "ED" },
		{ "inside-regions", 0 , POPT_ARG_STRING, 0, opt_regions_only,  "keep results inside annotated regions only", "FILE" },
		{ "outside-regions",0 , POPT_ARG_STRING, 0, opt_not_regions,   "keep results outside annotated regions only", "FILE" },
		{ "concat",        'c', POPT_ARG_NONE,   0, opt_concat,        "concatenate streams", 0 },
		{ "merge",         'm', POPT_ARG_NONE,   0, opt_merge,         "merge sorted streams", 0 },
		{ "join",          'j', POPT_ARG_NONE,   0, opt_join,          "join streams and retain best hits", 0 },
		{ "rmdup",         'd', POPT_ARG_INT,    0, opt_rmdup,         "remove PCR duplicates, clamp Q-scores to Q", "Q" },
		{ "add-alns",       0 , POPT_ARG_INT,    0, opt_add_alns,      "add alignments with N nt of context", "N" },
		{ "output",        'o', POPT_ARG_STRING, 0, opt_output,        "write native stream to file FILE", "FILE" },
		{ "output-text",    0 , POPT_ARG_STRING, 0, opt_output_text,   "write protobuf text stream to FILE", "FILE" },
		{ "output-sam",     0 , POPT_ARG_STRING, 0, opt_output_sam,    "write alignments in sam format to FILE", "FILE" },
		{ "output-glz",     0 , POPT_ARG_STRING, 0, opt_output_glz,    "write contigs in GLZ format to FILE", "FILE" },
		{ "output-fasta",   0 , POPT_ARG_STRING, 0, opt_output_fasta,  "write alignments(!) in fasta format to FILE", "FILE" },
		{ "output-fastq",   0 , POPT_ARG_STRING, 0, opt_output_fastq,  "write sequences(!) in fastq format to FILE", "FILE" },
		{ "output-table",   0 , POPT_ARG_STRING, 0, opt_output_table,  "write per-alignment stats to FILE", "FILE" },
		{ "output-wiggle",  0 , POPT_ARG_STRING, 0, opt_output_wiggle, "write depth-of-coverage in WIG format to FILE", "FILE" }, 
		{ "duct-tape",      0 , POPT_ARG_STRING, 0, opt_duct_tape,     "mock-assemble into contigs named NAME", "NAME" },

		{ "set-slope",      0 , POPT_ARG_DFLT,   &param.slope,      0, "set slope parameter to S", "S" },
		{ "set-intercept",  0 , POPT_ARG_DFLT,   &param.intercept,  0, "set length discount parameter to L", "L" },
		{ "set-genome",     0 , POPT_ARG_STRING, &param.genome,     0, "set interesting genome parameter to G", "G" },
		{ "clear-genome",   0 , POPT_ARG_VAL,    &param.genome,     0, "clear interesting genome parameter", 0 },

		{ "vmem",           0 , POPT_ARG_INT,    &core_limit,       0, "limit virtual memory to X megabytes", "X" },
		{ "quiet",         'q', POPT_ARG_VAL,    &console.loglevel, Console::error, "suppress most output", 0 },
		{ "verbose",       'v', POPT_ARG_VAL,    &console.loglevel, Console::info,  "produce more output", 0 },
		{ "debug",          0 , POPT_ARG_VAL,    &console.loglevel, Console::debug, "produce debugging output", 0 },
		{ "dry-run",       'n', POPT_ARG_NONE,   &dry_run,          1, "parse command line, then exit", 0 },
		{ "version",       'V', POPT_ARG_NONE,   0, opt_version,       "print version number and exit", 0 },
		POPT_AUTOHELP POPT_TABLEEND
	} ;

	FilterStack filters_initial ;
	FilterStack *filters_current = &filters_initial ;
	FilterParams< StreamBundle > merging_filter( ParamBlock(0,0,0,0), 0, mk_concat, desc_concat ) ;

	typedef std::deque< FilterStack > FilterStacks ;
	FilterStacks filters_terminal ;
	filters_terminal.push_back( FilterStack() ) ;

	poptContext pc = poptGetContext( "anfo", argc, argv, options, 0 ) ;
	poptSetOtherOptionHelp( pc, "[OPTION...] [sequence-file...]" ) ;
	poptReadDefaultConfig( pc, 0 ) ;

	if( argc <= 1 ) { poptPrintHelp( pc, stderr, 0 ) ; return 1 ; }
	for( int rc = poptGetNextOpt( pc ) ; rc != -1 ; rc = poptGetNextOpt(pc) )
	{
		if( rc == opt_version ) {
			std::cout << poptGetInvocationName(pc) << ", revision " << PACKAGE_VERSION << std::endl ;
			return 0 ;
		}
		else if( rc >= 0 && optdefs[rc].filter_maker )
		{
			filters_current->push_back( FilterParams< Stream >(
						param, poptGetOptArg1( pc ), optdefs[rc].filter_maker,
						optdefs[rc].description ) );
		}
		else if( rc >= 0 && optdefs[rc].merge_maker )
		{
			// make sure we are still creating input filters
			if( filters_current != &filters_initial )
				throw "merge-like commands cannot not follow merge- or output-like commands" ;

			merging_filter = FilterParams< StreamBundle >(
				param, poptGetOptArg1( pc ), optdefs[rc].merge_maker, optdefs[rc].description ) ;

			// from now on we build output filters
			filters_current = &filters_terminal.back() ;
		}
		else if( rc >= 0 && optdefs[rc].output_maker )
		{
			// make sure we are creating output filters
			if( filters_current == &filters_initial )
				filters_current = &filters_terminal.back() ;

			// create filter
			filters_current->push_back( FilterParams< Stream >(
				param, poptGetOptArg( pc ), optdefs[rc].output_maker, optdefs[rc].description ) ) ;

			// start new output stream
			filters_terminal.push_back( FilterStack() ) ;
			filters_current = &filters_terminal.back() ;
		}
		else
		{
			std::clog << poptGetInvocationName(pc) << ": " << poptStrerror( rc ) 
				<< ' ' << poptBadOption( pc, 0 ) << std::endl ;
			return 1 ; 
		}
	}

	if( core_limit ) {
		struct rlimit lim ;
		getrlimit( RLIMIT_AS, &lim ) ;
		lim.rlim_cur = 1024*1024 * (long)core_limit ;
		setrlimit( RLIMIT_AS, &lim ) ;
	}

	// iterate over non-option arguments, glob everything
	glob_t the_glob ;
	the_glob.gl_pathv = 0 ;
	the_glob.gl_pathc = 0 ;
	int glob_flag = GLOB_NOSORT | GLOB_NOCHECK ;

	while( const char* arg = poptGetArg( pc ) )
	{
		switch( glob( arg, glob_flag, 0, &the_glob ) ) 
		{
			case 0: break ;
			case GLOB_NOSPACE: throw "out of memory while globbing " + std::string(arg) ;
			case GLOB_ABORTED: throw "read error on " + std::string(arg) ;
			default: throw "strange error from glob()" ;
		}
		glob_flag |= GLOB_APPEND ;
	}

	// give a report of the filter stack (because it's so easy to get
	// wrong)
	console.output( Console::notice, "Input files:" ) ;
	if( the_glob.gl_pathc )
		for( char **arg = the_glob.gl_pathv ; arg != the_glob.gl_pathv + the_glob.gl_pathc ; ++arg )
		{
			if( !strcmp( *arg, "-" ) ) console.output( Console::notice, "  <stdin>" ) ;
			else console.output( Console::notice, "  " + string(*arg) ) ;
		}
	else console.output( Console::notice, "  <stdin>" ) ;

	if( !filters_initial.empty() ) 
	{
		console.output( Console::notice, "For every input file:" ) ;
		for( FilterStack::const_iterator i = filters_initial.begin() ; i != filters_initial.end() ; ++i )
		{
			stringstream s ;
			(i->describe)( s << "  ", *i ) ;
			console.output( Console::notice, s.str() ) ;
		}
	}
	{
		stringstream s ;
		merging_filter.describe( s, merging_filter ) ;
		console.output( Console::notice, s.str() ) ;
	}

	// last filter stack is not empty (== missing output filter) or only
	// one filter stack (== no output at all)?  --> add a writer for
	// stdout to last filter.  else remove the empty one
	if( !filters_terminal.back().empty() || filters_terminal.size() == 1 )
		filters_terminal.back().push_back( FilterParams< Stream >(
			param, 0, mk_output_text, desc_output_text ) ) ;
	else
		filters_terminal.pop_back() ;

	for( FilterStacks::const_iterator i = filters_terminal.begin() ; i != filters_terminal.end() ; ++i )
	{
		console.output( Console::notice, "Filter a copy of the result as follows:" ) ;
		for( FilterStack::const_iterator j = i->begin() ; j != i->end() ; ++j )
		{
			stringstream s ;
			(j->describe)( s << "  ", *j ) ;
			console.output( Console::notice, s.str() ) ;
		}
	}

	poptFreeContext( pc ) ;
	if( !dry_run ) 
	{
		Holder< StreamBundle > merging_stream( (merging_filter.maker)( merging_filter ) ) ;

		// iterate over glob results
		if( the_glob.gl_pathc )
		{
			vector< Holder< Compose > > cs ;
			for( char **arg = the_glob.gl_pathv ; arg != the_glob.gl_pathv + the_glob.gl_pathc ; ++arg )
			{
				Holder< Compose > c( new Compose ) ;
				c->add_stream( 
						strcmp( *arg, "-" ) ? new UniversalReader( *arg )
						: new UniversalReader( "<stdin>", new FileInputStream(0) ) ) ;
				for( FilterStack::const_iterator i = filters_initial.begin() ; i != filters_initial.end() ; ++i )
					c->add_stream( (i->maker)( *i ) ) ;
				cs.push_back( c ) ;
			}
			for( vector< Holder< Compose > >::const_iterator i = cs.begin(), ie = cs.end() ; i != ie ; ++i )
				merging_stream->add_stream( *i ) ;
		}
		else merging_stream->add_stream( new UniversalReader( "<stdin>", new FileInputStream(0) ) ) ;

		FanOut out ;
		for( FilterStacks::const_iterator i = filters_terminal.begin() ; i != filters_terminal.end() ; ++i )
		{
			Holder< Compose > c( new Compose ) ;
			for( FilterStack::const_iterator j = i->begin() ; j != i->end() ; ++j )
				c->add_stream( (j->maker)( *j ) ) ;
			out.add_stream( c ) ;
		}

		return transfer( *merging_stream, out ) ;
	}
	return 0 ;
}

