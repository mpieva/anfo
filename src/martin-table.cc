/* 
 * Finding SNP states for Martin:
 *
 * Input is Martin's list and two files with McAssemblies using hg18 and
 * pt2 as reference, respectively.  Anfo input is sorted, but not
 * indexed, Martin input is neither.  Therefore, scan Martin's stuff,
 * sort it by pt2 coordinate, scan first input file, store SNP states if
 * found.  Resort Martin's stuff by hg18 coordinate, scan second file,
 * store SNP states.  Print resulting table in same format as input,
 * it's automagically sorted on hg18, happiness ensues.
 */

#include "config.h"

#include "stream.h"
#include "output.pb.h"

#include <deque>
#include <istream>
#include <iostream>
#include <fstream>
#include <limits>
#include <map>
#include <ostream>
#include <sstream>
#include <string>

#include <malloc.h>

using namespace std ;
using namespace streams ;
using namespace output ;

static map< string, int > symtab ;
static deque< string > symbols ;

char lookup_sym( const string &s )
{
	int &i = symtab[s] ;
	if( !i )
	{
		symbols.push_back( s ) ;
		i = symbols.size() ;
		if( i == 256 ) {
			cerr << "symbol table full" << endl ;
			exit( 1 ) ;
		}
	}
	return i-1 ;
}

struct SnpRec {
	int hsa_pos ;
	int ptr_pos ;
	unsigned char hsa_strand : 1 ;
	unsigned char ptr_strand : 1 ;
	unsigned char gap_near_flag : 1 ;

	unsigned char hsa_chr ;
	unsigned char ptr_chr ;

	char hsa_base ;
	char ptr_base ;
	char out_base ;

	char nt_hsa_base ;
	char nt_ptr_base ;
} ;

template< typename C >
istream& read_martin_table( istream& s, C& d )
{
	char hsa_strand_code, ptr_strand_code ;
	string flags, hsa_chr, ptr_chr ;
	const int inf = numeric_limits<int>::max() ;
	int x = s.get() ;
	if( x == '#' ) s.ignore( inf, '\n' ) ; // drop header line
	else s.putback( x ) ;

	for(;;)
	{
		d.push_back( new SnpRec() ) ;
		SnpRec& r = *d.back() ;
		string line ;
		if( !getline( s, line ) ) break ;
		stringstream ss( line ) ;
		ss >> r.hsa_base >> hsa_chr >> hsa_strand_code >> r.hsa_pos
		   >> r.ptr_base >> ptr_chr >> ptr_strand_code >> r.ptr_pos
		   >> r.out_base >> flags ;
		if( !s ) break ;

		r.hsa_chr = lookup_sym( hsa_chr ) ;
		r.ptr_chr = lookup_sym( ptr_chr ) ;

		r.hsa_strand = hsa_strand_code == '+' ;
		r.ptr_strand = ptr_strand_code == '+' ;
		if( !flags.empty() ) cerr << "Flags field not empty at " << r.hsa_chr << ':' << r.hsa_pos << endl ;
		r.gap_near_flag = 0 ;
		r.nt_hsa_base = 0 ;
		r.nt_ptr_base = 0 ;
	}
	delete d.back() ;
	d.pop_back() ;
	return s ;
}

template< typename C >
ostream& write_martin_table( ostream& s, const C& d )
{
	s << "#HSA_Base	HSA_Chr	HSA_Strand	HSA_Pos	PAN_Base	PAN_Chr	PAN_Strand	PAN_Pos	OutBase	Flag	NEA_BaseH	NEA_BaseC\n" ;
	for( typename C::const_iterator i = d.begin(), i_end = d.end() ; i != i_end ; ++i )
	{
		const SnpRec &r = **i ;
		s << r.hsa_base << '\t' << symbols[r.hsa_chr] << '\t' << (r.hsa_strand ? '+' : '-') << '\t' << r.hsa_pos << '\t'
	      << r.ptr_base << '\t' << symbols[r.ptr_chr] << '\t' << (r.ptr_strand ? '+' : '-') << '\t' << r.ptr_pos << '\t'
		  << r.out_base << '\t' << (r.gap_near_flag ? "GC" : "") << '\t' ;
		if( r.nt_hsa_base ) s << r.nt_hsa_base ;
		s << '\t' ;
		if( r.nt_ptr_base ) s << r.nt_ptr_base ;
		s << '\n' ;
	}
	return s ;
}

struct ByHg18Coordinate {
	bool operator()( const SnpRec *l, const SnpRec *r ) {
		if( symbols[l->hsa_chr] < symbols[r->hsa_chr] ) return true ;
		if( symbols[r->hsa_chr] < symbols[l->hsa_chr] ) return false ;
		if( l->hsa_pos < r->hsa_pos ) return true ;
		if( r->hsa_pos < l->hsa_pos ) return false ;
		return false ;
	}
} ;

struct ByPt2Coordinate {
	bool operator()( const SnpRec *l, const SnpRec *r ) {
		if( l->ptr_chr < r->ptr_chr ) return true ;
		if( r->ptr_chr < l->ptr_chr ) return false ;
		if( l->ptr_pos < r->ptr_pos ) return true ;
		if( r->ptr_pos < l->ptr_pos ) return false ;
		return false ;
	}
} ;

int main_( int argc, char const **argv )
{
	console.loglevel = Console::debug ;

	if( argc != 3 ) {
		cerr << "Usage: " << argv[0] << " hg18.anfo pt2.anfo\n" 
			"  Reads a 'Martin table' from stdin, fills in information from\n"
			"  the two files on the command line assumed to contain McAssemblies\n"
			"  and writes an augmented 'Martin table' to stdout.\n\n" ;
		return 1 ;
	}
	std::deque<SnpRec*> mt ;
	if( !read_martin_table( cin, mt ).eof() ) {
		cerr << "Parse error in 'Martin table'.\n" ;
		return 1 ;
	}
	// cerr << mallinfo().arena << endl ;

	sort( mt.begin(), mt.end(), ByPt2Coordinate() ) ;
	Stream *pt2_file = make_input_stream( argv[2] ) ;
	deque<SnpRec*>::iterator next_snp = mt.begin() ;
	while( pt2_file->get_state() == Stream::have_output && next_snp != mt.end() )
	{
		Result res = pt2_file->fetch_result() ;
		// check for overlap, extract base
	}
	delete pt2_file ;

	sort( mt.begin(), mt.end(), ByHg18Coordinate() ) ;
	Stream *hg18_file = make_input_stream( argv[1] ) ;
	next_snp = mt.begin() ;
	while( hg18_file->get_state() == Stream::have_output && next_snp != mt.end() )
	{
		Result res = hg18_file->fetch_result() ;
		// check for overlap, extract base
	}
	delete hg18_file ;

	write_martin_table( cout, mt ) ;
	return 0 ;
}

