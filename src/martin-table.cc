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

#include "index.h"
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
#include <vector>

#include <error.h>
#include <malloc.h>

static const int gap_buffer = 5 ; // must be this far from a gap to not set a flag

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
		if( i == 256 ) error( 1, 0, "symbol table full" ) ;
	}
	return i-1 ;
}

struct SnpRec1 {
	int pos ;
	string nt_bases ;					// too bad we may observe more than a single base
	unsigned char length ;				// one for SNPs, zero for inserts(!), length for deletions --- XXX mustn't be too big
	unsigned char chr ;					// symbol index
	unsigned char base : 4 ;			// base as an ambiguity code
	unsigned char strand : 1 ;			// is this on the 'sense' strand?
	unsigned char seen : 1 ;			// set once we have an observation
	unsigned char gap_near_flag : 1 ;	// set if a gap was near (±5nt)
	unsigned char edge_near_flag : 1 ;	// set if contig end was near (±5nt)
} ;

struct SnpRec {
	SnpRec1 hsa, ptr ;
	char out_base ;
} ;

inline Ambicode maybe_compl( bool str, Ambicode a ) { return str ? a : complement(a) ; }

template< typename C >
istream& read_martin_table( istream& s, C& d )
{
	Chan ch ;

	char hsa_strand_code, ptr_strand_code, hsa_base, ptr_base ;
	string flags, hsa_chr, ptr_chr ;
	const int inf = numeric_limits<int>::max() ;
	int x = s.get() ;
	if( x == '#' ) s.ignore( inf, '\n' ) ; // drop header line
	else s.putback( x ) ;

	for(;;)
	{
		if( d.size() % 1024 == 0 )
		{
			stringstream s ;
			s << "reading SNPs: " << d.size() ; 
			ch( Console::info, s.str() ) ;
		}
		d.push_back( new SnpRec() ) ;
		SnpRec& r = *d.back() ;
		string line ;
		if( !getline( s, line ) ) break ;
		stringstream ss( line ) ;
		ss >> hsa_base >> hsa_chr >> hsa_strand_code >> r.hsa.pos
		   >> ptr_base >> ptr_chr >> ptr_strand_code >> r.ptr.pos
		   >> r.out_base >> flags ;
		if( !s ) break ;

		r.hsa.chr = lookup_sym( hsa_chr ) ;
		r.ptr.chr = lookup_sym( ptr_chr ) ;

		r.hsa.strand = hsa_strand_code == '+' ;
		r.ptr.strand = ptr_strand_code == '+' ;
		r.hsa.length = 1 ;
		r.ptr.length = 1 ;
		r.hsa.base = maybe_compl( r.hsa.strand, to_ambicode( hsa_base ) ) ;
		r.ptr.base = maybe_compl( r.ptr.strand, to_ambicode( ptr_base ) ) ;

		if( !flags.empty() ) error( 0, 0, "Flags field not empty at %s:%d", symbols[r.hsa.chr].c_str(), r.hsa.pos ) ;
		r.hsa.gap_near_flag = 0 ;
		r.ptr.gap_near_flag = 0 ;
		r.hsa.seen = 0 ;
		r.ptr.seen = 0 ;
	}
	delete d.back() ;
	d.pop_back() ;
	return s ;
}

inline ostream& write_half_record( ostream& s, const SnpRec1& r )
{
	return s 
		<< from_ambicode( maybe_compl( r.strand, r.base ) ) << '\t' 
		<< symbols[r.chr] << '\t' << (r.strand ? '+' : '-') << '\t'
		<< r.pos << '\t' ;
}
inline ostream& write_bases( ostream& s, const SnpRec1& r )
{
	s << '\t' ;
	if( r.seen && r.nt_bases.empty() ) s << '-' ;
	else if( r.seen ) s << r.nt_bases ;
	return s ;
}

template< typename C >
ostream& write_martin_table( ostream& s, const C& d )
{
	Chan ch ;
	s << "#HSA_Base	HSA_Chr	HSA_Strand	HSA_Pos	PAN_Base	PAN_Chr	PAN_Strand	PAN_Pos	OutBase	Flag	NEA_BaseH	NEA_BaseC\n" ;
	for( typename C::const_iterator i = d.begin() ; i != d.end() ; ++i )
	{
		if( (i-d.begin()) % 1024 == 0 )
		{
			stringstream s ;
			s << "writing SNPs: " << i-d.begin() << '/' << d.size() ;
			ch( Console::info, s.str() ) ;
		}

		const SnpRec &r = **i ;
		// if( !r.hsa_seen && !r.ptr_seen ) continue ;
		
		write_half_record( s, r.hsa ) ;
		write_half_record( s, r.ptr ) ;
		bool g = r.hsa.gap_near_flag || r.ptr.gap_near_flag ;
		bool e = r.hsa.edge_near_flag || r.ptr.edge_near_flag ;
		s << r.out_base << '\t' << (g ? e ? "GC:EC" : "GC" : e ? "EC" : "") ;
		write_bases( s, r.hsa ) ;
		write_bases( s, r.ptr ) ;
		s << '\n' ;
	}
	return s ;
}

inline bool operator < ( const SnpRec1& l, const SnpRec1& r )
{
	if( symbols[l.chr] < symbols[r.chr] ) return true ;
	if( symbols[r.chr] < symbols[l.chr] ) return false ;
	return l.pos < r.pos ;
}

struct ByHg18Coordinate {
	inline bool operator()( const SnpRec *l, const SnpRec *r ) {
		return l->hsa < r->hsa ;
	}
} ;

struct ByPt2Coordinate {
	inline bool operator()( const SnpRec *l, const SnpRec *r ) {
		return l->ptr < r->ptr ;
	}
} ;

struct GetHg18Rec {
	SnpRec1& operator()( vector<SnpRec*>::iterator i ) { return (*i)->hsa ; }
} ;
struct GetPt2Rec {
	SnpRec1& operator()( vector<SnpRec*>::iterator i ) { return (*i)->ptr ; }
} ;

template< typename T >
void scan_anfo_file( vector<SnpRec*> &mt, const char* fn, const char* genome, T get )
{
	Chan progress ;
	auto_ptr<Stream> anfo_file( make_input_stream( fn ) ) ;

	vector<SnpRec*>::iterator first_snp = mt.begin() ; // first SNP that hasn't been processed completely
	while( anfo_file->get_state() == Stream::have_output && first_snp != mt.end() )
	{
		Result res = anfo_file->fetch_result() ;
		// check for overlap, extract base
		// note: this is effectively broken for RC'ed alignments (but
		// that doesn't matter in *this* application).

		if( has_hit_to( res, genome ) )
		{
			const Hit &h = hit_to( res, genome ) ;
			unsigned char cur_chr = lookup_sym( h.sequence() ) ;

			// skip ahead to correct chromosome
			while( first_snp != mt.end() && get( first_snp ).chr != cur_chr )
				++first_snp ;

			// if( (first_snp - mt.begin()) % 1024 == 0 )
			{
				stringstream s ;
				s << "SNPs vs. " << genome << ": " << first_snp - mt.begin() << '/' << mt.size() ;
				progress( Console::info, s.str() ) ;
			}

			CompactGenome &g = Metagenome::find_sequence( h.genome_name(), h.sequence(), Metagenome::ephemeral ) ;
			DnaP ref = g.find_pos( h.sequence(), h.start_pos() ) ;
			int cigar_maj = 0, cigar_min = 0, qry_pos = 0, ref_pos = h.start_pos() ;

			while( qry_pos != res.read().sequence().size() )
			{
				while( cigar_maj != h.cigar_size() &&
						cigar_min == cigar_len( h.cigar(cigar_maj) ) )
				{
					cigar_min = 0 ;
					++cigar_maj ;
				}
				if( cigar_maj == h.cigar().size() ) break ; // shouldn't happen (but you never know)

				// skip SNPs that cannot possibly overlap current position
				// (note the +5 -- necessary for the GC flag)
				while( first_snp != mt.end() && get( first_snp ).chr == cur_chr
						&& get( first_snp ).pos + get( first_snp ).length + gap_buffer <= ref_pos )
					++first_snp ;

				// bail if we left the current chromosome or no SNPs are left
				if( first_snp == mt.end() || get( first_snp ).chr != cur_chr ) break ;

				for( vector<SnpRec*>::iterator cur_snp = first_snp ; cur_snp != mt.end() ; ++cur_snp )
				{
					SnpRec1 &snp = get( cur_snp ) ;
					Hit::Operation op = cigar_op( h.cigar(cigar_maj) ) ;
					if( snp.chr != cur_chr || snp.pos - gap_buffer >= ref_pos ) break ;

					// close to contig edge? set flag
					if( qry_pos < gap_buffer || qry_pos >= res.read().sequence().size() - gap_buffer )
						snp.edge_near_flag = 1 ;

					// sanity check: only possible if we're not looking
					// at an insert
					if( snp.pos == ref_pos && snp.base && *ref != snp.base ) switch( op )
					{
						case Hit::Delete:
						case Hit::Match:
						case Hit::Mismatch:
							error( 1, 0, "at %d, %c != %c -- wrong coordinate system?", 
									ref_pos, from_ambicode(*ref), from_ambicode( snp.base ) ) ;
					}
					
					if( snp.pos <= ref_pos && ref_pos < snp.pos + snp.length )
					{
						// SNP observed.  Anything to extract?
						snp.seen = 1 ;
						switch( op )
						{
							case Hit::Insert:
							case Hit::Match:
							case Hit::Mismatch:
								snp.nt_bases.push_back( res.read().sequence()[qry_pos] ) ;
						}
					}
					else switch( op )
					{
						// not observed, but close by.  Set gap flag?
						case Hit::Delete:
						case Hit::Insert:
							snp.gap_near_flag = 1 ;
					}
				}

				switch( cigar_op( h.cigar(cigar_maj) ) )
				{
					case Hit::Delete:
						++cigar_min ;
						++ref_pos ;
						++ref ;
						break ;

					case Hit::Insert:
						++cigar_min ;
						++qry_pos ;
						break ;

					case Hit::Match:
					case Hit::Mismatch:
						++cigar_min ;
						++qry_pos ;
						++ref_pos ;
						++ref ;
						break ;

					default: // anything else: can't do anything  In fact, isn't even supported.
						error( 1, 0, "Unexpected CIGAR operation %d.",  cigar_op( h.cigar(cigar_maj) ) ) ;
				}
			}
		}
	}
}


int main_( int argc, char const **argv )
{
	console.loglevel = Console::debug ;

	if( argc <= 3 ) {
		printf( "Usage: %s [hg18.anfo] [pt2.anfo] [martin-table.tsv...]\n" 
			"  Reads 'Martin tables' from files, fills in information from\n"
			"  the two Anfo files on the command line assumed to contain\n"
			"  McAssemblies and writes an augmented 'Martin table' to stdout.\n", argv[0] ) ;
		return 1 ;
	}
	std::vector<SnpRec*> mt ;
	for( char const **arg = argv+3 ; arg != argv+argc ; ++arg )
	{
		ifstream f( *arg ) ;
		if( !read_martin_table( f, mt ).eof() )
			error( 1, errno, "Parse error in 'Martin table' %s.", *arg ) ;
	}

	// cerr << mallinfo().arena << endl ;

	console.output( Console::info, "Sorting on pt2..." ) ;
	sort( mt.begin(), mt.end(), ByPt2Coordinate() ) ;
	console.output( Console::info, "Done." ) ;
	scan_anfo_file( mt, argv[2], "pt2", GetPt2Rec() ) ;

	console.output( Console::info, "Sorting on hg18..." ) ;
	sort( mt.begin(), mt.end(), ByHg18Coordinate() ) ;
	console.output( Console::info, "Done." ) ;
	scan_anfo_file( mt, argv[1], "hg18", GetHg18Rec() ) ;
	/*
	Stream *hg18_file = make_input_stream( argv[1] ) ;
	first_snp = mt.begin() ;
	while( hg18_file->get_state() == Stream::have_output && first_snp != mt.end() )
	{
		Result res = hg18_file->fetch_result() ;
		// check for overlap, extract base
	}
	delete hg18_file ;
	*/

	write_martin_table( cout, mt ) ;
	return 0 ;
}

