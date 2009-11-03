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

struct SnpRec {
	int hsa_pos ;
	int ptr_pos ;

	string nt_hsa_bases ; // too bad we may observe more than a single base
	string nt_ptr_bases ;

	unsigned char hsa_length ; // one for SNPs, zero for inserts(!), length for deletions
	unsigned char ptr_length ;

	unsigned char hsa_strand : 1 ;
	unsigned char ptr_strand : 1 ;
	unsigned char gap_near_flag : 1 ;
	unsigned char hsa_seen : 1 ; // set as soon as anything was seen
	unsigned char ptr_seen : 1 ;

	Ambicode hsa_chr ;	// index into symtab
	Ambicode ptr_chr ;	// index into symtab

	char hsa_base ;
	char ptr_base ;
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
		ss >> hsa_base >> hsa_chr >> hsa_strand_code >> r.hsa_pos
		   >> ptr_base >> ptr_chr >> ptr_strand_code >> r.ptr_pos
		   >> r.out_base >> flags ;
		if( !s ) break ;

		r.hsa_chr = lookup_sym( hsa_chr ) ;
		r.ptr_chr = lookup_sym( ptr_chr ) ;

		r.hsa_strand = hsa_strand_code == '+' ;
		r.ptr_strand = ptr_strand_code == '+' ;
		r.hsa_length = 1 ;
		r.ptr_length = 1 ;
		r.hsa_base = maybe_compl( r.hsa_strand, to_ambicode( hsa_base ) ) ;
		r.ptr_base = maybe_compl( r.ptr_strand, to_ambicode( ptr_base ) ) ;

		if( !flags.empty() ) error( 0, 0, "Flags field not empty at %s:%d", r.hsa_chr, r.hsa_pos ) ;
		r.gap_near_flag = 0 ;
		r.hsa_seen = 0 ;
		r.ptr_seen = 0 ;
	}
	delete d.back() ;
	d.pop_back() ;
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
		
		s << from_ambicode( maybe_compl( r.hsa_strand, r.hsa_base ) ) << '\t' 
		  << symbols[r.hsa_chr] << '\t' << (r.hsa_strand ? '+' : '-') << '\t' << r.hsa_pos << '\t'
		  << from_ambicode( maybe_compl( r.ptr_strand, r.ptr_base ) ) << '\t'
		  << symbols[r.ptr_chr] << '\t' << (r.ptr_strand ? '+' : '-') << '\t' << r.ptr_pos << '\t'
		  << r.out_base << '\t' << (r.gap_near_flag ? "GC" : "") << '\t' ;

		if( r.hsa_seen ) 
		{
			if( r.nt_hsa_bases.empty() ) s << '-' ;
			else s << r.nt_hsa_bases ;
		}
		s << '\t' ;
		if( r.ptr_seen )
		{
			if( r.nt_ptr_bases.empty() ) s << '-' ;
			else s << r.nt_ptr_bases ;
		}
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
		printf( "Usage: %s hg18.anfo pt2.anfo\n" 
			"  Reads a 'Martin table' from stdin, fills in information from\n"
			"  the two files on the command line assumed to contain McAssemblies\n"
			"  and writes an augmented 'Martin table' to stdout.\n", argv[0] ) ;
		return 1 ;
	}
	std::deque<SnpRec*> mt ;
	if( !read_martin_table( cin, mt ).eof() )
		error( 1, errno, "Parse error in 'Martin table'." ) ;
	// cerr << mallinfo().arena << endl ;

	{
		Chan progress ;

		sort( mt.begin(), mt.end(), ByPt2Coordinate() ) ;
		auto_ptr<Stream> pt2_file( make_input_stream( argv[2] ) ) ;
		deque<SnpRec*>::iterator first_snp = mt.begin() ; // first SNP that hasn't been processed completely
		while( pt2_file->get_state() == Stream::have_output && first_snp != mt.end() )
		{
			Result res = pt2_file->fetch_result() ;
			// check for overlap, extract base
			// note: this is effectively broken for RC'ed alignments (but
			// that doesn't matter in *this* application).

			if( has_hit_to( res, "pt2" ) )
			{
				const Hit &h = hit_to( res, "pt2" ) ;
				unsigned char cur_chr = lookup_sym( h.sequence() ) ;

				// skip ahead to correct chromosome
				while( first_snp != mt.end() && (*first_snp)->ptr_chr != cur_chr )
					++first_snp ;

				if( (first_snp - mt.begin()) % 1024 == 0 )
				{
					stringstream s ;
					s << "SNPs vs. pt2: " << first_snp - mt.begin() << '/' << mt.size() ;
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
					while( first_snp != mt.end() && (*first_snp)->ptr_chr == cur_chr
							&& (*first_snp)->ptr_pos + (*first_snp)->ptr_length + gap_buffer <= ref_pos )
						++first_snp ;

					// bail if we left the current chromosome or no SNPs are left
					if( first_snp == mt.end() || (*first_snp)->ptr_chr != cur_chr ) break ;

					switch( cigar_op( h.cigar(cigar_maj) ) )
					{
						case Hit::Delete: // nothing here, don't need to extract anything; just increment
							// However, sanity check if possible and set GC flag
							for( deque<SnpRec*>::iterator cur_snp = first_snp ;
									cur_snp != mt.end() && (*cur_snp)->ptr_chr == cur_chr &&
									(*cur_snp)->ptr_pos - gap_buffer < ref_pos ; ++cur_snp )
							{
								if( (*cur_snp)->ptr_pos == ref_pos && (*cur_snp)->ptr_base ) 
									if( *ref != (*cur_snp)->ptr_base ) 
										error( 1, 0, "at %d, %c != %c -- wrong coordinate system?", 
												ref_pos, from_ambicode(*ref), from_ambicode( (*cur_snp)->ptr_base ) ) ;

								if( (*cur_snp)->ptr_pos <= ref_pos && ref_pos < (*cur_snp)->ptr_pos + (*cur_snp)->ptr_length )
									(*cur_snp)->ptr_seen =  1 ;
								else 
									(*cur_snp)->gap_near_flag = 1 ;
							}

							++cigar_min ;
							++ref_pos ;
							++ref ;
							break ;

						case Hit::Insert: // inserted sequence: extract to current SNP if pos'n is right
							// Nothing to check sanity against, but set GC flag
							for( deque<SnpRec*>::iterator cur_snp = first_snp ;
									cur_snp != mt.end() && (*cur_snp)->ptr_chr == cur_chr &&
									(*cur_snp)->ptr_pos - gap_buffer < ref_pos ; ++cur_snp )
							{
								if( (*cur_snp)->ptr_pos <= ref_pos && ref_pos < (*cur_snp)->ptr_pos + (*cur_snp)->ptr_length )
								{
									(*cur_snp)->ptr_seen =  1 ;
									(*cur_snp)->nt_ptr_bases.push_back( res.read().sequence()[qry_pos] ) ;
								}
								else 
									(*cur_snp)->gap_near_flag = 1 ;
							}

							++cigar_min ;
							++qry_pos ;
							break ;

						case Hit::Match:
						case Hit::Mismatch: // matched up sequence: check sanity, extract, don't set any flags
							for( deque<SnpRec*>::iterator cur_snp = first_snp ;
									cur_snp != mt.end() && (*cur_snp)->ptr_chr == cur_chr &&
									(*cur_snp)->ptr_pos - gap_buffer <= ref_pos ; ++cur_snp )
							{
								if( (*cur_snp)->ptr_pos == ref_pos && (*cur_snp)->ptr_base ) 
									if( *ref != (*cur_snp)->ptr_base ) 
										error( 1, 0, "at %d, %c != %c -- wrong coordinate system?", 
												ref_pos, from_ambicode(*ref), from_ambicode( (*cur_snp)->ptr_base ) ) ;

								if( (*cur_snp)->ptr_pos <= ref_pos && ref_pos < (*cur_snp)->ptr_pos + (*cur_snp)->ptr_length )
								{
									(*cur_snp)->ptr_seen =  1 ;
									(*cur_snp)->nt_ptr_bases.push_back( res.read().sequence()[qry_pos] ) ;
								}
							}

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

	sort( mt.begin(), mt.end(), ByHg18Coordinate() ) ;
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

