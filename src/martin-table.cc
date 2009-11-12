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
#include <popt.h>
#include <sys/resource.h>

static int gap_buffer = 5 ; // must be this far from a gap to not set a flag

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
	string nt_bases, nt_qual ; 			// too bad we may observe more than a single base
	unsigned short length ;				// one for SNPs, zero for inserts(!), length for deletions
	unsigned char chr ;					// symbol index
	unsigned char base : 4 ;			// base as an ambiguity code
	unsigned char strand : 1 ;			// is this on the 'sense' strand?
	unsigned char seen : 1 ;			// set once we have an observation
	unsigned char gap_near_flag : 1 ;	// set if a gap was near (±5nt)
	unsigned char edge_near_flag : 1 ;	// set if contig end was near (±5nt)
} ;

struct SnpRec {
	SnpRec1 hsa, ptr ;
	string sequence ;		// for indels: the original sequence
	char out_base ;
} ;

inline Ambicode maybe_compl( bool str, Ambicode a ) { return str ? a : complement(a) ; }

template< typename C > istream& read_martin_table_snp( istream& s, C& d, const char* fn )
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
			s << "reading " << fn << ": " << d.size() ; 
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

template< typename C > istream& read_martin_table_indel( istream& s, C& d, const char* fn )
{
	Chan ch ;

	char hsa_strand_code, ptr_strand_code ;
	string flags, hsa_chr, ptr_chr, type ;
	int hsa_end, ptr_end ;
	const int inf = numeric_limits<int>::max() ;
	int x = s.get() ;
	if( x == '#' ) s.ignore( inf, '\n' ) ; // drop header line
	else s.putback( x ) ;

	for(;;)
	{
		if( d.size() % 1024 == 0 )
		{
			stringstream s ;
			s << "reading " << fn << ": " << d.size() ; 
			ch( Console::info, s.str() ) ;
		}
		d.push_back( new SnpRec() ) ;
		SnpRec& r = *d.back() ;
		string line ;
		if( !getline( s, line ) ) break ;
		stringstream ss( line ) ;
		ss >> type >> hsa_chr >> hsa_strand_code >> r.hsa.pos >> hsa_end
		   >> ptr_chr >> ptr_strand_code >> r.ptr.pos >> ptr_end 
		   >> r.sequence >> flags ;

		r.hsa.chr = lookup_sym( hsa_chr ) ;
		r.ptr.chr = lookup_sym( ptr_chr ) ;

		r.hsa.strand = hsa_strand_code == '+' ;
		r.ptr.strand = ptr_strand_code == '+' ;

		if( r.sequence.size() > 255 )
			error( 1, 0, "cannot represent long insert at %s:%d (%d)",
					hsa_chr.c_str(), r.hsa.pos, r.sequence.size() ) ;

		if( type == "deletion" )
		{
			if( r.hsa.pos != hsa_end )
				error( 1, 0, "parse error: hsa should have length 0 for deletion at %s:%d", hsa_chr.c_str(), r.hsa.pos ) ;

			// store (chimp) base, only for sanity check
			r.ptr.base = r.sequence.size() == 1 ? maybe_compl( r.ptr.strand, to_ambicode( r.sequence[0] ) ) : 0 ;
			r.hsa.base = 0 ;
		}
		else if( type == "insert" )
		{
			if( r.ptr.pos != ptr_end )
				error( 1, 0, "parse error: ptr should have length 0 for insert at %s:%d", ptr_chr.c_str(), r.ptr.pos ) ;

			// store (human) base, only for sanity check
			r.hsa.base = r.sequence.size() == 1 ? maybe_compl( r.hsa.strand, to_ambicode( r.sequence[0] ) ) : 0 ;
			r.ptr.base = 0 ;
		}
		else error( 1, 0, "parse error: unknown indel type %s at %s:%d", type.c_str(), ptr_chr.c_str(), r.ptr.pos ) ;

		r.hsa.length = hsa_end - r.hsa.pos ;
		r.ptr.length = ptr_end - r.ptr.pos ;

		if( !flags.empty() ) error( 0, 0, "Flags field not empty at %s:%d", hsa_chr.c_str(), r.hsa.pos ) ;
		r.hsa.edge_near_flag = 0 ;
		r.ptr.edge_near_flag = 0 ;
		r.hsa.gap_near_flag = 0 ;
		r.ptr.gap_near_flag = 0 ;
		r.hsa.seen = 0 ;
		r.ptr.seen = 0 ;
	}
	delete d.back() ;
	d.pop_back() ;
	return s ;
}

inline ostream& write_half_record_snp( ostream& s, const SnpRec1& r )
{
	return s 
		<< from_ambicode( maybe_compl( r.strand, r.base ) ) << '\t' 
		<< symbols[r.chr] << '\t' << (r.strand ? '+' : '-') << '\t'
		<< r.pos << '\t' ;
}
inline ostream& write_half_record_indel( ostream& s, const SnpRec1& r )
{
	return s 
		<< symbols[r.chr] << '\t' << (r.strand ? '+' : '-') << '\t'
		<< r.pos << '\t' << r.pos + r.length << '\t' ;
}
inline ostream& write_bases( ostream& s, const SnpRec1& r )
{
	s << '\t' ;
	// if( r.seen && r.nt_bases.empty() ) s << '-' ;
	if( r.seen ) s << '"' << r.nt_bases << '"' << '\t' << r.nt_qual ;
	else s << "n/a\tn/a" ;
	return s ;
}
inline const char* encode_flags( const SnpRec& r )
{
	bool g = r.hsa.gap_near_flag || r.ptr.gap_near_flag ;
	bool e = r.hsa.edge_near_flag || r.ptr.edge_near_flag ;
	return g ? e ? "GC:EC" : "GC" : e ? "EC" : "" ;
}

template< typename C > ostream& write_martin_table_snp( ostream& s, const C& d )
{
	Chan ch ;
	s << "#HSA_Base\tHSA_Chr\tHSA_Strand\tHSA_Pos\tPAN_Base\tPAN_Chr\tPAN_Strand\tPAN_Pos\tOutBase\tFlag\tNEA_BaseH\tNEA_QualH\tNEA_BaseC\tNEA_QualC\n" ;
	for( typename C::const_iterator i = d.begin() ; i != d.end() ; ++i )
	{
		if( (i-d.begin()) % 1024 == 0 )
		{
			stringstream s ;
			s << "writing SNPs: " << i-d.begin() << '/' << d.size() ;
			ch( Console::info, s.str() ) ;
		}

		const SnpRec &r = **i ;
		if( !r.sequence.empty() ) continue ; // this wasn't actually a SNP
		// if( !r.hsa_seen && !r.ptr_seen ) continue ;
		
		write_half_record_snp( s, r.hsa ) ;
		write_half_record_snp( s, r.ptr ) ;
		s << r.out_base << '\t' << encode_flags(r) ;
		write_bases( s, r.hsa ) ;
		write_bases( s, r.ptr ) ;
		s << '\n' ;
	}
	return s ;
}

template< typename C > ostream& write_martin_table_indel( ostream& s, const C& d )
{
	Chan ch ;
	s << "#Type\tHSA_Chr\tHSA_Strand\tHSA_Start\tHSA_End\tPAN_Chr\tPAN_Strand\tPAN_Start\tPAN_End\tSeq\tFlag\tNEA_SeqH\tNEA_QualH\tNEA_SeqC\tNEA_QualC\n" ;
	for( typename C::const_iterator i = d.begin() ; i != d.end() ; ++i )
	{
		if( (i-d.begin()) % 1024 == 0 )
		{
			stringstream s ;
			s << "writing indels: " << i-d.begin() << '/' << d.size() ;
			ch( Console::info, s.str() ) ;
		}

		const SnpRec &r = **i ;
		if( r.sequence.empty() ) continue ; // this wasn't actually an indel
		// if( !r.hsa_seen && !r.ptr_seen ) continue ;
		
		s << (r.hsa.length ? "insert\t" : "deletion\t") ;
		write_half_record_indel( s, r.hsa ) ;
		write_half_record_indel( s, r.ptr ) ;
		s << r.sequence << '\t' << encode_flags(r) ;
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
	int k = 0 ;
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
			while( first_snp != mt.end() && get( first_snp ).chr != cur_chr 
					&& symbols[ get( first_snp ).chr ] < h.sequence() )
				++first_snp ;

			if( ++k % 1024 == 0 )
			{
				stringstream s ;
				s   << "SNPs vs. " << genome << " (" << h.sequence() << '/' 
					<< symbols[ get( first_snp ).chr ] << "): " 
					<< first_snp - mt.begin() << '/' << mt.size() ;
				progress( Console::info, s.str() ) ;
			}

			if( get( first_snp ).chr != cur_chr ) continue ;

			CompactGenome &g = Metagenome::find_sequence( h.genome_name(), h.sequence(), Metagenome::ephemeral ) ;
			DnaP ref = g.find_pos( h.sequence(), h.start_pos() ) ;
			int cigar_maj = 0, ref_pos = h.start_pos() ;
			size_t cigar_min = 0, qry_pos = 0 ;

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
				// (keep gap_buffer in mind, we need to look at a bit more)
				while( first_snp != mt.end() && get( first_snp ).chr == cur_chr
						&& get( first_snp ).pos + get( first_snp ).length + gap_buffer <= ref_pos )
					++first_snp ;

				// bail if we left the current chromosome or no SNPs are left at all
				if( first_snp == mt.end() || get( first_snp ).chr != cur_chr ) break ;

				// iterate over SNPs sufficiently close to current position
				for( vector<SnpRec*>::iterator cur_snp = first_snp ; cur_snp != mt.end() ; ++cur_snp )
				{
					SnpRec1 &snp = get( cur_snp ) ;
					Hit::Operation op = cigar_op( h.cigar(cigar_maj) ) ;
					if( snp.chr != cur_chr || snp.pos - gap_buffer >= ref_pos ) break ;

					// sanity check: only possible if we're not looking at an insert
					if( snp.pos == ref_pos && snp.base && *ref != snp.base ) switch( op )
					{
						case Hit::Delete:
						case Hit::Match:
						case Hit::Mismatch:
							error( 1, 0, "at %d, %c != %c -- wrong coordinate system?", 
									ref_pos, from_ambicode(*ref), from_ambicode( snp.base ) ) ;
						default:
							break ;
					}
					

					// Anything to extract? (coordinates hit and not a deletion
					if( snp.pos <= ref_pos && ref_pos < snp.pos + snp.length ) switch( op )
					{
						case Hit::Insert:
						case Hit::Match:
						case Hit::Mismatch:
							{
								uint8_t q = res.read().quality()[qry_pos] ;
								snp.nt_bases.push_back( res.read().sequence()[qry_pos] ) ;
								snp.nt_qual.push_back( q < 126-33 ? q+33 : 126 ) ;
							}
						default:
							break ;
					}
					else switch( op )
					{
						// nodirect hit, but close by.  Set gap flag?
						case Hit::Delete:
						case Hit::Insert:
							snp.gap_near_flag = 1 ;
						default:
							break ;
					}

					// SNP observed?  (hit coordinates or observed left and right adjacent positions)
					if( (snp.pos <= ref_pos && ref_pos < snp.pos + snp.length) ||
							(ref_pos == snp.pos && snp.length == 0 && qry_pos > 0) )
					{
						// close to contig edge? if so, set flag
						if( (int)qry_pos < gap_buffer || qry_pos >= res.read().sequence().size() - gap_buffer )
							snp.edge_near_flag = 1 ;
						snp.seen = 1 ;
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
	char *ptr_file = 0 ;
	char *hsa_file = 0 ;
	char *snp_out_file = 0 ;
	char *indel_out_file = 0 ;
	int core_limit = 0 ;

	console.loglevel = Console::debug ;

	struct poptOption options[] = {
		{ "hsa-file",     0 , POPT_ARG_STRING, &hsa_file,           0,    "Neandertalized Human is in FILE", "FILE" },
		{ "ptr-file",     0 , POPT_ARG_STRING, &ptr_file,           0,    "Neandertalized Chimp is in FILE", "FILE" },
		{ "output-snp",   0 , POPT_ARG_STRING, &snp_out_file,       0,    "Write SNPs table to FILE", "FILE" },
		{ "output-indel", 0 , POPT_ARG_STRING, &indel_out_file,     0,    "Write indels table to FILE", "FILE" },
		{ "buffer",       0 , POPT_ARG_INT,    &gap_buffer,         0,    "Flag gaps closer than N", "N" },
		{ "vmem",         0 , POPT_ARG_INT,    &core_limit,         0, "limit virtual memory to X megabytes", "X" },
		POPT_AUTOHELP POPT_TABLEEND
	} ;

	poptContext pc = poptGetContext( "martin-table", argc, argv, options, 0 ) ;
	poptSetOtherOptionHelp( pc, "[OPTION...] [martin-table...]" ) ;
	poptReadDefaultConfig( pc, 0 ) ;

	int rc = poptGetNextOpt( pc ) ; 
	if( rc != -1 ) error( 1, 0, "%s %s", poptStrerror( rc ), poptBadOption( pc, 0 ) ) ;

	if( !poptPeekArg( pc ) ) error( 1, 0, "no input files (try --help)" ) ;

	if( core_limit ) {
		struct rlimit lim ;
		getrlimit( RLIMIT_AS, &lim ) ;
		lim.rlim_cur = 1024*1024 * (long)core_limit ;
		setrlimit( RLIMIT_AS, &lim ) ;
	}

	std::vector<SnpRec*> mt ;
	while( char const *arg = poptGetArg( pc ) )
	{
		ifstream f( arg ) ;
		if( strstr( arg, "SNP" ) ) read_martin_table_snp( f, mt, arg ) ;
		else if( strstr( arg, "indel" ) ) read_martin_table_indel( f, mt, arg ) ;
		else error( 1, errno, "cannot guess contents of %s", arg ) ;

		if( !f.eof() ) error( 1, errno, "Parse error in 'Martin table' %s.", arg ) ;
	}

	// cerr << mallinfo().arena << endl ;

	if( !ptr_file ) error( 1, 0, "missing Chimp alignments" ) ;
	if( !hsa_file ) error( 1, 0, "missing Human alignments" ) ;
	if( !snp_out_file && !indel_out_file ) error( 1, 0, "no output to write to" ) ;

	console.output( Console::info, "Sorting on pt2..." ) ;
	sort( mt.begin(), mt.end(), ByPt2Coordinate() ) ;
	console.output( Console::info, "Done." ) ;
	scan_anfo_file( mt, ptr_file, "pt2", GetPt2Rec() ) ;

	console.output( Console::info, "Sorting on hg18..." ) ;
	sort( mt.begin(), mt.end(), ByHg18Coordinate() ) ;
	console.output( Console::info, "Done." ) ;
	scan_anfo_file( mt, hsa_file, "hg18", GetHg18Rec() ) ;

	if( snp_out_file ) {
		ofstream out( snp_out_file ) ;
		write_martin_table_snp( out, mt ) ;
	}
	if( indel_out_file ) {
		ofstream out( indel_out_file ) ;
		write_martin_table_indel( out, mt ) ;
	}
	poptFreeContext( pc ) ;
	return 0 ;
}

