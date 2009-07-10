#include "sequence.h"
#include "logdom.h"
#include <cmath>
#include <limits>
#include <sstream>

using namespace std ;

namespace {

//! \brief test whether a line contains ASCII Q-scores
//! \internal
//! This is a small DFA that test whether a line contains only small
//! ACSII encoded numbers separated by spaces.
inline bool all_ascii_qscores( const std::string& s )
{
	int st = 0 ;
	for( size_t i = 0 ; i != s.length() ; ++i )
	{
		if( isspace(s[i]) ) st = 0 ;
		else if( st == 0 && s[i] == '-' ) st = 1 ;
		else if( st == 0 && isdigit(s[i]) ) st = 2 ;
		else if( st == 1 && isdigit(s[i]) ) st = 2 ;
		else if( st == 2 && isdigit(s[i]) ) st = 3 ;
		else if( st == 3 && isdigit(s[i]) ) st = 4 ;
		else return false ;
	}
	return true ;
}

// inline double sol_to_err_prob( int sol ) { return 1.0 / ( 1.0 + std::pow( 10.0, sol / 10.0 ) ) ; }
inline double phred_to_err_prob( int phred ) { return std::pow( 10.0, -phred / 10.0 ) ; }
inline int err_prob_to_phred( double p ) { return (int)( -10.0 * std::log( p ) / std::log( 10.0 ) ) ; }
inline int sol_to_phred( int sol ) { return sol + (int)( 0.5 + 10 * log1p( pow( 10.0, -sol / 10.0 ) ) / log( 10.0 ) ) ; }

int bits_in[] = { 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 } ;

class Reader {
	private: 
		const char *buf_ ;
		int len_ ;
		google::protobuf::io::ZeroCopyInputStream *zis_ ;
	public:
		Reader( google::protobuf::io::ZeroCopyInputStream *zis ) : zis_(zis) {
			const void *p ;
			do {
				buf_ = zis_->Next( &p, &len_ ) ? (char*)p : 0 ;
			} while( buf_ && !len_ ) ;
		}
		~Reader() { if( len_ && buf_ ) zis_->BackUp( len_ ) ; }

		operator const void* () const { return buf_ ; }
		char peek() const { return buf_ ? *buf_ : 0 ; }
		char get() { 
			if( !buf_ ) return 0 ;

			const void *p ;
			char r = *buf_++ ;
			--len_ ;

			while( buf_ && !len_ ) {
				buf_ = zis_->Next( &p, &len_ ) ? (char*)p : 0 ;
			}
			return r ;
		}
} ;

inline void getword( Reader& r, std::string& s )
{
	s.clear() ;
	while( r && !isspace( r.peek() ) ) s.push_back( r.get() ) ;
	while( r &&  isspace( r.peek() && r.peek() != '\n' ) ) r.get() ;
}
inline void getline( Reader& r, std::string& s )
{
	s.clear() ;
	while( r && r.peek() != '\n' ) s.push_back( r.get() ) ;
	if( r ) r.get() ;
}
inline void skipline( Reader& r) { while( r && r.get() != '\n' ) ; }

inline bool seq_continues( const Reader &s ) { return s && s.peek() != '@' && s.peek() != '>' && s.peek() != '+' && s.peek() != '*' ; }
inline bool descr_follows( const Reader &s ) { return s && s.peek() == ';' ; } 
}

QSequence::Base::Base( uint8_t a, int q ) : ambicode( a ), qscore( q )
{
	int bits = bits_in[ a ] ;
	Logdom prob = Logdom::from_phred( q ) ;
	for( int i = 0 ; i != 4 ; ++i )
	{
		// qualities[i] = a & (1<<i) ? (1.0-prob) / bits : prob / (4-bits) ;
		qscores[i] = a & (1<<i) ? ( (Logdom::from_float(1.0)-prob) / Logdom::from_float(bits) ).to_phred_byte()
			                    : ( prob / Logdom::from_float(4-bits) ).to_phred_byte() ;
		// qscores[i] = a & (1<<i) ? err_prob_to_phred( (1.0-prob)/bits )
		//	                    : q + (uint8_t)(0.5 + 10*std::log(4-bits)/std::log(10)) ;
	}
}

QSequence::QSequence( const char* p, int q_score ) 
		: seq_(), name_(), description_(), validity_( bases_only )
{
	seq_.push_back( Base() ) ;
	for( ; *p ; ++p )
		if( encodes_nuc( *p ) ) 
			seq_.push_back( Base( to_ambicode( *p ), q_score ) ) ;
	seq_.push_back( Base() ) ;
}

QSequence::QSequence( const char* p, const uint8_t* q, const std::string& name, const std::string &descr )
		: seq_(), name_(name), description_(descr), validity_( bases_with_quality )
{
	seq_.push_back( Base() ) ;
	for( ; *p ; ++p, ++q )
		if( encodes_nuc( *p ) ) 
			seq_.push_back( Base( to_ambicode( *p ), *q ) ) ;
	seq_.push_back( Base() ) ;
}
					
bool read_fastq( google::protobuf::io::ZeroCopyInputStream *zis, QSequence& qs, bool solexa_scores, char origin )
{
	Reader s( zis ) ;
	bool got_seq = false, got_qual = false, got_quals = false ;

	qs.seq_.clear() ;
	qs.seq_.push_back( QSequence::Base() ) ;
	// skip junk before sequence header
	while( seq_continues(s) ) skipline(s) ;

	// header follows, don't care for the delimiter, read name and
	// description
	s.get() ;
	getword( s, qs.name_ ) ;
	getline( s, qs.description_ ) ;

	// If at this point we have a valid stream, we definitely have a
	// sequence.  Bail out iff reading of the header failed.
	if( !s ) return s ;

	// if more description follows, read it in, dropping the delimiter
	while( descr_follows(s) ) 
	{
		s.get() ;
		string line ;
		getline( s, line ) ;
		qs.description_.push_back( '\n' ) ;
		qs.description_ += line ;
	}

	// read sequence while it continues
	// do not read gaps, ignore junk
	while( seq_continues(s) )
	{
		string line ;
		getline( s, line ) ;
		// If line starts with "SQ ", we ignore it.  Scan the rest for
		// nucleotide codes.
		for( size_t i = line.substr(0,3) == "SQ " ? 3 : 0 ; i != line.size() ; ++i )
		{
			if( encodes_nuc( line[i] ) )
			{
				got_seq = true ;
				qs.seq_.push_back( QSequence::Base( to_ambicode( line[i] ), 30 ) ) ;
			}
		}
	}
	qs.seq_.push_back( QSequence::Base() ) ;

	// if quality follows...
	if( s && s.peek() == '+' )
	{
		// skip delimiter, name, and description no additional
		// description lines can follow, since ';' is a valid Q-score
		skipline(s) ;

		// Q-scores must follow unless the sequence was empty or the stream ends
		if( s && qs.length() )
		{
			got_qual = true ;

			string line ;
			getline( s, line ) ;

			// Check one line; if it contains only spaces and numbers in
			// groups of no more than three with an optional sign, it's
			// an Alta Cyclic or U Rockefeller file and we decode ASCII
			// numbers until the next header
			if( all_ascii_qscores(line) ) 
			{
				for( int ix = 1 ; s ; )
				{
					stringstream ss( line ) ;
					for( int q = 0 ; ss >> q ; ++ix )
						qs.seq_[ix] = QSequence::Base(
								qs.seq_[ix].ambicode, solexa_scores ? sol_to_phred(q) : q
								) ; 
					if( !seq_continues(s) ) break ;
					getline( s, line ) ;
				}
			}
			// No ASCII coding.  Since delimiters aren't very useful
			// now, we'll take exactly one Q-score for each nucleotide,
			// ignoring LF and CR.
			else
			{
				size_t total = qs.length(), ix = 0 ; 
				while( ix != total && s )
				{
					for( size_t j = 0 ; ix != total && j != line.size() ; ++j )
					{
						int q = line[j] ;
						if( q != 13 ) // skip CRs
						{
							qs.seq_[ix+1] = QSequence::Base(
									qs.seq_[ix+1].ambicode, solexa_scores ? sol_to_phred( q-origin ) : q-origin ) ; 
							++ix ;
						}
					}
					if( ix != total ) getline( s, line ) ;
				}
				// There might be some junk left over; but it will be
				// eaten away in the next call.
			}
		}
	}
#if 0
	else
	{
		size_t pos[4] = {1,1,1,1} ;

		// We might get quality in 4Q format.  We'll handle it as
		// any number of optional quality lines, so FASTA degerates to a
		// special case of 4Q.
		while( s && s.peek() == '*' ) 
		{
			s.get() ;	// drop the star
			int tag = s.get() & ~32 ;
			tag = tag == 'A' ? 0 : tag == 'C' ? 1 : tag == 'T' ? 2 : tag == 'G' ? 3 : -1 ;
			while( s && isspace(s.peek()) && s.peek() != '\n' ) s.get() ;

			while( s && s.peek() != '\n' ) 
			{
				int q = s.get() ;
				// skip lines with unrecognized tag, skip CRs
				if( tag != -1 && q != 13 ) 
				{
					got_quals = true ;
					if( pos[tag] == qs.seq_.size()-1 ) qs.seq_.push_back( QSequence::Base() ) ;
					qs.seq_[pos[tag]].qualities[tag] = phred_to_err_prob( q-origin ) ;
					qs.seq_[pos[tag]].qscores[tag] = q-origin ;
				}
			}
		}
		for( size_t p = qs.seq_.size()-2 ; p>0 ; --p )
		{
			// simple basecall, in case we have Q scores, but no
			// sequence
			if( !qs.seq_[p].ambicode ) qs.seq_[p].ambicode = 1 << (
						std::max_element( qs.seq_[p].qualities, qs.seq_[p].qualities+4 ) 
						- qs.seq_[p].qualities ) ;
			// generate single Q-Score (take minimum, that's a good
			// approximation) X X X probably wrong
			qs.seq_[p].qscore = std::min( std::min( qs.seq_[p].qscores[0], qs.seq_[p].qscores[1] ),
					                      std::min( qs.seq_[p].qscores[2], qs.seq_[p].qscores[3] ) ) ;
		}
	}
#endif

	qs.validity_ = got_quals ? got_seq ? QSequence::bases_with_qualities : QSequence::qualities_only
		                     : got_qual ? QSequence::bases_with_quality : QSequence::bases_only ;
	// We did get a sequence, no matter the stream state now, so no
	// failure.
	return true ;
}

#if BUILD_TEST_HARNESS
// very stupid test harness, not normally compiled

#include <iostream>

int main()
{
	QSequence qs ;
	while( read_fastq( cin, qs ) ) {
		cout << qs.get_name() << endl 
			<< qs.get_descr() << endl ;
		for( int i = 0 ; i != qs.length() ; ++i )
			cout << (int)qs.qual(i) << ' ' ;
		cout << endl ;
	}
}
#endif

