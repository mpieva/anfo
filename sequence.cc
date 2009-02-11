#include <sequence.h>
#include <cmath>
#include <limits>
#include <sstream>

using namespace std ;

static inline bool seq_continues( std::istream& s )
{ return s && s.peek() != '@' && s.peek() != '>' && s.peek() != '+' ; }

static inline bool descr_follows( std::istream& s )
{ return s && s.peek() == ';' ; }

static inline bool all_acsii_qscores( const std::string& s )
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

static inline int sol_to_phred( int sol )
{
	return (int)( 0.5 + 10.0/log(10.0) * log( 1.0 + exp( sol*log(10.0)/10.0 ))) ;
}

istream& read_fastq( istream& s, QSequence& qs, bool solexa_scores )
{
	qs.seq_.clear() ;
	qs.seq_.push_back( 0 ) ;
	// skip junk before sequence header
	while( seq_continues(s) ) s.ignore( std::numeric_limits<int>::max(), '\n' ) ;

	// header follows, don't care for the delimiter, read name and
	// description
	s.get() ;
	s >> qs.name_ ;
	getline( s, qs.description_ ) ;

	// If at this point we have a valid stream, we definitely have a
	// sequence.  Bail out if reading the header failed.
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
		for( size_t i = 0 ; i != line.size() ; ++i )
			if( encodes_nuc( line[i] ) )
				qs.seq_.push_back( 0x2800 | to_ambicode( line[i] ) ) ;
	}
	qs.seq_.push_back( 0 ) ;

	// if quality follows...
	if( s && s.peek() == '+' )
	{
		// skip delimiter, name, and description no additional
		// description lines can follow, since ';' is a valid Q-score
		s.ignore( std::numeric_limits<int>::max(), '\n' ) ;

		// Q-scores must follow unless the sequence was empty or the stream ends
		if( s && qs.length() )
		{
			qs.has_quality_ = true ;

			string line ;
			getline( s, line ) ;

			// Check one line; if it contains only spaces and numbers in
			// groups of no more than three with an optional sign, it's
			// an Alta Cyclic or U Rockefeller file and we decode ASCII
			// numbers until the next header
			if( all_acsii_qscores(line) ) 
			{
				for( int ix = 0 ; s ; )
				{
					stringstream ss( line ) ;
					for( int q = 0 ; ss >> q ; ++ix )
						qs.qual( ix, solexa_scores ? sol_to_phred(q) : q ) ;
					if( !seq_continues(s) ) break ;
					getline( s, line ) ;
				}
			}
			// No ASCII coding.  Since delimiters aren't very useful
			// now, we'll take exactly one Q-score for each nucleotide,
			// ignoring line feeds.
			else
			{
				size_t total = qs.length(), ix = 0 ; 
				while( ix != total && s )
				{
					for( size_t j = 0 ; ix != total && j != line.size() ; ++j, ++ix )
					{
						int q = line[j] ;
						qs.qual( ix, solexa_scores ? sol_to_phred( q-64 ) : q-33 ) ;
					}
					if( ix != total ) getline( s, line ) ;
				}
				// There might be some junk left over; it will be eaten
				// away in the next call.
			}
		}
	}

	// We did get a sequence, no matter the stream state now, so no
	// failure.  If we reached EOF, the flags must remain.
	s.clear( s.rdstate() & ~istream::failbit ) ;
	return s ;
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

