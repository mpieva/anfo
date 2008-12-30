#include <sequence.h>
#include <sstream>

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

std::istream& read_fastq( std::istream& s, QSequence& qs )
{
	qs.seq.clear() ;
	qs.seq.push_back( 0 ) ;
	while( seq_continues(s) ) s.ignore( INT_MAX, '\n' ) ;

	// header follows, don't care for the delimiter, read name and
	// description
	s.get() ;
	s >> qs.name ;
	std::getline( s, qs.description ) ;

	// If at this point we have a valid stream, we definitely have a
	// sequence.  Bail out if reading the header failed.
	if( !s ) return s ;

	// if more description follows, read it in, dropping the delimiter
	while( descr_follows(s) ) 
	{
		s.get() ;
		std::string line ;
		getline( s, line ) ;
		qs.description.push_back( '\n' ) ;
		qs.description += line ;
	}

	// read sequence while it continues
	while( seq_continues(s) )
	{
		std::string line ;
		getline( s, line ) ;
		for( size_t i = 0 ; i != line.size() ; ++i )
			qs.seq.push_back( 0x2800 | to_ambicode( line[i] ) ) ;
	}
	qs.seq.push_back( 0 ) ;

	// if quality follows...
	// XXX any hint how to see if it's a Solexa file?
	if( s && s.peek() == '+' )
	{
		// don't care for the delimiter, but get name and description
		// (again)
		s.get() ;
		std::string descr_, name_ ;
		if( s.peek() != '\n' ) s >> name_ ;
		std::getline( s, descr_ ) ;

		// if more description follows, read it in, dropping the delimiter
		while( descr_follows(s) ) 
		{
			s.get() ;
			std::string line ;
			getline( s, line ) ;
			descr_.push_back( '\n' ) ;
			descr_ += line ;
		}

		// were name and descr. repeated?  might be a useful hint...
		bool rep_name = name_ == qs.name ;
		bool rep_descr = descr_ == qs.description ;

		if( seq_continues(s) ) {
			std::string line ;
			getline( s, line ) ;

			// Check one line; if it contains only spaces and numbers in
			// groups of no more than three with an optional sign, it's
			// an Alta Cyclic or U Rockefeller file and we decode ACSII
			// numbers until the next header
			if( all_acsii_qscores(line) ) 
			{
				for( int ix = 0 ; s ; )
				{
					std::stringstream ss( line ) ;
					// XXX assume phred quality scores.  Need a
					// solution for Solexa nonsense...
					for( int q = 0 ; ss >> q ; ++ix ) qs.qual( ix, q ) ;
					if( !seq_continues(s) ) break ;
					getline( s, line ) ;
				}
			}
			// No ASCII coding.  Since delimiters aren't very useful
			// now, we'll take exactly one Q-score for each nucleotide,
			// ignoring line feeds.
			else
			{
				size_t total = qs.length() ; 
				for( size_t ix = 0 ; s && ix != total ; ++ix )
				{
					for( size_t j = 0 ; ix != total && j != line.size() ; ++j, ++ix )
					{
						// XXX assume phred quality scores.  Need a
						// solution for Solexa nonsense...
						uint8_t q = line[j] - 33 ;
						qs.qual( ix, q ) ;
					}
					if( !seq_continues(s) ) break ;
					getline( s, line ) ;
				}
				// There might be some junk left over; it will be eaten
				// away in the next call.
			}
		}
	}

	// We did get a sequence, no matter the stream state now, so no
	// failure.  If we reached EOF, the flags must remain.
	s.clear( s.rdstate() & ~std::istream::failbit ) ;
	return s ;
}

#if 0
// very stupid test harness, not normally compiled

#include <iostream>

int main()
{
	QSequence qs ;
	while( read_fastq( std::cin, qs ) ) {
		std::cout << qs.get_name() << std::endl 
			<< qs.get_descr() << std::endl ;
		for( int i = 0 ; i != qs.length() ; ++i )
			std::cout << (int)qs.qual(i) << ' ' ;
		std::cout << std::endl ;
	}
}
#endif

