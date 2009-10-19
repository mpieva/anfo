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

inline bool seq_continues( const Reader &s ) { char c = s.peek() ; return s && c != '@' && c != '>' && c != '+' && c != '*' ; }
inline bool descr_follows( const Reader &s ) { return s && s.peek() == ';' ; } 
inline void skipline     (       Reader &r ) { while( r && r.get() != '\n' ) ; }
}

QSequence::Base::Base( uint8_t a, int q ) : ambicode( a ), qscore( q )
{
	int bits = bits_in[ a ] ;
	Logdom prob = Logdom::from_phred( q ) ;
	for( int i = 0 ; i != 4 ; ++i )
	{
		qscores[i] = a & (1<<i) ? ( (1-prob) / Logdom::from_float(bits) ).to_phred_byte()
			                    : ( prob / Logdom::from_float(4-bits) ).to_phred_byte() ;
	}
}

QSequence::QSequence( const output::Read& r, int default_q )
{
	std::string::const_iterator p = r.sequence().begin() + r.trim_left(),
		pe = r.has_trim_right() ? r.sequence().begin() + r.trim_right() : r.sequence().end() ;

	seq_.push_back( Base() ) ;
	if( r.has_quality() ) 
	{
		std::string::const_iterator q = r.quality().begin() + r.trim_left(),
			qe = r.has_trim_right() ? r.quality().begin() + r.trim_right() : r.quality().end() ;
		for( ; p != pe ; ++p, ++q )
			seq_.push_back( Base( to_ambicode( *p ), *q ) ) ;
	}
	else for( ; p != pe ; ++p )
			seq_.push_back( Base( to_ambicode( *p ), default_q ) ) ;

	seq_.push_back( Base() ) ;
}
					
bool read_fastq( google::protobuf::io::ZeroCopyInputStream *zis, output::Read& r, bool solexa_scores, char origin )
{
	Reader s( zis ) ;

	// skip junk before sequence header
	while( seq_continues(s) ) skipline(s) ;

	// header follows, don't care for the delimiter, read name and
	// description
	s.get() ;
	getword( s, *r.mutable_seqid() ) ;
	getline( s, *r.mutable_description() ) ;
	r.clear_sequence() ;

	// If at this point we have a valid stream, we definitely have a
	// sequence.  Bail out iff reading of the header failed.
	if( !s ) return s ;

	// if more description follows, read it in, dropping the delimiter
	while( descr_follows(s) ) 
	{
		s.get() ;
		string line ;
		getline( s, line ) ;
		r.mutable_description()->push_back( '\n' ) ;
		*r.mutable_description() += line ;
	}

	// read sequence while it continues
	// do not read gaps, ignore junk
	while( seq_continues(s) )
	{
		string line ;
		getline( s, line ) ;
		// If line starts with "SQ ", we ignore it (4Q support).  Scan
		// the rest for nucleotide codes.
		for( size_t i = line.substr(0,3) == "SQ " ? 3 : 0 ; i != line.size() ; ++i )
		{
			if( encodes_nuc( line[i] ) )
			{
				r.mutable_sequence()->push_back( line[i] ) ;
			}
		}
	}

	// if quality follows...
	if( s && s.peek() == '+' )
	{
		// skip delimiter, name, and description no additional
		// description lines can follow, since ';' is a valid Q-score
		skipline(s) ;

		// Q-scores must follow unless the sequence was empty or the stream ends
		if( s && r.sequence().length() )
		{
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
						r.mutable_quality()->push_back( solexa_scores ? sol_to_phred(q) : q ) ;
					if( !seq_continues(s) ) break ;
					getline( s, line ) ;
				}
			}
			// No ASCII coding.  Since delimiters aren't very useful
			// now, we'll take exactly one Q-score for each nucleotide,
			// ignoring LF and CR.
			else
			{
				size_t total = r.sequence().length(), ix = 0 ; 
				while( ix != total && s )
				{
					for( size_t j = 0 ; ix != total && j != line.size() ; ++j )
					{
						int q = line[j] ;
						if( q != 13 ) // skip CRs
						{
							r.mutable_quality()->push_back( solexa_scores ? sol_to_phred( q-origin ) : q-origin ) ; 
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

	// We did get a sequence, no matter the stream state now, so no
	// failure.
	return true ;
}
