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

#ifndef INCLUDED_STREAM_H
#define INCLUDED_STREAM_H

#include "compress_stream.h"
#include "logdom.h"
#include "output.pb.h"
#include "sequence.h"
#include "util.h"

#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>

#include <algorithm>
#include <deque>
#include <fstream>
#include <memory>
#include <vector>

#include <sys/types.h>
#include <sys/wait.h>

// Including Elk defines some macros that collide with protobuf.  We
// undefine them (and hope they aren't needed...).

#if HAVE_ELK_SCHEME_H
#include <elk/scheme.h>
#undef Print
#undef MAX_TYPE
#endif

namespace streams {
	using namespace output ;
	using namespace std ;

inline uint32_t mk_msg_tag( uint32_t i ) { return i << 3 | 2 ; }

//! \defgroup outputfile Convenience functions to handle output files
//! The output file is a series of protocol buffer messages.  This makes
//! it compact and extensible.  Instead of writing a serial file, the
//! very same messages could be entered into a Berkeley-DB type
//! database.  Since conversion to all sorts of text files is a no
//! brainer, the binary format is not considered undue hardship on the
//! users of scripting languages.
//!
//! (If you are a user of a glorified logfile munching language and want
//! an ASCII-based table, define a sensible layout and ask for it.  Or
//! write the single function yourself, it isn't that hard.)
//!
//! The protobuf format itself is not self-delimiting, therefore we
//! prepend each message with the number of bytes it occupies, in
//! variable protobuf format.  The first message in a file
//! one is of type output::Header and should contain all sorts of
//! meta information, maybe even configuration of the mapper that was
//! used to create it (an option would be to include the full set of
//! policies).  Subsequent messages are of type output::Result.
//!
//! \note The output file is essentially a single, giant protobuf message.
//!       However, it is impossible to read it in one go (due to memory
//!       constraints), but it's also impossible to read it through a single
//!       CodedInputStream.  That's why methods in here tend to construct and
//!       destruct a fresh CodedInputStream on each invocation.
//!
//! @{

template< typename Msg >
void write_delimited_message( google::protobuf::io::CodedOutputStream& os, int tag, const Msg& m )
{
	os.WriteTag( mk_msg_tag( tag ) ) ;
	os.WriteVarint32( m.ByteSize() ) ;
	if( !m.SerializeToCodedStream( &os ) )
		throw "error while serializing" ;
}

template< typename Msg >
void write_delimited_message( google::protobuf::io::ZeroCopyOutputStream *os, int tag, const Msg& m )
{
	google::protobuf::io::CodedOutputStream o( os ) ;
	write_delimited_message( o, tag, m ) ;
}

template< typename Msg >
bool read_delimited_message( google::protobuf::io::CodedInputStream& is, Msg &m )
{
	uint32_t size ;
	std::string code ;
	if( !is.ReadVarint32( &size ) ) return false ;
	int lim = is.PushLimit( size ) ;
	if( !m.ParseFromCodedStream( &is ) ) throw "error while deserializing" ;
	is.PopLimit( lim ) ;
	return true ;
}

void sanitize( Header& ) ;
void sanitize( Read& ) ;
void merge_sensibly( output::Header& lhs, const output::Header& rhs ) ;
void merge_sensibly( output::Footer& lhs, const output::Footer& rhs ) ;
void merge_sensibly( output::Result& lhs, const output::Result& rhs ) ;


//! \brief returns the hit to some genome
//! If an empty genome is asked for, returns the best hit.  Behaviour is
//! undefined if no suitable hit exists.
const output::Hit* hit_to( const output::Result& ) ;
const output::Hit* hit_to( const output::Result&, const string& ) ;
inline const output::Hit* hit_to( const output::Result& r, const char* g ) { return g ? hit_to( r ) : hit_to( r, string(g) ) ; }

template< typename I > const output::Hit* hit_to( const output::Result& r, I begin, I end )
{
	for( ; begin != end ; ++begin )
		if( const output::Hit* h = hit_to( r, *begin ) ) return h ;
	return 0 ;
}

//! \brief returns the mutable hit to some genome
//! If an empty genome is asked for, returns the best hit.  If no
//! suitable hit exists, a new one is created.
output::Hit* mutable_hit_to( output::Result* ) ;
output::Hit* mutable_hit_to( output::Result*, const string& ) ;

//! \brief computes (trimmed) query length from CIGAR line
template< typename C > unsigned len_from_bin_cigar( const C& cig )
{
	unsigned l = 0 ;
	for( typename C::const_iterator i = cig.begin(), e = cig.end() ; i != e ; ++i )
	{
		switch( cigar_op( *i ) )
		{
			case output::Hit::Match:
			case output::Hit::Mismatch:
			case output::Hit::Insert:
			case output::Hit::SoftClip:
				l += cigar_len( *i ) ;
				break ;

			case output::Hit::Delete:
			case output::Hit::Skip:
			case output::Hit::HardClip:
			case output::Hit::Pad:
				break ;
		}
	}
	return l ;
}

inline output::Hit::Operation cigar_op( uint32_t c ) { return (output::Hit::Operation)(c & 0xf) ; }
inline uint32_t cigar_len( uint32_t c ) { return c >> 4 ; }
inline uint32_t mk_cigar( output::Hit::Operation op, uint32_t len ) { return len << 4 | op ; }

inline void push_op( std::vector<unsigned>& s, unsigned m, output::Hit::Operation op )
{
	if( m && !s.empty() && (streams::cigar_op( s.back() ) == op) && streams::cigar_len( s.back() ) ) 
		s.back() = streams::mk_cigar( op, m + streams::cigar_len( s.back() ) ) ;
	else if( m ) s.push_back( streams::mk_cigar( op, m ) ) ;
}
inline void push_m( std::vector<unsigned>& s, unsigned m ) { push_op( s, m, output::Hit::Match ) ; }
inline void push_M( std::vector<unsigned>& s, unsigned m ) { push_op( s, m, output::Hit::Mismatch ) ; }
inline void push_i( std::vector<unsigned>& s, unsigned i ) { push_op( s, i, output::Hit::Insert ) ; }
inline void push_d( std::vector<unsigned>& s, unsigned d ) { push_op( s, d, output::Hit::Delete ) ; }

//! @}

//! \brief stream of result messages
//! Each stream is a header, followed by many results, followed by a
//! single footer.  The header will be cached internally, so it can be
//! asked for repeatedly.  Results are forgotten once read, the footer
//! is only available after all results have been read, and then it is
//! stored and can be read repeatedly.
//!
//! This class is both an input stream, an output stream and a stream
//! transducer.  Concrete implementations decide what sets of methods to
//! actually implement.  Streams behave as state machine, cycling
//! between need_input, have_output and end_of_stream.  It's undefined
//! behaviour to call methods in states where they don't make sense.

class Stream
{
	public:
		int refcount_ ;
		enum state { invalid, end_of_stream, need_input, have_output } ;

		static void cleanup( Stream* p ) { delete p ; }

	protected:
		// internal state---not strictly necessary here, but used almost
		// everywhere and therefore simply convenient
		Header hdr_ ;
		Result res_ ;
		Footer foot_ ;
		state state_ ;

		virtual ~Stream() {} 				// want control over instantiation

	private:
		Stream( const Stream& ) ; 			// must not copy
		void operator = ( const Stream& ) ;	// must not copy

	public:
		Stream() : refcount_(0), state_( end_of_stream ) {}

		//! \brief returns stream state
		//! The state determines which methods may be called: get_header
		//! can be called unless the satet is \c invalid, get_result can
		//! be called in state have_output, put_result only in
		//! need_input, get_footer only in end_of_stream.  put_footer
		//! can only be called in need_input and it will signal the end
		//! of the input stream (if there is such a thing).  put_header
		//! can only be called in state invalid, it may change the state
		//! to either need_input or have_output.  Particular streams may
		//! need special initialization before the state bacomes
		//! something other than invalid.
		virtual state get_state() { return state_ ; }

		//! \brief returns the header
		//! The header can be requested any time, unless the stream is
		//! in the invalid state.
		virtual Header fetch_header() { return hdr_ ; }

		//! \brief reads the next result
		//! Every result can only be read once, internal iterator style.
		//! If the stream state is have_output, this method returns the
		//! next result.  Otherwise the behavior is undefined.
		virtual Result fetch_result() { state_ = need_input ; return res_ ; }

		//! \brief returns the footer
		//! Only after all results have been consumed is the footer
		//! available.  If anything goes wrong internally, the LSB
		//! should be set in \c exit_code.
		virtual Footer fetch_footer() { return foot_ ; }

		//! \brief sets the stream header
		//! Output streams and stream filters need the header to become
		//! valid streams.
		virtual void put_header( const Header& h ) { hdr_ = h ; state_ = need_input ; }

		//! \brief outputs one result
		//! This method can only be called in state need_input, and it
		//! will insert a result record into the stream.  Else the
		//! behaviour is undefined.
		virtual void put_result( const Result& r ) { res_ = r ; state_ = have_output ; }

		//! \brief sets the footer
		//! This method can only be called in state need_input, and it
		//! signals that no more input is available.  A filter will then
		//! flush internal buffers and signal end_of_stream.
		virtual void put_footer( const Footer& f ) { foot_ = f ; state_ = end_of_stream ; }

		//! \brief get the 'summary' of having processed this stream
		//! This functionality is dependent on Elk being present:  since
		//! the meaning of what the 'summary' is differs from stream to
		//! stream, the result is simply an Elk object.  The default is
		//! to return the exit code contained in the footer.
#if HAVE_ELK_SCHEME_H
		virtual Object get_summary() const ;
#endif

	protected:
		// doesn't belong here, but it's convenient
		void read_next_message( google::protobuf::io::CodedInputStream&, const std::string& ) ;
} ;

typedef ::Holder<Stream> StreamHolder ;

int transfer( Stream& in, Stream& out ) ;

//! \brief base class of streams that read from many streams
class StreamBundle : public Stream
{
	protected:
		std::deque< StreamHolder > streams_ ;
		typedef std::deque< StreamHolder >::const_iterator citer ;
		typedef std::deque< StreamHolder >::const_reverse_iterator criter ;

	public:
		void add_stream( StreamHolder s ) { streams_.push_back( s ) ; }
#if HAVE_ELK_SCHEME_H
		virtual Object get_summary() const ;
#endif
} ;

//! \internal
//! Tracked to avoid bumping into the file descriptor limit
//! (mostly important for merge sorting and mega-merge).
extern int anfo_reader__num_files_ ;

struct ParseError : public Exception {
	std::string msg_ ;
	ParseError( const std::string& m ) : msg_(m) {}
	virtual void print_to( std::ostream& s ) const { s << "AnfoReader: " << msg_ ; }
} ;


//! \brief presents ANFO files as series of messages
//! This class will read a stream of results formatted as a continuous
//! protobuf message.
class AnfoReader : public Stream
{
	private:
		std::auto_ptr< google::protobuf::io::ZeroCopyInputStream > is_ ;
		std::string name_ ;

		virtual ~AnfoReader() { --anfo_reader__num_files_ ; }

	public: 
		AnfoReader( std::auto_ptr< google::protobuf::io::ZeroCopyInputStream > is, const std::string& name ) ;
		virtual Result fetch_result() ;

		//! \internal
		static unsigned num_open_files() { return anfo_reader__num_files_ ; }
} ;

// \brief reader for all supported formats
// Here we take care not to open files before the header is requested.
// This is necessary to allow merging thousands of files without
// directly running into a filedescriptor limit.
// To this end, the UniversalReader can be initialized with or without a
// stream object.  If the stream exists, we take care not to read from
// it until the header is needed, and the name given serves just for
// informational purposes.  If no stream exists, we create a
// FileInputStream from the name (which must be a filename, obviously)
// when the header is requested.  At this point we also inspect the
// stream to determine its format and create the appropriate filters to
// decode it.
class UniversalReader : public Stream
{
	private:
		std::auto_ptr< google::protobuf::io::ZeroCopyInputStream > is_ ;
		std::string name_ ;
		StreamHolder str_ ;

		bool solexa_scores_ ;
		int origin_ ;

	public: 
		UniversalReader(
				const std::string& name,
				google::protobuf::io::ZeroCopyInputStream* is = 0,
				bool solexa_scores = false,
				int origin = 33 
				) ;

		virtual Header fetch_header() ;

		virtual state get_state() { return str_ ? str_->get_state() : invalid ; }
		virtual Result fetch_result() { if( get_state() == have_output ) return str_->fetch_result() ; throw "calling sequence violated" ; }
		virtual Footer fetch_footer() { return str_->fetch_footer() ; }
#if HAVE_ELK_SCHEME_H
		virtual Object get_summary() const { return str_->get_summary() ; }
#endif
} ;


//! \brief stream that writes result in native (ANFO) format
//! The file will be in a format that can be read in by streams::AnfoReader.
class AnfoWriter : public Stream
{
	private:
		std::auto_ptr< google::protobuf::io::ZeroCopyOutputStream > zos_ ;
		google::protobuf::io::CodedOutputStream o_ ;
		Chan chan_ ;
		std::string name_ ;
		int64_t wrote_ ;

	public:
		AnfoWriter( google::protobuf::io::ZeroCopyOutputStream*, const char* = "<pipe>" ) ;
		AnfoWriter( int fd, const char* = "<pipe>", bool expensive = false ) ;
		AnfoWriter( const char* fname, bool expensive = false ) ;

		virtual void put_header( const Header& h ) { write_delimited_message( o_, 1, h ) ; Stream::put_header( h ) ; } 
		virtual void put_footer( const Footer& f ) { write_delimited_message( o_, 3, f ) ; Stream::put_footer( f ) ; }
		virtual void put_result( const Result& r ) ;
} ;

//! \brief new blocked native format
//! Writes in a format that can be read by stream::ChunkedReader.  The
//! file is made up of individually compressed blocks so that
//! near-random access is possible... in principle.
class ChunkedWriter : public Stream
{
	public:
		enum method { none, fastlz, gzip, bzip } ;

	private:
		std::auto_ptr< google::protobuf::io::ZeroCopyOutputStream > zos_ ;	// final output
		std::vector< char > buf_ ;											// in-memory buffer
		std::auto_ptr< google::protobuf::io::ArrayOutputStream > aos_ ;		// output to buffer
		Chan chan_ ;
		std::string name_ ;
		int64_t wrote_ ;
		uint8_t method_, level_ ;

		virtual ~ChunkedWriter() ;
		void flush_buffer( unsigned needed = 0 ) ;
		void init() ;

	public:
		static uint8_t method_of( int l ) {
#if HAVE_LIBBZ2 && HAVE_BZLIB_H
			if( l >= 75 ) return bzip ;
#endif
#if HAVE_LIBZ && HAVE_ZLIB_H
			if( l >= 50 ) return gzip ;
#endif
			if( l >= 10 ) return fastlz ;
			return none ;
		}
		static uint8_t level_of( int l ) {
			if( l >= 65 ) return 9 ;	// thorough gzip or bzip
			if( l >= 30 ) return 2 ;	// thorough fastlz or fast gzip
			return 1 ;					// fast fastlz or none 
		}

		ChunkedWriter( const pair< google::protobuf::io::ZeroCopyOutputStream*, string >&, int ) ;
		ChunkedWriter( int fd, int, const char* = "<pipe>" ) ;
		ChunkedWriter( const char* fname, int ) ;

		virtual void put_header( const Header& h ) ;
		virtual void put_result( const Result& r ) ;
		virtual void put_footer( const Footer& f ) ;
} ;

class ChunkedReader : public Stream
{
	private:
		std::auto_ptr< google::protobuf::io::ZeroCopyInputStream > is_ ;
		std::vector< char > buf_ ;											// in-memory buffer
		std::auto_ptr< google::protobuf::io::ArrayInputStream > ais_ ;		// output to buffer
		std::string name_ ;

		virtual ~ChunkedReader() { --anfo_reader__num_files_ ; }
		bool get_next_chunk() ;

	public: 
		ChunkedReader( std::auto_ptr< google::protobuf::io::ZeroCopyInputStream > is, const std::string& name ) ;
		virtual Result fetch_result() ;
		virtual Footer fetch_footer() ;

		//! \internal
		static unsigned num_open_files() { return anfo_reader__num_files_ ; }
} ;

//! \brief filters that drop or modify isolated records
//! Think "mapMaybe".
//! \todo Maybe some stats could be gathered into some sort of a result.

class Filter : public Stream
{
	public:
		virtual bool xform( Result& ) = 0 ;
		virtual void put_result( const Result& res ) { res_ = res ; if( xform( res_ ) ) state_ = have_output ; }
} ;

//! \brief filters that drop some alignments
//! Filtering only applies to hits to the specified genome(s), or to all
//! hits if no genomes are specified.  Other hits pass through.
//! \todo Maybe some stats could be gathered into some sort of a result.
class HitFilter : public Stream
{
	private:
		vector<string> gs_ ;

	public:
		HitFilter() {}
		HitFilter( const vector<string> &gs ) : gs_(gs) {}

		virtual void put_header( const Header& h ) {
			Stream::put_header( h ) ;
			hdr_.clear_is_sorted_by_coordinate() ;
			hdr_.clear_is_sorted_by_all_genomes() ;
		}

		virtual bool keep( const Hit& ) = 0 ;
		virtual void put_result( const Result& res ) ; 
} ;

namespace {
	template < typename C, typename V > bool contains( const C& c, const V& v )
	{ return find( c.begin(), c.end(), v ) != c.end() ; }
}

//! \brief deletes some alignments
//! A list of sequences and a list of genomes can be given.  An
//! alignment is deleted if a) the list of genomes is empty or the hit
//! genome is a member of that list and b) the list of sequences is
//! empty or the hit sequence is a member of that list.
class IgnoreHit : public HitFilter
{
	private:
		vector< string > ss_ ;

	public:
		IgnoreHit( const vector< string > &gs, const vector< string > &ss ) : HitFilter( gs ), ss_( ss ) {}
		virtual bool keep( const Hit& h ) { return !ss_.empty() && !contains( ss_, h.sequence() ) ; }
} ;

//! \brief deletes hits to uninteresting genomes
//! Hits to the specified genomes pass through, all others are dropped.
class OnlyGenome : public Filter
{
	private:
		vector< string > gs_ ;

	public:
		OnlyGenome( const vector< string > &gs ) : Filter(), gs_( gs ) {}
		virtual bool xform( Result& r ) ;
} ;

//! \brief stream that filters for a given score
//! Alignments that exceed the score are deleted, but the sequences are
//! kept.  If genomes are given, only alignments to that genome are
//! tested and removed, else all are considered individually and
//! selectively removed.  If no alignments remains, \c reason is set to
//! \c bad_alignment.  Most interesting downstream operations will
//! ignore sequences without alignments.  Scores are in Phred scale now,
//! experience tells that a typical butoff between good alignments and
//! junk is \f$ 7.5 \cdot ( \mbox{length} - 20 ) \f$

class ScoreFilter : public HitFilter
{
	private:
		double slope_ ;
		double intercept_ ;

	public:
		ScoreFilter( double s, double i, const vector<string> &gs ) : HitFilter(gs), slope_(s), intercept_(i) {}
		virtual bool keep( const Hit& ) ;
} ;

//! \brief stream that filters for minimum mapping quality
//! All alignments of sequences where the difference to the next hit is
//! too low are deleted.  Filtering can be restricted to some genomes.
//! The sequence itself is always kept.
class MapqFilter : public HitFilter
{
	private:
		int minmapq_ ; ;

	public:
		MapqFilter( const vector<string> &gs, int q ) : HitFilter(gs), minmapq_(q) {}
		virtual bool keep( const Hit& ) ;
} ;

//! \brief stream that filters for minimum sequence length
//! All alignments of sequences that are too short are deleted, relying
//! on down stream filters to completely get rid of the sequences
//! themselves.
class LengthFilter : public Filter
{
	private:
		int minlength_ ;

	public:
		LengthFilter( int l ) : minlength_(l) {}
		virtual bool xform( Result& ) ;
} ;

//! \brief a stream that removes sequences without a specific hit
//! This is intended to shrink a file by removing junk that didn't
//! align at all.
class RequireHit : public Filter
{
	private:
		vector< string > gs_, ss_ ;

	public:
		RequireHit( const vector<string> &gs, const vector<string> &ss ) : gs_(gs), ss_(ss) {}
		virtual bool xform( Result& ) ;
} ;

//! \brief a stream that removes sequences without a specific best hit
//! This is intended to shrink a file by removing junk that didn't
//! align to the expected genome.
class RequireBestHit : public Filter
{
	private:
		vector< string > gs_, ss_ ;

	public:
		RequireBestHit( const vector<string> &gs, const vector<string> &ss ) : gs_(gs), ss_(ss) {}
		virtual bool xform( Result& ) ;
} ;

class Subsample : public Filter
{
	private:
		const float f_ ;

	public:
		Subsample( float f ) : f_(f) {}
		virtual bool xform( Result& ) ;
} ;
//! \brief filters for minimum multiplicity
//! Only results that stem from duplicate removal with a minimum
//! multiplicity are retained.  Intended for the analysis of libraries
//! sequenced with high redundancy to lower sequencing error.
class MultiFilter : public Filter
{
	private:
		int n_ ;

	public:
		MultiFilter( int n ) : n_(n) {}
		virtual bool xform( Result& r ) { return r.member_size() >= n_ ; }
} ;

//! \brief stream that filters out low quality bases
//! Bases with insufficient quality are replaced by N symbols.
//! (Originally I used gap symbols, which doesn't make sense and
//! actually confused the legacy tool downstream.  Ns should be fine and
//! are actually the correct symbol in a certain sense.)
//! suppress the counting of low quality bases.
class QualFilter : public Filter
{
	private:
		int q_ ;

	public:
		QualFilter( int q ) : q_(q) {}
		virtual bool xform( Result& ) ;
} ;

/*! \brief stream with PCR duplicates removed

    A set of duplicates is (more or less by definition) a set of reads
    that map to the same position.
   
	Any set of duplicates is merged (retaining the original reads in an
	auxilliary structure), and a consensus is called with new quality
	scores.  A quality score (for base A) is defined as \f[ Q := -10
	\log_{10} P( \bar{A} | \omega ), \f] that is, a representation of
	the probability of the base being wrong given some observation \f$
	\omega. \f$  It follows that \f[ Q = -10 \log_{10} \frac{
	P(\omega|\bar{A}) P(\bar{A}) }{ P(\omega) } = -10 \log_{10} P( \omega
	| \bar{A} ) \f] by assuming a uniform base composition in the sample
	(\f$ P(A)=\frac{1}{4} \f$) and after base calling (\f$
	P(\omega)=\frac{1}{4} \f$).

	When combining observations \f$ \omega_1, \ldots, \omega_n \f$
	conditional on the same base, they are all independent, hence the
	probabilities conditioned on the actual base can be multiplied.  To
	get the final quality when calling an A:
	\f[ P(A|\bigcap_n \omega_n) = \frac{ P(\bigcap_n \omega_n | A) P(A) }{ P(\bigcap_n \omega_n) }
	                            = \frac{ P(A) \prod_n P(\omega_n | A) }{ P(\bigcap_n \omega_n) } \f]

	The unconditional probabilities in the denominator are not
	independent (for they collectively depend on the true base), we have
	to replace them by total probabilites:
	\f[ = \frac{ P(A) \prod_n P(\omega_n | A) }{ \sum_N P(\bigcap_n \omega_n | N) P(N) } \f]

	Again, assuming a uniform base composition, the priors are all equal
	and cancel, and now the conditional probabilities are again
	independent.
	\f[ = \frac{ \prod_n P(\omega_n | A) }{ \sum_N \prod_n P(\omega_n | N) } \f]

	Tracking this incrementally is easy, we need to store the four
	products.  Computation needs to be done in the log domain to avoid
	loss of precision.  Likewise, computation of a quality must not involve \f$ P(A|\omega) \f$,
	since that is too close to one to be representable, so we calculate the quality as
	\f[ Q := -10 \log_{10} P( \bar{A} | \omega ) = -10 \log_{10} ( P(C | \omega ) + P(G | \omega ) + P(T | \omega ) ) \f]

	What remains is how to get an estimate of \f$
	P(\omega | N) \f$.  We easily get \f$ P(\omega | A) \f$ (see above),
	which will normally be very close to one, and the sum of the other
	three.  To estimate the summands, we set
	\f[ P(\omega|N) = \frac{P(N|\omega) * P(\omega)}{P(N)} =
	P(N|\omega) = P(N|\bar{A}) * P(\bar{A}|\omega) \f]
	and estimate \f$ P(N|\bar{A}) \f$ from misclassfication statistics
	of the base caller.
 */

class RmdupStream : public Stream
{
	private:
		output::Result cur_ ;
		std::vector< Logdom > quals_[4] ;
		double slope_ ;
		double intercept_ ;
		vector<string> gs_ ;
		int maxq_ ;

		// XXX double err_prob_[4][4] ; // get this from config or
		// something?

		bool is_duplicate( const Result& , const Result& ) ;
		void add_read( const Result& ) ;
		void call_consensus() ;

	public:
		//! \brief sets parameters
		//! \param s Slope of score function, bad alignments are disregarded.
		//! \param i Intercept of score function.
		//! \param q (Assumed) quality of the polymerase.
		RmdupStream( double s, double i, int q ) : slope_(s), intercept_(i), maxq_( std::min(q,127) ) {}

		virtual void put_header( const Header& h )
		{
			if( !h.has_is_sorted_by_all_genomes() && !h.is_sorted_by_coordinate_size() )
				throw "RmdupStream: need sorted stream to remove duplicates" ;
			gs_.assign( h.is_sorted_by_coordinate().begin(), h.is_sorted_by_coordinate().end() ) ;
			Stream::put_header( h ) ;
		}

		virtual void put_result( const Result& ) ;
		virtual void put_footer( const Footer& ) ;
		virtual Result fetch_result() ;
} ;

//! \brief a stream that concatenates its input streams
//! The headers and footers are merged sensibly (plain merge with removal of
//! redundant information), then the results are simply concatenated.
class ConcatStream : public StreamBundle
{
	public:
		virtual Header fetch_header() ;
		virtual Result fetch_result() ;
} ;

class FastqReader : public Stream
{
	private:
		std::auto_ptr< google::protobuf::io::ZeroCopyInputStream > is_ ;
		bool sol_scores_ ;
		char origin_ ;

		void read_next_message() {
			state_ = read_fastq( is_.get(), *res_.mutable_read(), sol_scores_, origin_ )
				? have_output : end_of_stream ;
			sanitize( *res_.mutable_read() ) ;
		}

	public: 
		FastqReader( std::auto_ptr< google::protobuf::io::ZeroCopyInputStream > is, bool solexa_scores, char origin ) ;
		virtual Result fetch_result() { Result r ; std::swap( r, res_ ) ; read_next_message() ; return r ; }
} ;

class SffReader : public Stream
{
	private:
		std::auto_ptr< google::protobuf::io::ZeroCopyInputStream > is_ ;
		string name_ ;
		unsigned remaining_ ;
		unsigned number_of_flows_ ;

		const void* buf_ ;
		int buf_size_ ;

		uint8_t read_uint8() ;
		uint16_t read_uint16() ;
		uint32_t read_uint32() ;
		void read_string( unsigned, string* ) ;
		void skip( int ) ;

	public:
		SffReader( auto_ptr< google::protobuf::io::ZeroCopyInputStream > is, const string& name ) : 
			is_( is ), name_( name ), buf_(0), buf_size_(0) {}

		virtual Header fetch_header() ;
		virtual Result fetch_result() ;
} ;

} // namespace streams

class PipeOutputStream : public google::protobuf::io::FileOutputStream 
{
	private:
		pid_t cpid_ ;

	public:
		PipeOutputStream( int fd, pid_t cpid ) : google::protobuf::io::FileOutputStream( fd ), cpid_( cpid ) {} 
		virtual ~PipeOutputStream() { Close() ; throw_errno_if_minus1( waitpid( cpid_, 0, 0 ), "waiting for pipe process" ) ; }
} ;

std::pair< PipeOutputStream*, std::string > make_PipeOutputStream( const std::string& ) ;

//! \brief adapts a ZeroCopyOutputStream into a streambuf
class zero_copy_output_buf : public std::streambuf {
	private:
		std::auto_ptr< google::protobuf::io::ZeroCopyOutputStream > os_ ;

	public:
		zero_copy_output_buf( google::protobuf::io::ZeroCopyOutputStream* os ) : os_(os) {}
		virtual ~zero_copy_output_buf() ;
		virtual int sync() ;
		virtual int_type overflow( int_type ) ;
} ;

class zero_copy_ostream : public std::ostream {
	private:
		zero_copy_output_buf b_ ;

	public:
		zero_copy_ostream( google::protobuf::io::ZeroCopyOutputStream* os ) 
			: std::ostream( &b_ ), b_( os ) {}
} ;

#endif
