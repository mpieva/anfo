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

namespace streams {
	using namespace output ;

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

void sanitize( Header& hdr ) ;
void merge_sensibly( output::Header& lhs, const output::Header& rhs ) ;
void merge_sensibly( output::Footer& lhs, const output::Footer& rhs ) ;
void merge_sensibly( output::Result& lhs, const output::Result& rhs ) ;

//! \brief checks if a genome was hit
//! If an empty genome is asked for, checks for any hit.
bool has_hit_to( const output::Result&, const char* ) ;

//! \brief returns the hit to some genome
//! If an empty genome is asked for, returns the best hit.  Behaviour is
//! undefined if no suitable hit exists.
const output::Hit& hit_to( const output::Result&, const char* ) ;

//! \brief returns the mutable hit to some genome
//! If an empty genome is asked for, returns the best hit.  If no
//! suitable hit exists, a new one is created.
output::Hit* mutable_hit_to( output::Result*, const char* ) ;

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
	if( !s.empty() && (streams::cigar_op( s.back() ) == op) && streams::cigar_len( s.back() ) ) 
		s.back() = streams::mk_cigar( op, m + streams::cigar_len( s.back() ) ) ;
	else s.push_back( streams::mk_cigar( op, m ) ) ;
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
		enum state { invalid, end_of_stream, need_input, have_output } ;

	protected:
		// internal state---not strictly necessary here, but used almost
		// everywhere and therefore simply convenient
		Header hdr_ ;
		Result res_ ;
		Footer foot_ ;
		state state_ ;

	public:
		Stream() : state_( end_of_stream ) {}
		virtual ~Stream() {}

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
		state get_state() const { return state_ ; }

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
} ;

void transfer( Stream& in, Stream& out ) ;

//! \brief base class of streams that read from many streams
class StreamBundle : public Stream
{
	protected:
		std::deque< Stream* > streams_ ;
		typedef std::deque< Stream* >::const_iterator citer ;
		typedef std::deque< Stream* >::const_reverse_iterator criter ;

	public:
		virtual ~StreamBundle() { for_each( streams_.begin(), streams_.end(), delete_ptr<Stream>() ) ; }
		
		//! \brief adds an input stream
		//! The stream is taken ownership of and freed when the StreamBundle
		//! is destroyed.
		virtual void add_stream( Stream* ) = 0 ;
} ;

//! \brief presents ANFO files as series of messages
//! This class will read a possibly compressed result stream and present
//! it using an iteration interface.  Normally, gzip and bzip2
//! decompression will transparently be done.
class AnfoReader : public Stream
{
	private:
		std::auto_ptr< google::protobuf::io::ZeroCopyInputStream > is_ ;
		std::string name_ ;

		void read_next_message( google::protobuf::io::CodedInputStream& ) ;

		//! \internal
		//! Tracked to avoid bumping into the file descriptor limit
		//! (mostly important for merge sorting and mega-merge).
		static int num_files_ ;

	public: 
		struct ParseError : public Exception {
			std::string msg_ ;
			ParseError( const std::string& m ) : msg_(m) {}
			virtual void print_to( std::ostream& ) const ;
		} ;

		AnfoReader( google::protobuf::io::ZeroCopyInputStream *is, const std::string& name ) ;

		virtual ~AnfoReader() { --num_files_ ; }
		virtual Result fetch_result() ;

		//! \internal
		static unsigned num_open_files() { return num_files_ ; }
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

//! \brief filters that drop or modify isolated records
//! Think "mapMaybe".
//! \todo Maybe the fact that some filtering was done should be recorded
//! in a header, maybe some stats could be put into the footer.

class Filter : public Stream
{
	public:
		virtual bool xform( Result& ) = 0 ;
		virtual void put_result( const Result& res ) { res_ = res ; if( xform( res_ ) ) state_ = have_output ; }
} ;

//! \brief stream that filters for a given score
//! Alignments that exceed the score are deleted, but the sequences are
//! kept.  If a genome is given, only alignments to that genome are
//! removed, else all are considered individually and selectively
//! removed.  If no alignments remains, \c reason is set to \c
//! bad_alignment.  Most interesting downstream operations will ignore
//! sequences without alignments.  Scores are in Phred scale now,
//! experience tells that a typical butoff between good alignments and
//! junk is \f$ 7.5 \cdot ( \mbox{length} - 20 ) \f$

class ScoreFilter : public Filter
{
	private:
		double slope_ ;
		double intercept_ ;
		const char *genome_ ;

	public:
		ScoreFilter( double s, double i, const char* g ) : slope_(s), intercept_(i), genome_(g) {}
		virtual ~ScoreFilter() {}
		virtual bool xform( Result& ) ;
} ;

//! \brief stream that filters for minimum mapping quality
//! All alignments of sequences where the differenc eto the next hit is
//! too low are deleted.
class MapqFilter : public Filter
{
	private:
		const char* g_ ;
		int minmapq_ ; ;

	public:
		MapqFilter( const char* g, int q ) : g_(g), minmapq_(q) {}
		virtual ~MapqFilter() {}
		virtual bool xform( Result& ) ;
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
		virtual ~LengthFilter() {}
		virtual bool xform( Result& ) ;
} ;

//! \brief a stream that removes sequences without a hit
//! This is intended to shrink a file by removing junk that didn't
//! align.  It's not needed for conversion to SAM or similar, since
//! those filters drop unaligned sequences internally.
class HitFilter : public Filter
{
	private:
		const char* g_ ;
		const char* s_ ;

	public:
		HitFilter( const char* g, const char* s ) : g_(g), s_(s) {}
		virtual ~HitFilter() {}
		virtual bool xform( Result& ) ;
} ;

class Subsample : public Filter
{
	private:
		const float f_ ;

	public:
		Subsample( float f ) : f_(f) {}
		virtual ~Subsample() {}
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
		virtual ~MultiFilter() {}
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
		virtual ~QualFilter() {}
		virtual bool xform( Result& ) ;
} ;

/*! \brief stream with PCR duplicates removed
    A set of duplicate is (more or less by definition) a set of reads
    that map to the same position.
   
    \todo Refine the definition of 'duplicate'.  It appears sensible to
		  forbid alignments that hang off of a contig, though maybe the
		  requirement for the same length will catch most of them.
	\todo Check what happens when sequences are trimmed (compare effective lengths?)
   
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

	When combining observations \f$ \omega_1, \ldots, \omega_n \f$, they
	are all independent, hence the probabilities conditioned on the
	actual base can be multiplied.  To get the final quality when
	calling an A:
	\f[ P(A|\bigcap_n \omega_n) = \frac{ P(\bigcap_n \omega_n | A) P(A) }{ P(\bigcap_n \omega_n) }
	                            = \frac{ P(A) \prod_n P(\omega_n | A) }{ P(\bigcap_n \omega_n) } \f]

	The unconditional probabilities in the denominator are not
	independent, we have to replace them by total probabilites:
	\f[ = \frac{ P(A) \prod_n P(\omega_n | A) }{ \sum_N P(\bigcap_n \omega_n | N) P(N) } \f]

	Again, assuming a uniform base composition, the priors are all equal
	and cancel, and now the conditional probabilities are again
	independent.
	\f[ = \frac{ \prod_n P(\omega_n | A) }{ \sum_N \prod_n P(\omega_n | N) } \f]

	Tracking this incrementally is easy, we need to store the four
	products.  Computation needs to be done in the log domain to avoid
	loss of precision.  What remains is how to get an estimate of \f$
	P(\omega | N) \f$.  We easily get \f$ P(\omega | A) \f$ (see above),
	which will normally be very close to one, and the sum of the other
	three.  To distribute the sum, we set
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
		const char *g_ ;
		int maxq_ ;

		// XXX double err_prob_[4][4] ; // get this from config or
		// something?

		bool is_duplicate( const Result& , const Result& ) ;
		void add_read( const Result& ) ;
		void call_consensus() ;

	public:
		RmdupStream( double s, double i, int maxq ) : slope_(s), intercept_(i), g_(0), maxq_( std::min(maxq,127) ) {}
		virtual ~RmdupStream() { free( const_cast<char*>(g_) ) ; }

		virtual void put_header( const Header& h )
		{
			if( !h.has_is_sorted_by_coordinate() ) throw "RmdupStream: need sorted stream to remove duplicates" ;
			g_ = h.is_sorted_by_coordinate().empty() ? 0 : strdup( h.is_sorted_by_coordinate().c_str() ) ;
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
	private:
		std::deque< output::Result > rs_ ;

	public:
		virtual void add_stream( Stream* ) ;
		virtual Result fetch_result() ;
} ;

Stream* make_input_stream( const char* name, bool solexa_scores = false, char origin = 33 ) ;
Stream* make_input_stream( int fd, const char* name = "<pipe>", bool solexa_scores = false, char origin = 33 ) ;
Stream* make_input_stream( google::protobuf::io::ZeroCopyInputStream *is, const char* name = "<pipe>", int64_t total = -1, bool solexa_scores = false, char origin = 33 ) ;


class FastqReader : public Stream
{
	private:
		std::auto_ptr< google::protobuf::io::ZeroCopyInputStream > is_ ;
		bool sol_scores_ ;
		char origin_ ;

		void read_next_message() {
			state_ = read_fastq( is_.get(), *res_.mutable_read(), sol_scores_, origin_ )
				? have_output : end_of_stream ;
		}

	public: 
		FastqReader( google::protobuf::io::ZeroCopyInputStream *is, bool solexa_scores, char origin ) ;
		virtual Result fetch_result() { Result r ; std::swap( r, res_ ) ; read_next_message() ; return r ; }
} ;

} // namespace streams

#endif
