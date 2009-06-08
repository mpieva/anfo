#ifndef INCLUDED_STREAM_H
#define INCLUDED_STREAM_H

#include "compress_stream.h"
#include "logdom.h"
#include "output.pb.h"
#include "util.h"

#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>

#include <deque>
#include <fstream>
#include <memory>

namespace streams {
	using namespace output ;

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
	if(
			!os.WriteTag( (tag << 3) | 2 ) ||
			!os.WriteVarint32( m.ByteSize() ) ||
			!m.SerializeToCodedStream( &os )
	  )
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

void merge_sensibly( output::Header& lhs, const output::Header& rhs ) ;
void merge_sensibly( output::Footer& lhs, const output::Footer& rhs ) ;
void merge_sensibly( output::Result& lhs, const output::Result& rhs ) ;

//! \brief computes (trimmed) query length from CIGAR line
unsigned len_from_bin_cigar( const std::string& cigar ) ;

//! @}

struct MissingMethod : public Exception {
	std::string n_ ;
	MissingMethod( const std::string& n ) : n_(n) {}
	virtual void print_to( std::ostream& s ) const { s << "missing method: " << n_ ; }
} ;

//! \brief stream of result messages
//! Each stream is a header, followed by many results, followed by a
//! single footer.  The header will be cached internally, so it can be
//! asked for repeatedly.  Results are forgotten once read, the footer
//! is only available after all results have been read, and then it is
//! stored and can be read repeatedly.  It is undefined behaviour to
//! request the footer without first consuming all results.
//!
//! XXX
//!
//! Trying a unified design: this class is both an input stream, an
//! output stream and a stream transducer.  Concrete implementations
//! decide what sets of methods to actually implement.  

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
		Stream() : state_( invalid ) {}
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
		virtual Result fetch_result() { throw MissingMethod(__FUNCTION__) ; }

		//! \brief returns the footer
		//! Only after all results have been consumed is the footer
		//! available.  If anything goes wrong internally, the LSB
		//! should be set in \c exit_code.
		virtual Footer fetch_footer() { return foot_ ; }

		//! \brief sets the stream header
		//! Output streams and stream filters need the header to become
		//! valid streams.
		virtual void put_header( const Header& ) { throw MissingMethod(__FUNCTION__) ; }

		//! \brief outputs one result
		//! This method can only be called in state need_input, and it
		//! will insert a result record into the stream.  Else the
		//! behaviour is undefined.
		virtual void put_result( const Result& ) { throw MissingMethod(__FUNCTION__) ; }

		//! \brief sets the footer
		//! This method can only be called in state need_input, and it
		//! signals that no more input is available.  A filter will then
		//! flush internal buffers and signal end_of_stream.
		virtual void put_footer( const Footer& ) { throw MissingMethod(__FUNCTION__) ; }
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
		google::protobuf::io::FileInputStream iis_ ;
		std::auto_ptr<google::protobuf::io::ZeroCopyInputStream> zis_ ;
		std::string name_ ;
		bool quiet_ ;

		void initialize() ;
		void read_next_message( google::protobuf::io::CodedInputStream& ) ;

		//! \internal
		//! Tracked to avoid bumping into the file descriptor limit
		//! (mostly important for merge sorting and mega-merge).
		static int num_files_ ;

	public: 
		//! \brief opens the named file
		AnfoReader( const std::string& name, bool quiet = false ) ;

		//! \brief uses a given name and filedescriptor
		//! No actual file is touched, the name is for informational
		//! purposes only.
		AnfoReader( int fd, const std::string& name = "<unknown>", bool quiet = false ) ;

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
		std::auto_ptr< google::protobuf::io::FileOutputStream > fos_ ;
		std::auto_ptr< google::protobuf::io::ZeroCopyOutputStream > zos_ ;
		google::protobuf::io::CodedOutputStream o_ ;

	public:
		AnfoWriter( google::protobuf::io::ZeroCopyOutputStream* ) ;
		AnfoWriter( int fd, bool expensive = false ) ;
		AnfoWriter( const char* fname, bool expensive = false ) ;


		virtual void put_header( const Header& h ) { write_delimited_message( o_, 1, h ) ; state_ = need_input ; } 
		virtual void put_result( const Result& r ) { write_delimited_message( o_, 2, r ) ; }
		virtual void put_footer( const Footer& f ) { write_delimited_message( o_, 3, f ) ; state_ = end_of_stream ; }
} ;

//! \brief filters that drop or modify isolated records
//! Think "mapMaybe".

class Filter : public Stream
{
	public:
		virtual bool xform( Result& ) = 0 ;

		virtual void put_header( const Header& hdr ) { hdr_ = hdr ; state_ = need_input ; }// XXX leave some traces? 
		virtual void put_footer( const Footer& foot ) { foot_ = foot ; state_ = end_of_stream ; }// XXX leave some traces? 

		virtual void put_result( const Result& res ) { res_ = res ; if( xform( res_ ) ) state_ = have_output ; }
		virtual Result fetch_result() { state_ = need_input ; return res_ ; }
} ;

//! \brief stream that filters for a given score
//! Alignments that exceed the score are deleted, but the sequences are
//! kept.  (Most interesting operations will ignore sequences without
//! alignments).
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

	public:
		HitFilter( const char* g ) : g_(g) {}
		virtual ~HitFilter() {}
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
//! Bases with insufficient quality are replaced by gap symbols.
//! Strictly speaking gap symbols in a sequence don't make sense and
//! would confuse many tools, but the intention here is to apply this
//! filter only when producing legacy FASTA out where it serves to
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
	\todo Include a score cutoff so only good alignments are taken into account.
	\todo Supply the genome the coordinates of which we want to look at.
   
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
	\f[ P(\omega|N) = \frac{P(N|\omega) * P(\omega)}{\P(N)} =
	P(N|\omega) = P(N|\bar{A}) * P(\bar{A}|\omega) \f]
	and estimate \f$ P(N|\bar{A}) \f$ from misclassfication statistics
	of the base caller.
 */

// XXX control flow is messed up
#if 0
class RmdupStream : public Stream
{
	private:
		output::Result cur_ ;
		std::vector< Logdom > quals_[4] ;
		bool good_ ;
		// XXX double err_prob_[4][4] ;

	public:
		RmdupStream( Input *str ) : str_( str )
		{
			assert( str_->get_header().is_sorted_by_coordinate() ) ;
			good_ = str_->read_result( cur_ ) ;
		}

		virtual const output::Header& get_header() { return str_->get_header() ; }
		virtual const output::Footer& get_footer() { return str_->get_footer() ; }
		virtual bool read_result( output::Result& ) ;

	private:
		static bool is_duplicate( const output::Result& lhs, const output::Result& rhs ) ;
		void add_read( const output::Result& rhs ) ;
} ;
#endif

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

} // namespace streams

#endif
