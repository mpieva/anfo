#include "config.h"

#include "compress_stream.h"
#include "logdom.h"
#include "outputfile.h"
#include "util.h"

#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>

#include <algorithm>
#include <cstdlib>
#include <deque>
#include <iostream>

#include <glob.h>

#if HAVE_ALLOCA_H
#  include <alloca.h>
#endif

#if HAVE_UNISTD_H
#  include <unistd.h>
#endif

using namespace google::protobuf::io ;
using namespace std ;
using namespace output ;

static unsigned max_que_size = 256 ;
static unsigned max_arr_size = 1024*1024*1024 ;

//! \page anfo-sort
//! \brief reads ANFO files, sorts and merges them.
//!
//! Outline: We convert a list of files (use glob matching here, no need
//! to rely on the shell) into a single file.  First we need to collect
//! files that operated on the same sequences, but different indices and
//! keep the best hit from each of them.  Then we keep a queue of sorted
//! streams, whereby a stream can be a sorted temporary file or a
//! virtual stream that sorts on the fly.  Streams are ordered according
//! to estimated size, whenever we run the danger of running out of file
//! descriptors, we merge small streams into a temporary file containing
//! a bigger one.  In the end, everything is merged to output.
//!
//! Are parameters necessary?  Certainly an output and many input files;
//! the number of different indices is encoded in the header.  We run
//! into problem if multiple distinct indices were used, it might be
//! useful to have an override for that case.  A directory for temporary
//! files is useful, but we take that from the environment.  As far as
//! genome files are needed, we can take the configuration from the
//! input files.
//!
//! \todo implement copnsensus calling of PCR duplicates
//! \todo implement score cutoff (we don't want crap alignments, and we
//!       certainly don't want to include crap sequences in the
//!       consensus calling)


//! \brief compares hits by smallest genome coordinate
//! Comparison is first done lexically on the subject name, then on the
//! smallest coordinate that's part of the alignment, then on the
//! alignment length.
struct by_genome_coordinate {
	bool operator() ( const Result *a, const Result *b ) {
		if( a->has_best_to_genome() && !b->has_best_to_genome() ) return true ;
		if( !a->has_best_to_genome() ) return false ;
		const Hit& u = a->best_to_genome(), v = b->best_to_genome() ;
		if( u.genome_name() < v.genome_name() ) return true ;
		if( v.genome_name() < u.genome_name() ) return false ;
		if( u.sequence() < v.sequence() ) return true ;
		if( v.sequence() < u.sequence() ) return false ;
		if( u.start_pos() < v.start_pos() ) return true ;
		if( v.start_pos() < u.start_pos() ) return false ;
		return u.aln_length() < v.aln_length() ;
	}
} ;


template< typename T > struct delete_ptr { void operator()( T* p ) const { delete p ; } } ;

class MergeStream : public Stream
{
	private:
		vector< Stream* > streams_ ;
		deque< Result > rs_ ;
		output::Header hdr_ ;
		output::Footer foot_ ;

	public:
		MergeStream() {}
		virtual ~MergeStream() { for_each( streams_.begin(), streams_.end(), delete_ptr<Stream>() ) ; }

		void add_stream( auto_ptr< Stream > s )
		{
			assert( s->get_header().is_sorted_by_coordinate() ) ;
			merge_sensibly( hdr_, s->get_header() ) ;

			Result r ;
			if( s->read_result( r ) ) 
			{
				rs_.push_back( r ) ;
				streams_.push_back( s.release() ) ;
			}
			else
			{
				merge_sensibly( foot_, s->get_footer() ) ;
			}
		}

		virtual const output::Header& get_header() { return hdr_ ; }
		virtual const output::Footer& get_footer() { return foot_ ; }
		virtual bool read_result( output::Result& ) ;
} ;

bool MergeStream::read_result( output::Result& res ) 
{
	if( rs_.empty() ) return false ;
	int min_idx = 0 ;
	for( size_t i = 1 ; i != rs_.size() ; ++i ) 
		if( by_genome_coordinate()( &rs_[ i ], &rs_[ min_idx ] ) )
			min_idx = i ;

	res = rs_[ min_idx ] ;
	Stream *sm = streams_[ min_idx ] ;
	if( !sm->read_result( rs_[ min_idx ] ) ) 
	{
		merge_sensibly( foot_, sm->get_footer() ) ;
		delete sm ;
		streams_.erase( streams_.begin() + min_idx ) ;
		rs_.erase( rs_.begin() + min_idx ) ;
	}
	return true ;
}


//! \brief merges multiple streams by taking the best hit
//! If asked for a result, this class takes results from each of the
//! streams in turn, until one sequence has received an entry from each.
//! This sequence is then delivered.  Everything works fine, as long as
//! the sequence names come in in roughly the same order.  Else it still
//! works, but eats memory.
class BestHitStream : public Stream
{
	private:
		vector< Stream* > streams_ ;
		output::Header hdr_ ;
		output::Footer foot_ ;

		typedef map< string, pair< size_t, output::Result > > Buffer ;
		Buffer buffer_ ;
		size_t cur_input_, nread_, nwritten_ ;

	public:
		BestHitStream() : streams_(), cur_input_(0), nread_(0), nwritten_(0) {}
		virtual ~BestHitStream() {for_each(streams_.begin(),streams_.end(),delete_ptr<Stream>());}

		//! \brief reads a stream's header and adds the stream as input
		void add_stream( auto_ptr< Stream > s ) {
			merge_sensibly( hdr_, s->get_header() ) ;
			streams_.push_back( s.release() ) ;
		}

		//! \brief checks if enough inputs are known
		//! Checks if at least as many inputs are known as requested.
		//! If that number isn't known, it assumed only one index was
		//! used and takes the number of slices declared for it.  
		bool enough_inputs( size_t n ) const { return streams_.size() >= n ; }
		bool enough_inputs() const
		{
			if( streams_.empty() ) return false ;
			if(  hdr_.config().policy_size() == 0 ) return true ;
			if(  hdr_.config().policy(0).use_compact_index_size() == 0 ) return true ;
			if( !hdr_.config().policy(0).use_compact_index(0).has_number_of_slices() ) return true ;
			return enough_inputs( hdr_.config().policy(0).use_compact_index(0).number_of_slices() ) ;
		}

		//! \brief returns the merged headers
		//! Redundant information is removed from the headers (e.g.
		//! repeated paths), the rest is merged as best as possible.
		virtual const output::Header& get_header() { return hdr_ ; }
		
		//! \brief reads a footer
		//! All footers are merged, the exit codes are logically OR'ed,
		//! and the LSB is set if something wen't wrong internally.
		virtual const output::Footer& get_footer() { return foot_ ; }


		//! \brief reports on one sequence
		//! Enough information is read until one sequence has been
		//! described by each input stream.  Everything else is buffered
		//! if necessary.
		virtual bool read_result( output::Result& res )
		{
			for( bool cont = true ; cur_input_ || cont ; )
			{
				if( !cur_input_ ) cont = false ;
				Result r ;
				bool good = streams_[cur_input_]->read_result( r ) ;
				if( !good ) merge_sensibly( foot_, streams_[cur_input_]->get_footer() ) ;

				if( ++cur_input_ == streams_.size() ) cur_input_ = 0 ;
				if( good )
				{
					cont = true ;
					pair< size_t, Result > &p = buffer_[ r.seqid() ] ;
					++nread_ ;
					++p.first ;
					if( p.second.has_seqid() ) merge_sensibly( p.second, r ) ;
					else p.second = r ;

					if( p.first == streams_.size() ) 
					{
						++nwritten_ ;
						res = p.second ;
						buffer_.erase( r.seqid() ) ;

						if( ((nread_ + nwritten_) & 0xFFFF) == 0 ) 
							clog << "\033[KRead " << nread_ << ", delivered "
								<< nwritten_ << ", buffering "
								<< buffer_.size() << '\r' << flush ;

						return true ;
					}
				}
			}

			if( buffer_.empty() ) return false ;
			res = buffer_.begin()->second.second ;
			buffer_.erase( buffer_.begin()->first ) ;
			return true ;
		}
} ;

template< typename I > class StreamAdapter : public Stream
{
	private:
		Header hdr_ ;
		I cur_, end_ ;
		Footer foot_ ;
	public:
		StreamAdapter( const Header &hdr, I begin, I end, const Footer &foot ) 
			: hdr_( hdr ), cur_( begin ), end_( end ), foot_( foot ) {}

		virtual const Header& get_header() { return hdr_ ; }
		virtual const Footer& get_footer() { return foot_ ; }
		virtual bool read_result( Result& r ) {
			if( cur_ == end_ ) return false ;
			r = *( *cur_++ ) ;
			return true ;
		}
} ;


/*! \brief stream with PCR duplicates removed
    A set of duplicate is (more or less by definition) a set of reads
    that map to the same position.
   
    \todo Refine the definition of 'duplicate'.  It appears sensible to
          also require the sequences to have the same length
          (effectively always the case on Solexa, unless the first base
          is lost, which already kills the algorithm; it's not so easy
          on 454).  It also appears worthwhile to forbid alignments that
          hang off of a contig, though maybe the requirement for the
          same length will catch them.
   
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

class RmdupStream : public Stream
{
	private:
		Stream* str_ ;
		output::Result cur_ ;
		std::vector< Logdom > quals_[4] ;
		bool good_ ;

	public:
		RmdupStream( Stream *str ) : str_( str )
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

bool RmdupStream::is_duplicate( const output::Result& lhs, const output::Result& rhs ) 
{
	if( lhs.has_sequence() && rhs.has_sequence() && lhs.sequence() != rhs.sequence() ) return false ;
	if( !lhs.has_best_to_genome() || !rhs.has_best_to_genome() ) return false ;
	const output::Hit &l = lhs.best_to_genome(), &r = rhs.best_to_genome() ;

	bool z = l.genome_name() == r.genome_name() && l.sequence() == r.sequence()
		&& l.start_pos() == r.start_pos() && l.aln_length() == r.aln_length() ;
	return z ;
}

void RmdupStream::add_read( const output::Result& rhs ) 
{
	output::Read *rd = cur_.add_member() ;
	if( rhs.has_seqid() ) rd->set_seqid( rhs.seqid() ) ;
	if( rhs.has_description() ) rd->set_description( rhs.description() ) ;
	if( rhs.has_sequence() ) rd->set_sequence( rhs.sequence() ) ;
	if( rhs.has_quality() ) rd->set_quality( rhs.quality() ) ;

	// XXX How do we deal with ambiguity codes?  What's the meaning of
	// their quality scores?
	for( int i = 0 ; i != abs( rhs.best_to_genome().aln_length() ) ; ++i )
	{
		int base = -1 ;
		switch( rhs.sequence()[i] ) {
			case 'a': case 'A': base = 0 ; break ;
			case 'c': case 'C': base = 1 ; break ;
			case 't': case 'T':
			case 'u': case 'U': base = 2 ; break ;
			case 'g': case 'G': base = 3 ; break ;
		}

		Logdom qual = Logdom::from_phred( rhs.has_quality() ? rhs.quality()[i] : 30 ) ;
		for( int j = 0 ; j != 4 ; ++j )
			// XXX distribute errors sensibly
			quals_[j].at(i) *= j != base 
				? qual / Logdom::from_float( 3 ) 
				: Logdom::from_float(1) - qual ;
	}
}

bool RmdupStream::read_result( output::Result& r ) 
{
	if( !good_ ) return false ;

	// Get next result from input stream; it can either be merged or
	// not.  If merged, continue.  Else return whatever we accumulated
	// and store the one that didn't fit.

	output::Result next ;
	while( (good_ = str_->read_result( next )) && is_duplicate( cur_, next ) )
	{
		// if cur is a plain result, turn it into a degenerate merged one
		if( cur_.member_size() == 0 )
		{
			for( size_t i = 0 ; i != 4 ; ++i )
			{
				quals_[i].clear() ;
				quals_[i].resize( abs( cur_.best_to_genome().aln_length() ) ) ;
			}

			add_read( cur_ ) ;
			cur_.set_seqid( "C_" + cur_.seqid() ) ;
			cur_.clear_description() ;
		}
		add_read( next ) ;
	}

	// input stream ended or alignments don't match --> return accumulator
	// re-basecall iff we actually did merge
	if( cur_.member_size() ) 
	{
		cur_.clear_sequence() ;
		cur_.clear_quality() ;
		for( int i = 0 ; i != abs( cur_.best_to_genome().aln_length() ) ; ++i )
		{
			size_t m = 0 ;
			for( size_t j = 1 ; j != 4 ; ++j )
				// select base with highest quality
				if( quals_[j].at( i ) > quals_[m].at( i ) ) m = j ;
			// but calculate the error probability from _all_others_ to
			// retain precision

			for( size_t j = 0 ; j != 4 ; ++j )
				std::cerr << ' ' << (j==m ? '*' : ' ') << ' ' << quals_[j].at( i ).to_float() << std::endl ;

			Logdom num = quals_[m?0:1].at(i) ;
			Logdom denom = quals_[0].at(i) ;
			for( size_t j = 1 ; j != 4 ; ++j )
			{
				denom += quals_[j].at( i ) ;
				// std::cerr << "     " << denom.to_float() << std::endl ;
				if( j != m ) num += quals_[j].at( i ) ;
			}

			int qscore = (num/denom).to_phred() ;
			if( qscore > 127 ) qscore = 127 ;

			std::cerr << num.to_float() << '/' << denom.to_float() 
				<< " == " << (num/denom).to_phred() << std::endl ;
			cur_.mutable_sequence()->push_back( "ACTG"[m] ) ;
			cur_.mutable_quality()->push_back( qscore ) ;
		}
	}

	r = cur_ ;
	cur_ = next ;
	return true ;
}

// keep queues of streams ordered and separate, so we don't merge the
// same stuff repeatedly
typedef deque< Stream* > MergeableQueue ;
typedef map< unsigned, MergeableQueue > MergeableQueues ;
MergeableQueues mergeable_queues ;

typedef deque< Result* > ScratchSpace ;
ScratchSpace scratch_space ;
Header scratch_header ;
Footer scratch_footer ;
int64_t total_scratch_size = 0 ;

void sort_scratch() {
    if( scratch_space.size() > 1 )
		clog << "\033[KSorting " << scratch_space.size() << " results in memory" << endl ;
    sort( scratch_space.begin(), scratch_space.end(), by_genome_coordinate() ) ;
}

int mktempfile( string &name )
{
	const char *suffix = "/anfo_sort_XXXXXX" ;
	const char *base = getenv("ANFO_TEMP") ;
	if( !base ) base = getenv("TMPDIR") ;
	if( !base ) base = getenv("TEMP") ;
	if( !base ) base = getenv("TMP") ;
	if( !base ) base = "." ;

	char *n1 = (char*)alloca( strlen(base) + strlen(suffix) + 1 ) ;
	char *n2 = n1 ;
	while( *base ) *n2++ = *base++ ;
	while( *suffix ) *n2++ = *suffix++ ;
	*n2 = 0 ;
    int fd = throw_errno_if_minus1( mkstemp( n1 ), "making temp file" ) ;
	throw_errno_if_minus1( unlink( n1 ), "unlinking temp name" ) ;
	name = n1 ;
	return fd ;
}

void enqueue_stream( auto_ptr<Stream>, int = 0 ) ;

void flush_scratch() {
	sort_scratch() ;
	string tempname ;
	int fd = mktempfile( tempname ) ;

	scratch_header.set_is_sorted_by_coordinate( true ) ;
	StreamAdapter< deque< Result* >::const_iterator >
		sa( scratch_header, scratch_space.begin(), scratch_space.end(), scratch_footer ) ;
	clog << "\033[KWriting to tempfile " << tempname << endl ;
	write_stream_to_file( fd, sa ) ;
	throw_errno_if_minus1( lseek( fd, 0, SEEK_SET ), "seeking in ", tempname.c_str() ) ;
	
	auto_ptr<Stream> fs( new AnfoFile( fd, tempname.c_str() ) ) ;
	enqueue_stream( fs ) ;

	for_each( scratch_space.begin(), scratch_space.end(), delete_ptr<Result>() ) ;
	scratch_space.clear() ;
	total_scratch_size = 0 ;
}

void enqueue_stream( auto_ptr<Stream> s, int level ) 
{
	Header h = s->get_header() ;
	if( h.is_sorted_by_coordinate() ) {
		mergeable_queues[ level ].push_back( s.release() ) ;

		if( AnfoFile::num_open_files() > max_que_size ) {
			// get the biggest bin, we'll merge everything below that
			unsigned max_bin = 0 ;
			for( MergeableQueues::const_iterator i = mergeable_queues.begin() ;
					i != mergeable_queues.end() ; ++i ) 
				if( i->second.size() > mergeable_queues[max_bin].size() )
					max_bin = i->first ;

			unsigned total_inputs = 0 ;
			for( MergeableQueues::iterator i = mergeable_queues.begin() ; i->first <= max_bin ; ++i ) 
				total_inputs += i->second.size() ;

			// we must actully make progress, and it must be more than
			// just a single stream to avoid quadratic behaviour (only
			// important in a weird corner case)
			if( total_inputs > 2 ) {
				string fname ;
				int fd = mktempfile( fname ) ;
				clog << "\033[KMerging bins 0.." << max_bin << " to tempfile " << fname << endl ;
				{
					MergeStream ms ;
					for( MergeableQueues::iterator i = mergeable_queues.begin() ; i->first <= max_bin ; ++i ) 
					{
						for( size_t j = 0 ; j != i->second.size() ; ++j )
							ms.add_stream( auto_ptr<Stream>( i->second[j] ) ) ;
						i->second.clear() ;
					}
					write_stream_to_file( fd, ms ) ;
				}

				throw_errno_if_minus1( lseek( fd, 0, SEEK_SET ), "seeking in ", fname.c_str() ) ;
				auto_ptr<Stream> fs( new AnfoFile( fd, fname.c_str() ) ) ;
				enqueue_stream( fs, max_bin + 1 ) ;
			}
		}
	}
	else
	{
		scratch_header.MergeFrom( h ) ;
		for( Result r ; s->read_result( r ) ; ) {
			scratch_space.push_back( new Result( r ) ) ;
			total_scratch_size += r.SpaceUsed() ;
			if( total_scratch_size >= max_arr_size ) flush_scratch() ;
	    }
		scratch_footer.MergeFrom( s->get_footer() ) ;
	}
}

int main_( int argc, const char **argv )
{
	if( argc < 2 ) return 0 ;

    if( !strcmp( argv[1], "-q" ) ) {
        clog.rdbuf(0) ;
        ++argv ;
        --argc ;
    }

	// iterate over command line, glob everything
	glob_t the_glob ;
	glob( argv[1], GLOB_NOSORT, 0, &the_glob ) ;
    for( const char **arg = argv+2 ; arg != argv+argc ; ++arg )
		glob( *arg, GLOB_NOSORT | GLOB_APPEND, 0, &the_glob ) ;

	// iterate over glob results, collect into BestHitStreams
	map< int, BestHitStream* > stream_per_slice ;
	for( char **arg = the_glob.gl_pathv ; 
			arg != the_glob.gl_pathv + the_glob.gl_pathc ; ++arg )
	{
		clog << "\033[KReading from file " << *arg << endl ;
		auto_ptr<Stream> s( new AnfoFile( *arg ) ) ;

		Header h = s->get_header() ;
		if( h.has_sge_slicing_stride() ) {
			BestHitStream* &bhs = stream_per_slice[ h.sge_slicing_index(0) ] ;
			if( !bhs ) bhs = new BestHitStream ;
			bhs->add_stream( s ) ;

			// if a BestHitStream is ready, add it to the queue of
			// mergeable streams
			if( bhs->enough_inputs() ) {
				clog << "\033[KGot everything for slice " << h.sge_slicing_index(0) << endl ;
				enqueue_stream( auto_ptr<Stream>( bhs ) ) ;
				stream_per_slice.erase( h.sge_slicing_index(0) ) ;
			}
		}
		// if no slicing was done, add to queue of mergeable streams
		else enqueue_stream( s ) ;
	}
	if( !stream_per_slice.empty() ) 
		cerr << "\033[KWARNING: input appears to be incomplete" << endl ;

	for( map< int, BestHitStream* >::iterator l = stream_per_slice.begin(), r = stream_per_slice.end() ;
			l != r ; ++l ) enqueue_stream( auto_ptr<Stream>( l->second ) ) ;

	MergeStream final_stream ;

	// at the end, merge everything
	if( scratch_space.begin() != scratch_space.end() ) {
		clog << "\033[Kfinal sort" << endl ;
		sort_scratch() ;
		clog << "\033[Kdone sorting" << endl ;
		auto_ptr<Stream> as( new StreamAdapter< deque< Result* >::const_iterator >(
					scratch_header, scratch_space.begin(), scratch_space.end(), scratch_footer ) ) ;
		final_stream.add_stream( as ) ;
	}

	for( MergeableQueues::const_iterator i = mergeable_queues.begin() ; i != mergeable_queues.end() ; ++i )
		for( MergeableQueue::const_iterator j = i->second.begin() ; j != i->second.end() ; ++j )
			final_stream.add_stream( auto_ptr<Stream>( *j ) ) ;
	
	if( final_stream.get_header().has_version() ) {
		clog << "\033[KMerging everything to output." << endl ;
		RmdupStream rs( &final_stream ) ;
		// int r = write_stream_to_file( 1, final_stream, true ) ;
		int fd = open( "/var/tmp/x.anfo", O_CREAT | O_TRUNC | O_WRONLY ) ;
		int r = write_stream_to_file( fd, rs, true ) ;
		for_each( scratch_space.begin(), scratch_space.end(), delete_ptr<Result>() ) ;
		clog << "\033[KDone" << endl ;
		return r ;
	} else {
		clog << "\033[KInsufficient input, cannot write." << endl ;
		return 1 ;
	}
}


