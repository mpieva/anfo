#include "config.h"

#include "compress_stream.h"
#include "output_streams.h"
#include "stream.h"
#include "util.h"

#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/text_format.h>

#include <algorithm>
#include <deque>
#include <iostream>

#include <glob.h>
#include <popt.h>

using namespace google::protobuf ;
using namespace google::protobuf::io ;
using namespace std ;
using namespace output ;

class NotImplemented : public Exception
{
	private:
		const char *method_ ;

	public:
		NotImplemented( const char* m ) : method_( m ) {}
		virtual void print_to( ostream& s ) const 
		{ s << "method not implemented: " << method_ ; }
} ;

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

namespace streams {

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


//! \brief merges sorted streams into a sorted stream
//! \todo needs to be parameterized with a comparison function (we
//!       already need two).
//! \todo Should be adapted to also combine streams, taking some
//!       functionality from BestHitStream, requiring sorted input, but
//!       not requiring unfiltered input.
class MergeStream : public StreamBundle
{
	private:
		deque< Result > rs_ ;

	public:
		MergeStream() {}
		virtual ~MergeStream() {}

		virtual void add_stream( Stream* s )
		{
			Header h = s->fetch_header() ;
			assert( h.is_sorted_by_coordinate() ) ;
			merge_sensibly( hdr_, h ) ;

			if( s->get_state() == have_output )
			{
				rs_.push_back( s->fetch_result() ) ;
				streams_.push_back( s ) ;
				state_ = have_output ;
			}
			else
			{
				if( state_ == invalid ) state_ = end_of_stream ;
				merge_sensibly( foot_, s->fetch_footer() ) ;
				delete s ;
			}
		}

		virtual Header fetch_header() { return hdr_ ; }
		virtual Footer fetch_footer() { return foot_ ; }
		virtual Result fetch_result() ;
} ;

Result MergeStream::fetch_result() 
{
	int min_idx = 0 ;
	for( size_t i = 1 ; i != rs_.size() ; ++i ) 
		if( by_genome_coordinate()( &rs_[ i ], &rs_[ min_idx ] ) )
			min_idx = i ;

	Result res = rs_[ min_idx ] ;
	Stream *sm = streams_[ min_idx ] ;
	if( sm->get_state() == have_output )
	{
		rs_[ min_idx ] = sm->fetch_result() ; 
	}
	else
	{
		merge_sensibly( foot_, sm->fetch_footer() ) ;
		delete sm ;
		streams_.erase( streams_.begin() + min_idx ) ;
		rs_.erase( rs_.begin() + min_idx ) ;
		if( rs_.empty() ) state_ = end_of_stream ;
	}
	return res ;
}


//! \brief merges multiple streams by taking the best hit
//! If asked for a result, this class takes results from each of the
//! streams in turn, until one sequence has received an entry from each.
//! This sequence is then delivered.  Everything works fine, as long as
//! the sequence names come in in roughly the same order.  Else it still
//! works, but eats memory.
class BestHitStream : public StreamBundle
{
	private:
		typedef map< string, pair< size_t, Result > > Buffer ;
		Buffer buffer_ ;
		size_t nread_, nwritten_, nstreams_ ;
		deque< Stream* >::iterator cur_input_ ;

	public:
		BestHitStream() : nread_(0), nwritten_(0), nstreams_(0) {}
		virtual ~BestHitStream() {}

		//! \brief reads a stream's header and adds the stream as input
		virtual void add_stream( Stream* s ) {
			merge_sensibly( hdr_, s->fetch_header() ) ;
			++nstreams_ ;
			if( s->get_state() == have_output )
			{
				streams_.push_back( s ) ;
				state_ = have_output ;
			}
			else
			{
				if( state_ == invalid ) state_ = end_of_stream ;
				merge_sensibly( foot_, s->fetch_footer() ) ;
				delete s ;
			}
			cur_input_ = streams_.begin() ;
		}

		//! \brief checks if enough inputs are known
		//! Checks if at least as many inputs are known as needed.  It
		//! assumed only one genome index was used and takes the number
		//! of slices declared for it.  (This is mostly useful in
		//! MegaMergeStream to clean up the mess left by a heavily
		//! sliced grid job.)
		bool enough_inputs() const
		{
			if( hdr_.config().policy_size() == 0 ||
					hdr_.config().policy(0).use_compact_index_size() == 0 ) return false ;
			if( !hdr_.config().policy(0).use_compact_index(0).has_number_of_slices() ) return true ;
			return nstreams_ >= hdr_.config().policy(0).use_compact_index(0).number_of_slices() ;
		}

		//! \brief reports final results for one sequence
		//! Enough information is read until one sequence has been
		//! reported on by each input stream.  Everything else is buffered
		//! if necessary.  This works fine as long as the inputs come in
		//! in roughly the same order and are complete.
		virtual Result fetch_result() ;
} ;


Result BestHitStream::fetch_result()
{
	while( !streams_.empty() ) {
		if( cur_input_ == streams_.end() ) cur_input_ = streams_.begin() ;

		Stream *s = *cur_input_ ;
		assert( s->get_state() == have_output ) ;

		Result r = s->fetch_result() ;
		if( s->get_state() == end_of_stream )
		{
			merge_sensibly( foot_, s->fetch_footer() ) ;
			delete s ;
			cur_input_ = streams_.erase( cur_input_ ) ;
		}
		else ++cur_input_ ;

		pair< size_t, Result > &p = buffer_[ r.seqid() ] ;
		++nread_ ;
		++p.first ;
		if( p.second.IsInitialized() ) merge_sensibly( p.second, r ) ;
		else p.second = r ;

		if( p.first == nstreams_ )
		{
			++nwritten_ ;
			Result r1 ;
			swap( r1, p.second ) ;
			buffer_.erase( r.seqid() ) ;

			if( buffer_.empty() && streams_.empty() )
				state_ = end_of_stream ;

			// if( ((nread_ + nwritten_) & 0xFFFF) == 0 ) 
			clog << "\033[KRead " << nread_ << ", delivered "
				<< nwritten_ << ", buffering "
				<< buffer_.size() << '\r' << flush ;

			return r1 ;
		}
	}

	Result r1 ;
	swap( r1, buffer_.begin()->second.second ) ;
	buffer_.erase( buffer_.begin()->first ) ;
	if( buffer_.empty() ) state_ = end_of_stream ;
	return r1 ;
}


//! \brief presents a container as an input stream
//! The container must be of pointers to \c Result, the stream acts as
//! input stream.  A suitable header and footer are supplied at
//! construction time.
template< typename I > class ContainerStream : public Stream
{
	private:
		I cur_, end_ ;

	public:
		ContainerStream( const Header &hdr, I begin, I end, const Footer &foot ) 
			: cur_( begin ), end_( end )
		{
			hdr_ = hdr ;
			foot_ = foot ;
			state_ = cur_ == end_ ? end_of_stream : have_output ;
		}

		virtual Header fetch_header() { return hdr_ ; }
		virtual Footer fetch_footer() { return foot_ ; }
		virtual Result fetch_result() {
			Result r = **cur_ ;
			++cur_ ;
			if( cur_ == end_ ) state_ = end_of_stream ;
			return r ;
		}
} ;

//! \brief stream filter that sorts its input
//! This stream is intended to sort large amounts of data.  To do that,
//! it performs a quick sort on blocks that fit into memory (the maximum
//! size is configured at construction time), writes them out to
//! temporary files, then reads them back in and does a merge sort.  If
//! too many files are open (as counted by \c AnfoReader, again
//! configurable at construction time), some temporary files are merge
//! sorted into a bigger one.
//!
//! \todo We need to parameterize the way to sort (by name, by
//!       coordinate, anything else).

class SortingStream : public Stream
{
	private:
		typedef deque< streams::Stream* > MergeableQueue ;
		typedef map< unsigned, MergeableQueue > MergeableQueues ;
		typedef deque< Result* > ScratchSpace ;

		//! We keep multiple queues of streams ordered and separated by
		//! the number of times they have been merged.  This way we
		//! don't merge the stuff over and over, and instead tend to do
		//! a merge that gives the greatest reduction in number of open
		//! files with least amount of IO.  The key is the number of
		//! times a file has already been merged with others, the value
		//! is just a set of streams.
		MergeableQueues mergeable_queues_ ;
		MergeStream final_stream_ ;

		//! Storage to perform quicksort in.
		ScratchSpace scratch_space_ ;
		// Header scratch_header_ ;
		// Footer scratch_footer_ ;
		int64_t total_scratch_size_ ;

		unsigned max_que_size_, max_arr_size_ ;


		//! \brief quicksort the scratch area
		//! \internal
		void sort_scratch() {
			if( scratch_space_.size() > 1 )
				clog << "\033[KSorting " << scratch_space_.size() << " results in memory" << endl ;
			sort( scratch_space_.begin(), scratch_space_.end(), streams::by_genome_coordinate() ) ;
		}

		void enqueue_stream( streams::Stream*, int = 0 ) ;
		void flush_scratch() ;

	public:
		SortingStream( unsigned as = 256*1024*1024, unsigned qs = 256)
		// SortingStream( unsigned as = 32*1024*1024, unsigned qs = 32 )
			: total_scratch_size_(0), max_que_size_( qs ), max_arr_size_( as )
		{ foot_.set_exit_code( 0 ) ; }

		virtual ~SortingStream()
		{
			for_each( scratch_space_.begin(), scratch_space_.end(), delete_ptr<Result>() ) ;
			for( MergeableQueues::iterator i = mergeable_queues_.begin(), e = mergeable_queues_.end() ; i != e ; ++i )
				for_each( i->second.begin(), i->second.end(), delete_ptr<Stream>() ) ;
		}

		virtual void put_header( const Header& h ) { hdr_ = h ; hdr_.set_is_sorted_by_coordinate(true) ; state_ = need_input ; }
		virtual void put_footer( const Footer& ) ;
		virtual void put_result( const Result& r ) {
			scratch_space_.push_back( new Result( r ) ) ;
			total_scratch_size_ += scratch_space_.back()->SpaceUsed() ;
			if( total_scratch_size_ >= max_arr_size_ ) flush_scratch() ;
		}

		virtual Header fetch_header() { return hdr_ ; }
		virtual Result fetch_result() { Result r = final_stream_.fetch_result() ; state_ = final_stream_.get_state() ; return r ; }
		virtual Footer fetch_footer() { merge_sensibly( foot_, final_stream_.fetch_footer() ) ; return foot_ ; }
} ;


void SortingStream::flush_scratch()
{
	sort_scratch() ;
	string tempname ;
	int fd = mktempfile( &tempname ) ;
	{
		// scratch_header_.set_is_sorted_by_coordinate( true ) ;
		ContainerStream< deque< Result* >::const_iterator >
			sa( hdr_, scratch_space_.begin(), scratch_space_.end(), foot_ ) ;
		AnfoWriter out( fd ) ;
		transfer( sa, out ) ;
		clog << "\033[KWriting to tempfile " << tempname << endl ;
	}
	throw_errno_if_minus1( lseek( fd, 0, SEEK_SET ), "seeking in ", tempname.c_str() ) ;
	enqueue_stream( new AnfoReader( fd, tempname.c_str() ), 1 ) ;

	for_each( scratch_space_.begin(), scratch_space_.end(), delete_ptr<Result>() ) ;
	scratch_space_.clear() ;
	total_scratch_size_ = 0 ;
}

void SortingStream::enqueue_stream( streams::Stream* s, int level ) 
{
	Header h = s->fetch_header() ;
	assert( h.is_sorted_by_coordinate() ) ;
	mergeable_queues_[ level ].push_back( s ) ;

	if( AnfoReader::num_open_files() > max_que_size_ ) {
		// get the biggest bin, we'll merge everything below that
		unsigned max_bin = 0 ;
		for( MergeableQueues::const_iterator i = mergeable_queues_.begin() ;
				i != mergeable_queues_.end() ; ++i ) 
			if( i->second.size() > mergeable_queues_[max_bin].size() )
				max_bin = i->first ;

		unsigned total_inputs = 0 ;
		for( MergeableQueues::iterator i = mergeable_queues_.begin() ; i->first <= max_bin ; ++i ) 
			total_inputs += i->second.size() ;

		// we must actually make progress, and more than just a
		// single stream must be merged to avoid quadratic behaviour
		// (only important in a weird corner case)
		if( total_inputs > 2 ) {
			string fname ;
			int fd = mktempfile( &fname ) ;
			clog << "\033[KMerging bins 0.." << max_bin << " to tempfile " << fname << endl ;
			{
				streams::MergeStream ms ;
				for( MergeableQueues::iterator i = mergeable_queues_.begin() ; i->first <= max_bin ; ++i ) 
				{
					for( size_t j = 0 ; j != i->second.size() ; ++j )
						ms.add_stream( i->second[j] ) ;
					i->second.clear() ;
				}
				AnfoWriter out( fd ) ;
				transfer( ms, out ) ;
			}
			throw_errno_if_minus1( lseek( fd, 0, SEEK_SET ), "seeking in ", fname.c_str() ) ;
			enqueue_stream( new streams::AnfoReader( fd, fname.c_str() ), max_bin + 1 ) ;
		}
	}
}

//! \brief ends the input, initiates sorting
//! Only when the input ends can we completely sort it, so setting the
//! footer switches to output mode.  Here we collect the temporary files
//! we've written and become a \c MergeStream.
void SortingStream::put_footer( const Footer& f ) 
{
	foot_ = f ;

	// if anything's buffered, sort it and include it in the merge
	if( scratch_space_.begin() != scratch_space_.end() ) {
		clog << "\033[Kfinal sort" << endl ;
		sort_scratch() ;
		clog << "\033[Kdone sorting" << endl ;
		final_stream_.add_stream( new ContainerStream< deque< Result* >::const_iterator >(
					hdr_, scratch_space_.begin(), scratch_space_.end(), foot_ ) ) ;
	}

	// add any streams that have piled up
	for( MergeableQueues::const_iterator i = mergeable_queues_.begin() ; i != mergeable_queues_.end() ; ++i )
		for( MergeableQueue::const_iterator j = i->second.begin() ; j != i->second.end() ; ++j )
			final_stream_.add_stream( *j ) ;
	mergeable_queues_.clear() ;
	
	state_ = final_stream_.get_state() ;

	// sanity check if we got a complete header (not sure if this can
	// possibly go wrong any more; probably only gets triggered when
	// there was no input at all)
	// if( !hdr_.has_version() ) throw "SortingStream: insufficient input, cannot produce output" ;
}


class MegaMergeStream : public ConcatStream
{
	private:
		map< int, BestHitStream* > stream_per_slice_ ;

	public:
		MegaMergeStream() {}
		virtual ~MegaMergeStream() {}

		virtual void add_stream( Stream* ) ;
		virtual Header fetch_header() ;
} ;

void MegaMergeStream::add_stream( Stream* s ) 
{
	Header h = s->fetch_header() ;
	if( h.has_sge_slicing_stride() )
	{
		BestHitStream* &bhs = stream_per_slice_[ h.sge_slicing_index(0) ] ;
		if( !bhs ) bhs = new BestHitStream ;
		bhs->add_stream( s ) ;

		if( bhs->enough_inputs() ) {
			clog << "\033[KGot everything for slice " << h.sge_slicing_index(0) << endl ;
			stream_per_slice_.erase( h.sge_slicing_index(0) ) ;
			ConcatStream::add_stream( bhs ) ;
		}
	}
	else ConcatStream::add_stream( s ) ;
}

Header MegaMergeStream::fetch_header()
{
	if( !stream_per_slice_.empty() ) 
		cerr << "\033[KWARNING: input appears to be incomplete" << endl ;

	for( map< int, BestHitStream* >::iterator
			l = stream_per_slice_.begin(),
			r = stream_per_slice_.end() ; l != r ; ++l )
		ConcatStream::add_stream( l->second ) ;

	return ConcatStream::fetch_header() ;
}


class RepairHeaderStream : public Stream
{
	private:
		const char* editor_ ;

	public:
		RepairHeaderStream( const char* e ) : editor_( e ) {}
		virtual ~RepairHeaderStream() {}

		virtual void put_header( const Header& ) ;
		virtual Header fetch_header() { return hdr_ ; }

		virtual void put_result( const Result& r ) { res_ = r ; state_ = have_output ; }
		virtual Result fetch_result() { state_ = need_input ; return res_ ; }

		virtual void put_footer( const Footer& f ) { foot_ = f ; state_ = end_of_stream ; }
		virtual Footer fetch_footer() { return foot_ ; }
} ;

void RepairHeaderStream::put_header( const Header& h ) 
{
	char tmpname[] = "/tmp/anfo_header_XXXXXX" ;
	int fd = mkstemp( tmpname ) ;
	{
		FileOutputStream fos( fd ) ;
		TextFormat::Print( h, &fos ) ;
	}

	string cmd = string( editor_ ? editor_ : getenv("EDITOR") ) + " " + tmpname ;
	for(;;) {
		system( cmd.c_str() ) ;
		lseek( fd, 0, SEEK_SET ) ;
		FileInputStream fis( fd ) ;
		if( TextFormat::Parse( &fis, &hdr_ ) ) break ;
	} 
	throw_errno_if_minus1( unlink( tmpname ), "unlinking", tmpname ) ;
	state_ = need_input ;
}

class FanOut : public StreamBundle
{
	public:
		FanOut() {}
		virtual ~FanOut() {}

		virtual void add_stream( Stream* s ) { streams_.push_back( s ) ; }

		virtual void put_header( const Header& ) ;
		virtual void put_result( const Result& ) ;
		virtual void put_footer( const Footer& ) ;
} ;

void FanOut::put_header( const Header& h )
{
	for( citer i = streams_.begin() ; i != streams_.end() ; ++i )
		(*i)->put_header( h ) ;
	state_ = streams_.front()->get_state() ;
}

void FanOut::put_result( const Result& r )
{
	for( citer i = streams_.begin() ; i != streams_.end() ; ++i )
		(*i)->put_result( r ) ;
	state_ = streams_.front()->get_state() ;
}

void FanOut::put_footer( const Footer& f )
{
	for( citer i = streams_.begin() ; i != streams_.end() ; ++i )
		(*i)->put_footer( f ) ;
	state_ = streams_.front()->get_state() ;
}

class Compose : public StreamBundle
{
	private:
		void update_status() ;

	public:
		Compose() {}
		virtual ~Compose() {}

		virtual void add_stream( Stream* s ) { streams_.push_back( s ) ; }

		virtual void put_header( const Header& ) ;
		virtual void put_result( const Result& ) ;
		virtual void put_footer( const Footer& ) ;

		virtual Header fetch_header() ;
		virtual Result fetch_result() ;
		virtual Footer fetch_footer() ;
} ;

void Compose::put_header( const Header& h_ )
{
	Header h = h_ ;
	citer i = streams_.begin(), e = streams_.end() ;
	for( e-- ; i != e ; ++i )
	{
		(*i)->put_header( h ) ;
		h = (*i)->fetch_header() ;
	}
	(*i)->put_header( h ) ;
	update_status() ;
}

Header Compose::fetch_header() 
{
	citer i = streams_.begin(), e = streams_.end() ;
	Header h = (*i)->fetch_header() ;
	for( ++i ; i != e ; ++i )
	{
		(*i)->put_header( h ) ;
		h = (*i)->fetch_header() ;
	}
	update_status() ;
	return h ;
}

void Compose::put_result( const Result& r )
{
	streams_.front()->put_result( r ) ;
	update_status() ;
}

Result Compose::fetch_result()
{
	Result r = streams_.back()->fetch_result() ;
	update_status() ;
	return r ;
}

void Compose::put_footer( const Footer& f )
{
	streams_.front()->put_footer( f ) ;
	update_status() ;
}

Footer Compose::fetch_footer() 
{
	Footer f = streams_.back()->fetch_footer() ;
	update_status() ;
	return f ;
}

void Compose::update_status()
{
	// look at a stream at a time, starting from the end
	for( criter i = streams_.rbegin() ;; )
	{
		// if we fell of the far end, we need more input
		if( i == streams_.rend() ) { state_ = need_input ; return ; }
		// else consider the state
		state s = (*i)->get_state() ;
		if( s == need_input ) {
			// input comes from previous stream
			++i ;
		}
		else if( s == have_output ) {
			// output's available, either to the outside or to
			// downstream filters
			if( i == streams_.rbegin() ) { state_ = have_output ; return ; }
			Result r = (*i)->fetch_result() ;
			--i ;
			(*i)->put_result( r ) ;
		}
		else if( s == end_of_stream ) {
			// nothing left, pass the footer to see if some data is
			// buffered
			if( i == streams_.rbegin() ) { state_ = end_of_stream ; return ; }
			Footer f = (*i)->fetch_footer() ;
			--i ;
			(*i)->put_footer( f ) ;
		}
		else {
			// easy: something's broken
			state_ = invalid ;
			return ;
		}
	}
}

class StatStream : public Stream
{
	private:
		const char* fn_ ;
		string name_ ;

		unsigned total_, mapped_, mapped_u_ ;
		uint64_t bases_, bases_gc_, bases_m_, bases_gc_m_ ;
		uint64_t bases_squared_, bases_m_squared_ ; 

		void printout( ostream& ) ;

	public:
		StatStream( const char* fn )
			: fn_(fn), total_(0), mapped_(0), mapped_u_(0)
		    , bases_(0), bases_gc_(0), bases_m_(0), bases_gc_m_(0)
		    , bases_squared_(0), bases_m_squared_(0)
		{ state_ = need_input ; }

		virtual ~StatStream() {}

		void put_header( const Header& ) {}
		void put_result( const Result& ) ;
		void put_footer( const Footer& ) ;
} ;

void StatStream::put_result( const Result& r )
{
	unsigned bases = r.sequence().size() ;
	unsigned gc = count( r.sequence().begin(), r.sequence().end(), 'G' )
		        + count( r.sequence().begin(), r.sequence().end(), 'C' ) ;
	++total_ ;
	bases_ += bases ;
	bases_gc_ += gc ;
	bases_squared_ += bases*bases ;
	// XXX use genome?
	if( r.has_best_to_genome() )
	{
		++mapped_ ;
		bases_m_ += bases ;
		bases_m_squared_ += bases*bases ;
		bases_gc_m_ += gc ;
		if( !r.has_diff_to_next() || r.diff_to_next() >= 140 )
		{
			++mapped_u_ ;
		}
	}
	if( name_.empty() ) name_ = r.seqid() ;
}

void StatStream::put_footer( const Footer& )
{
	if( fn_ && strcmp( fn_, "-" ) )
	{
		ofstream s( fn_ ) ;
		printout( s ) ;
	}
	else printout( cout ) ;
}

static inline float std_dev( int n, uint64_t m1, uint64_t m2 )
{ return sqrt( (m2 - m1*m1/n) / (n-1) ) ; }

//! \brief prints statistics
//! The fields are:
//! -# some read name (as an audit trail)
//! -# total number of reads
//! -# number of mapped reads
//! -# number of uniquely mapped reads
//! -# GC content in raw data
//! -# GC content in mapped data
//! -# average length of raw data
//! -# standard deviation of raw data
//! -# average length of mapped data
//! -# standard deviation of length of mapped data
void StatStream::printout( ostream& s )
{
	s << name_ << '\t'
	  << total_ << '\t' << mapped_ << '\t' << mapped_u_ << '\t' 	// reads, mapped reads, uniquely mapped reads
	  << 100*(float)mapped_/(float)total_ << '\t'					// percent hominid
      << 100*(float)bases_gc_ / (float)bases_ << '\t'				// GC content
	  << 100*(float)bases_gc_m_ / (float)bases_m_ << '\t'			// GC content, mapped only
	  << bases_ / total_ << '\t'									// avg. length
	  << std_dev( total_, bases_, bases_squared_ ) << '\t'       	// std.dev. length
	  << bases_m_ / mapped_ << '\t'  								// avg. mapped length
	  << std_dev( mapped_, bases_m_, bases_m_squared_ ) << endl ;   // std.dev. mapped length
} 

} // namespace

//! \page anfo_stream stream-like operations on ANFO files
//!
//! This program builds a chain out of a variety of filters for ANFO
//! streams, allowing various filters, merging and format conversions to
//! be applied.  The options on the command line are read from left to
//! right.  Those that name filters will build a chain of filters to be
//! applied to every input file until a merge operation (merge,
//! mega-merge, or concat) is reached, then they build up a chain to be
//! applied to the result of the merge.  Some filters will write a file
//! containing a subset of the alignments.  Afterwards, all the input
//! files are "sucked" through the filter pipeline, driven by the
//! calculation of a few statistics, which are printed at the end.

using namespace streams ;

template< typename S > struct FilterParams {
	typedef S* (*F)( float, float, const char*, const char* ) ;

	float slope ;
	float intercept ;
	const char* genome ;
	const char* arg ;
	F maker ;
} ;

typedef std::vector< FilterParams< Stream > > FilterStack ;

Stream* mk_sort_by_pos( float, float, const char* genome, const char* arg )
{ return new SortingStream( (arg ? atoi( arg ) : 1024) * 1024 * 1024 ) ; } // XXX use genome 

Stream* mk_sort_by_name( float, float, const char*, const char* )
{ throw NotImplemented( __PRETTY_FUNCTION__ ) ; } // XXX stream is missing

Stream* mk_filter_by_length( float, float, const char*, const char* arg )
{ return new LengthFilter( atoi(arg) ) ; }

Stream* mk_filter_by_score( float s, float i, const char* g, const char* )
{ return new ScoreFilter( s, i, g ) ; }

Stream* mk_filter_by_hit( float, float, const char* genome, const char* arg )
{ return new HitFilter( genome, arg ) ; }

Stream* mk_filter_qual( float, float, const char*, const char* arg )
{ return new QualFilter( atoi( arg ) ) ; }

Stream* mk_filter_multi( float, float, const char*, const char* arg )
{ return new MultiFilter( atoi( arg ) ) ; }

Stream* mk_subsample( float, float, const char*, const char* arg )
{ return new Subsample( atof( arg ) ) ; }

Stream* mk_edit_header( float, float, const char*, const char* arg )
{ return new RepairHeaderStream( arg ) ; }

Stream* mk_rmdup( float, float, const char*, const char* )
{ return new RmdupStream() ; } // XXX use genome? how?

StreamBundle* mk_merge( float, float, const char*, const char* )
{ return new MergeStream() ; } // XXX use genome?

StreamBundle* mk_join( float, float, const char*, const char* )
{ return new BestHitStream() ; } // XXX use genome?

StreamBundle* mk_mega_merge( float, float, const char*, const char* )
{ return new MegaMergeStream() ; }

StreamBundle* mk_concat( float, float, const char*, const char* )
{ return new ConcatStream() ; }

Stream* mk_output      ( float, float, const char*, const char* fn )
{ return 0 == strcmp( fn, "-" ) ? new AnfoWriter( 1, true ) : new AnfoWriter( fn, true ) ; } 

Stream* mk_output_text ( float, float, const char*, const char* fn )
{ return 0 == strcmp( fn, "-" ) ? new TextWriter( 1 ) : new TextWriter( fn ) ; } 

Stream* mk_output_sam  ( float, float, const char*, const char* fn )
{ return 0 == strcmp( fn, "-" ) ? new SamWriter( cout.rdbuf() ) : new SamWriter( fn ) ; } 

Stream* mk_output_fasta( float, float, const char*, const char* fn )
{ return 0 == strcmp( fn, "-" ) ? new FastaWriter( cout.rdbuf() ) : new FastaWriter( fn ) ; } 

Stream* mk_output_table( float, float, const char*, const char* fn )
{ return 0 == strcmp( fn, "-" ) ? new TableWriter( cout.rdbuf() ) : new TableWriter( fn ) ; }

Stream* mk_output_glz  ( float, float, const char*, const char* arg )
{ throw NotImplemented( __PRETTY_FUNCTION__ ) ; } // XXX

Stream* mk_stats       ( float, float, const char*, const char* arg )
{ return new StatStream( arg ) ; }

int main_( int argc, const char **argv )
{
	GOOGLE_PROTOBUF_VERIFY_VERSION ;
	enum { opt_none, opt_sort_pos, opt_sort_name, opt_filter_length,
		opt_filter_score, opt_filter_hit, opt_filter_qual, opt_subsample, opt_filter_multi, opt_edit_header, opt_merge, opt_join,
		opt_mega_merge, opt_concat, opt_rmdup, opt_output, opt_output_sam,
		opt_output_fasta, opt_output_table, opt_output_glz, opt_stats, opt_version, opt_MAX } ;

	FilterParams<Stream>::F filter_makers[opt_MAX] = {
		0, mk_sort_by_pos, mk_sort_by_name, mk_filter_by_length,
		mk_filter_by_score, mk_filter_by_hit, mk_filter_qual, mk_subsample, mk_filter_multi, mk_edit_header, 0, 0,
		0, 0, mk_rmdup, 0, 0,
		0, 0, 0, 0, 0 } ;

	FilterParams<StreamBundle>::F merge_makers[opt_MAX] = {
		0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, mk_merge, mk_join,
		mk_mega_merge, mk_concat, 0, 0, 0,
		0, 0, 0, 0, 0 } ;

	FilterParams<Stream>::F output_makers[opt_MAX] = {
		0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, mk_output, mk_output_sam,
		mk_output_fasta, mk_output_table, mk_output_glz, mk_stats, 0 } ;


	float param_slope = 0, param_intercept = 0 ;
	const char *param_genome = 0 ;
	int param_noise = 1 ;

	int POPT_ARG_OSTR = POPT_ARG_STRING | POPT_ARGFLAG_OPTIONAL ;
	int POPT_ARG_OINT = POPT_ARG_INT | POPT_ARGFLAG_OPTIONAL ;

	struct poptOption options[] = {
		{ "sort-pos",      's', POPT_ARG_OINT,   0, opt_sort_pos,      "sort by alignment position", 0 },
		{ "sort-name",     'S', POPT_ARG_NONE,   0, opt_sort_name,     "sort by read name", 0 },
		{ "filter-length", 'l', POPT_ARG_INT,    0, opt_filter_length, "filter for length of at least L", "L" },
		{ "filter-score",  'f', POPT_ARG_NONE,   0, opt_filter_score,  "filter for max score", 0 },
		{ "filter-hit",    'h', POPT_ARG_OSTR,   0, opt_filter_hit,    "filter for hitting SEQ/anything in G", "SEQ" },
		{ "filter-qual",    0 , POPT_ARG_INT,    0, opt_filter_qual,   "delete bases with quality below Q", "Q" },
		{ "subsample",      0,  POPT_ARG_FLOAT,  0, opt_subsample,     "subsample a fraction F of the results", "F" },
		{ "multiplicity",   0 , POPT_ARG_INT,    0, opt_filter_multi,  "keep reads with multiplicity above N", "N" },
		{ "edit-header",    0 , POPT_ARG_NONE,   0, opt_edit_header,   "invoke editor on the stream's header", 0 },
		{ "concat",        'c', POPT_ARG_NONE,   0, opt_concat,        "concatenate streams", 0 },
		{ "merge",         'm', POPT_ARG_NONE,   0, opt_merge,         "merge sorted streams", 0 },
		{ "join",          'j', POPT_ARG_NONE,   0, opt_join,          "join streams and retain best hits", 0 },
		{ "mega-merge",     0 , POPT_ARG_NONE,   0, opt_mega_merge,    "merge many streams, e.g. from grid", 0 },
		{ "rmdup",         'd', POPT_ARG_NONE,   0, opt_rmdup,         "remove PCR duplicates", 0 },
		{ "output",        'o', POPT_ARG_STRING, 0, opt_output,        "write native stream to file FILE", "FILE" },
		{ "output-sam",     0 , POPT_ARG_STRING, 0, opt_output_sam,    "write alns in sam format to FILE", "FILE" },
		{ "output-fasta",   0 , POPT_ARG_STRING, 0, opt_output_fasta,  "write alignments in fasta format to FILE", "FILE" },
		{ "output-table",   0 , POPT_ARG_STRING, 0, opt_output_table,  "write a table with simple stats to FILE", "FILE" },
		{ "duct-tape",      0 , POPT_ARG_STRING, 0, opt_output_glz,    "not-quite-assemble in glz format into FILE", "FILE" },
		{ "stats",          0,  POPT_ARG_STRING, 0, opt_stats,         "write simple statistics to FILE", "FILE" },

		{ "set-slope",      0 , POPT_ARG_FLOAT,  &param_slope,      0, "set slope for subsequent filters to S", "S" },
		{ "set-intercept",  0 , POPT_ARG_FLOAT,  &param_intercept,  0, "set length intercept for filters to L", "L" },
		{ "set-genome",     0 , POPT_ARG_STRING, &param_genome,     0, "set interesting genome for most operations to G", "G" },
		{ "clear-genome",   0 , POPT_ARG_VAL,    &param_genome,     0, "clear interesting genome", 0 },

		{ "quiet",         'q', POPT_ARG_VAL,    &param_noise,      0, "suppress most output", 0 },
		{ "verbose",       'v', POPT_ARG_VAL,    &param_noise,      2, "produce more output", 0 },
		{ "version",       'V', POPT_ARG_NONE,   0, opt_version,       "print version number and exit", 0 },
		POPT_AUTOHELP POPT_TABLEEND
	} ;

	FilterStack filters_initial ;
	FilterStack *filters_current = &filters_initial ;

	std::auto_ptr< StreamBundle > merging_stream( new ConcatStream ) ;

	typedef std::deque< FilterStack > FilterStacks ;
	FilterStacks filters_terminal ;
	filters_terminal.push_back( FilterStack() ) ;

	poptContext pc = poptGetContext( "anfo", argc, argv, options, 0 ) ;
	poptSetOtherOptionHelp( pc, "[OPTION...] [sequence-file...]" ) ;
	poptReadDefaultConfig( pc, 0 ) ;

	if( argc <= 1 ) { poptPrintHelp( pc, stderr, 0 ) ; return 1 ; }
	for( int rc = poptGetNextOpt( pc ) ; rc != -1 ; rc = poptGetNextOpt(pc) )
	{
		if( rc == opt_version ) {
			std::cout << poptGetInvocationName(pc) << ", revision " << PACKAGE_VERSION << std::endl ;
			return 0 ;
		}
		else if( rc >= 0 && filter_makers[rc] )
		{
			FilterParams< Stream > fp = {
				param_slope, param_intercept, param_genome, 
				poptGetOptArg( pc ), filter_makers[rc]
			} ;
			filters_current->push_back( fp ) ;
		}
		else if( rc >= 0 && merge_makers[rc] )
		{
			// make sure we are still creating input filters
			if( filters_current != &filters_initial )
				throw "merge-like commands cannot not follow merge- or output-like commands" ;

			merging_stream.reset( (merge_makers[rc])(
						param_slope, param_intercept, param_genome, poptGetOptArg( pc ) ) ) ;

			// from now on we build output filters
			filters_current = &filters_terminal.back() ;
		}
		else if( rc >= 0 && output_makers[rc] )
		{
			// make sure we are creating output filters
			if( filters_current == &filters_initial )
				filters_current = &filters_terminal.back() ;

			// create filter
			FilterParams< Stream > fp = {
				param_slope, param_intercept, param_genome,
				poptGetOptArg( pc ), output_makers[rc]
			} ;
			filters_current->push_back( fp ) ;

			// start new output stream
			filters_terminal.push_back( FilterStack() ) ;
			filters_current = &filters_terminal.back() ;
		}
		else
		{
			std::clog << poptGetInvocationName(pc) << ": " << poptStrerror( rc ) 
				<< ' ' << poptBadOption( pc, 0 ) << std::endl ;
			return 1 ; 
		}
	}

	// iterate over non-option arguments, glob everything
	glob_t the_glob ;
	int glob_flag = GLOB_NOSORT ;

	while( const char* arg = poptGetArg( pc ) )
	{
		glob( arg, glob_flag, 0, &the_glob ) ;
		glob_flag |= GLOB_APPEND ;
	}

	// iterate over glob results
	for( char **arg = the_glob.gl_pathv ; arg != the_glob.gl_pathv + the_glob.gl_pathc ; ++arg )
	{
		auto_ptr< Compose > c( new Compose ) ;
		c->add_stream( new AnfoReader( *arg ) ) ;
		for( FilterStack::const_iterator i = filters_initial.begin() ; i != filters_initial.end() ; ++i )
			c->add_stream( (i->maker)( i->intercept, i->slope, i->genome, i->arg ) ) ;
		merging_stream->add_stream( c.release() ) ;
	}

	FanOut out ;
	// only one output and that one is empty?  add a writer for stdout
	if( filters_terminal.size() == 1 && filters_terminal[0].empty() )
	{
		out.add_stream( new TextWriter( 1 ) ) ;
	}
	else for( FilterStacks::const_iterator i = filters_terminal.begin() ; i != filters_terminal.end() ; ++i )
	{
		if( !i->empty() ) {
			auto_ptr< Compose > c( new Compose ) ;
			for( FilterStack::const_iterator j = i->begin() ; j != i->end() ; ++j )
				c->add_stream( (j->maker)( j->intercept, j->slope, j->genome, j->arg ) ) ;
			out.add_stream( c.release() ) ;
		}
	}

	transfer( *merging_stream, out ) ;
	poptFreeContext( pc ) ;
	return 0 ;
}

