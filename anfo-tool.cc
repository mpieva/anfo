#include "config.h"

#include "compress_stream.h"
#include "ducttape.h"
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

namespace streams {

//! \brief compares hits by smallest genome coordinate
//! Comparison is first done lexically on the subject name, then on the
//! smallest coordinate that's part of the alignment, then on the
//! alignment length.  A particular genome can be selected, if this
//! isn't done, we sort on the best hit and compare the genome name even
//! before the subject name.
struct by_genome_coordinate {
	const char *g_ ;
	by_genome_coordinate( const char *g ) : g_(g) {}

	bool operator() ( const Result *a, const Result *b ) {
		if( has_hit_to( *a, g_ ) && !has_hit_to( *b, g_ ) ) return true ;
		if( !has_hit_to( *a, g_ ) ) return false ;

		const Hit& u = hit_to( *a, g_ ), v = hit_to( *b, g_ ) ;
		if( !g_ ) {
			if( u.genome_name() < v.genome_name() ) return true ;
			if( v.genome_name() < u.genome_name() ) return false ;
		}
		if( u.sequence() < v.sequence() ) return true ;
		if( v.sequence() < u.sequence() ) return false ;
		if( u.start_pos() < v.start_pos() ) return true ;
		if( v.start_pos() < u.start_pos() ) return false ;
		return u.aln_length() < v.aln_length() ;
	}

	void tag_header( output::Header& h ) {
		h.clear_is_sorted_by_name() ;
		h.set_is_sorted_by_coordinate( g_ ? g_ : "" ) ;
	}
} ;

struct by_seqid {
	bool operator() ( const Result *a, const Result *b ) {
		return a->read().seqid() < b->read().seqid() ;
	}
	void tag_header( output::Header& h ) {
		h.clear_is_sorted_by_coordinate() ;
		h.set_is_sorted_by_name( true ) ;
	}
} ;


//! \brief merges sorted streams into a sorted stream
//! What to compare on is read from the input streams' header.  If they
//! are unsorted, we fail.  Else we check that they are sorted in the
//! same way and merge accordingly.
class MergeStream : public StreamBundle
{
	private:
		deque< Result > rs_ ;
		enum { unknown, by_name, by_coordinate } mode_ ;
		const char *g_ ;

	public:
		MergeStream() : mode_( unknown ), g_(0) {}
		virtual ~MergeStream() { free( const_cast<char*>( g_ ) ) ; }

		virtual void add_stream( Stream* s )
		{
			Header h = s->fetch_header() ;
			if( streams_.empty() ) hdr_ = h ;
			else merge_sensibly( hdr_, h ) ;

			if( h.is_sorted_by_name() ) {
				if( mode_ == unknown ) mode_ = by_name ;
				else if( mode_ != by_name ) 
					throw "MergeStream: inconsistent sorting of input" ;
			}
			else if( h.has_is_sorted_by_coordinate() ) {
				if( mode_ == unknown ) {
					mode_ = by_coordinate ;
					g_ = strdup( h.is_sorted_by_coordinate().c_str() ) ;
				}
				else if( mode_ != by_coordinate || g_ != h.is_sorted_by_coordinate() )
					throw "MergeStream: inconsistent sorting of input" ;
			}

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

		virtual Header fetch_header()
		{
			if( mode_ == unknown ) throw "MergeStream: don't know what to merge on" ;
			return hdr_ ;
		}
		virtual Result fetch_result() ;
} ;

Result MergeStream::fetch_result() 
{
	int min_idx = 0 ;
	for( size_t i = 1 ; i != rs_.size() ; ++i ) 
		if( ( mode_ == by_coordinate && by_genome_coordinate( *g_ ? g_ : 0 )( &rs_[ i ], &rs_[ min_idx ] ) )
				|| ( mode_ == by_name && by_seqid()( &rs_[ i ], &rs_[ min_idx ] ) ) )
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
		Chan progress_ ;

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

		pair< size_t, Result > &p = buffer_[ r.read().seqid() ] ;
		++nread_ ;
		++p.first ;
		if( p.second.IsInitialized() ) merge_sensibly( p.second, r ) ;
		else p.second = r ;

		if( p.first == nstreams_ )
		{
			++nwritten_ ;
			Result r1 ;
			swap( r1, p.second ) ;
			buffer_.erase( r.read().seqid() ) ;

			if( buffer_.empty() && streams_.empty() )
				state_ = end_of_stream ;

			if( ((nread_ + nwritten_) & 0xFFFF) == 0 )
			{
				stringstream s ;
				s << "BestHitStream: " << nread_ << "in "
				<< nwritten_ << "out "
				<< buffer_.size() << "buf" ;
				progress_( Console::info, s.str() ) ;
			}
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

		virtual Result fetch_result() {
			Result r = **cur_ ;
			++cur_ ;
			if( cur_ == end_ ) state_ = end_of_stream ;
			return r ;
		}
} ;

unsigned SortingStream__ninstances = 0 ;

//! \brief stream filter that sorts its input
//! This stream is intended to sort large amounts of data.  To do that,
//! it performs a quick sort on blocks that fit into memory (the maximum
//! size is configured at construction time), writes them out to
//! temporary files, then reads them back in and does a merge sort.  If
//! too many files are open (as counted by \c AnfoReader, again
//! configurable at construction time), some temporary files are merge
//! sorted into a bigger one.

template <class Comp> class SortingStream : public Stream
{
	private:
		typedef deque< streams::Stream* > MergeableQueue ;
		typedef map< unsigned, MergeableQueue > MergeableQueues ;
		typedef deque< Result* > ScratchSpace ;

		//! Storage to perform quicksort in.
		ScratchSpace scratch_space_ ;

		//! We keep multiple queues of streams ordered and separated by
		//! the number of times they have been merged.  This way we
		//! don't merge the stuff over and over, and instead tend to do
		//! a merge that gives the greatest reduction in number of open
		//! files with least amount of IO.  The key is the number of
		//! times a file has already been merged with others, the value
		//! is just a set of streams.
		MergeableQueues mergeable_queues_ ;
		MergeStream final_stream_ ;

		int64_t total_scratch_size_ ;

		unsigned max_que_size_, max_arr_size_ ;
		Comp comp_ ;

		//! \brief quicksort the scratch area
		//! \internal
		void sort_scratch() {
			if( scratch_space_.size() > 1 )
			{
				std::stringstream s ;
				s << "SortingStream: qsorting " << scratch_space_.size() << " results" ; 
				console.output( Console::notice, s.str() ) ;
			}
			sort( scratch_space_.begin(), scratch_space_.end(), comp_ ) ;
		}

		void enqueue_stream( streams::Stream*, int = 0 ) ;
		void flush_scratch() ;

	public:
		SortingStream( unsigned as = 256*1024*1024, unsigned qs = 256, Comp comp = Comp() )
			: total_scratch_size_(0), max_que_size_( qs ), max_arr_size_( as ), comp_( comp )
		{ foot_.set_exit_code( 0 ) ; ++SortingStream__ninstances ; }

		virtual ~SortingStream()
		{
			for_each( scratch_space_.begin(), scratch_space_.end(), delete_ptr<Result>() ) ;
			for( MergeableQueues::iterator i = mergeable_queues_.begin(), e = mergeable_queues_.end() ; i != e ; ++i )
				for_each( i->second.begin(), i->second.end(), delete_ptr<Stream>() ) ;
			--SortingStream__ninstances ;
		}

		virtual void put_header( const Header& h ) { Stream::put_header( h ) ; comp_.tag_header( hdr_ ) ; }
		virtual void put_footer( const Footer& ) ;
		virtual void put_result( const Result& r ) {
			scratch_space_.push_back( new Result( r ) ) ;
			total_scratch_size_ += scratch_space_.back()->SpaceUsed() ;
			if( total_scratch_size_ >= max_arr_size_ ) flush_scratch() ;
		}

		virtual Result fetch_result() { Result r = final_stream_.fetch_result() ; state_ = final_stream_.get_state() ; return r ; }
		virtual Footer fetch_footer() { merge_sensibly( foot_, final_stream_.fetch_footer() ) ; return foot_ ; }
} ;

template < typename Comp > void SortingStream<Comp>::flush_scratch()
{
	sort_scratch() ;
	string tempname ;
	int fd = mktempfile( &tempname ) ;
	{
		ContainerStream< deque< Result* >::const_iterator >
			sa( hdr_, scratch_space_.begin(), scratch_space_.end(), foot_ ) ;
		AnfoWriter out( fd, tempname.c_str() ) ;
		console.output( Console::notice, "SortingStream: Writing to tempfile " + tempname ) ;
		transfer( sa, out ) ;
	}
	throw_errno_if_minus1( lseek( fd, 0, SEEK_SET ), "seeking in ", tempname.c_str() ) ;
	enqueue_stream( make_input_stream( /*new AnfoReader(*/ fd, tempname.c_str() ), 1 ) ;

	for_each( scratch_space_.begin(), scratch_space_.end(), delete_ptr<Result>() ) ;
	scratch_space_.clear() ;
	total_scratch_size_ = 0 ;
}

template < typename Comp > void SortingStream<Comp>::enqueue_stream( streams::Stream* s, int level ) 
{
	Header h = s->fetch_header() ;
	assert( h.has_is_sorted_by_coordinate() ) ;
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
			std::stringstream s ;
			s << "SortingStream: Merging bins 0.." << max_bin << " to tempfile " << fname ;
			console.output( Console::notice, s.str() ) ;
			{
				streams::MergeStream ms ;
				for( MergeableQueues::iterator i = mergeable_queues_.begin() ; i->first <= max_bin ; ++i ) 
				{
					for( size_t j = 0 ; j != i->second.size() ; ++j )
						ms.add_stream( i->second[j] ) ;
					i->second.clear() ;
				}
				AnfoWriter out( fd, fname.c_str() ) ;
				transfer( ms, out ) ;
			}
			throw_errno_if_minus1( lseek( fd, 0, SEEK_SET ), "seeking in ", fname.c_str() ) ;
			enqueue_stream( make_input_stream( /*new streams::AnfoReader(*/ fd, fname.c_str() ), max_bin + 1 ) ;
		}
	}
}

//! \brief ends the input, initiates sorting
//! Only when the input ends can we completely sort it, so setting the
//! footer switches to output mode.  Here we collect the temporary files
//! we've written and become a \c MergeStream.
template < typename Comp > void SortingStream<Comp>::put_footer( const Footer& f ) 
{
	Stream::put_footer( f ) ;

	// We have to be careful about buffering; if more than one
	// SortingStream is active, we could run out of RAM.  Therefore, if
	// we're alone, we sort and add a a stream.  Else we flush to
	// temporary storage.

	if( scratch_space_.begin() != scratch_space_.end() )
	{
		if( SortingStream__ninstances > 1 ) flush_scratch() ; 
		else {
			console.output( Console::notice, "SortingStream: final sort" ) ;
			sort_scratch() ;
			final_stream_.add_stream( new ContainerStream< deque< Result* >::const_iterator >(
						hdr_, scratch_space_.begin(), scratch_space_.end(), foot_ ) ) ;
		}
	}

	// add any streams that have piled up
	for( MergeableQueues::const_iterator i = mergeable_queues_.begin() ; i != mergeable_queues_.end() ; ++i )
		for( MergeableQueue::const_iterator j = i->second.begin() ; j != i->second.end() ; ++j )
			final_stream_.add_stream( *j ) ;
	mergeable_queues_.clear() ;
	console.output( Console::notice, "SortingStream: merging everything to output" ) ;
	
	final_stream_.fetch_header() ;
	state_ = final_stream_.get_state() ;
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
			std::stringstream s ;
			s << "MegaMergeStream: got everything for slice " << h.sge_slicing_index(0) ;
			console.output( Console::info, s.str() ) ;
			stream_per_slice_.erase( h.sge_slicing_index(0) ) ;
			ConcatStream::add_stream( bhs ) ;
		}
	}
	else ConcatStream::add_stream( s ) ;
}

Header MegaMergeStream::fetch_header()
{
	if( !stream_per_slice_.empty() ) 
		console.output( Console::warning, "MegaMergeStream: input appears to be incomplete" ) ;

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
		if( system( cmd.c_str() ) ) break ;;
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
	Stream::put_footer( f ) ;
	for( citer i = streams_.begin() ; i != streams_.end() ; ++i )
		(*i)->put_footer( f ) ;
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
		const char* g_ ;
		string name_ ;

		unsigned total_, mapped_, mapped_u_, different_ ;
		uint64_t bases_, bases_gc_, bases_m_, bases_gc_m_ ;
		uint64_t bases_squared_, bases_m_squared_ ; 

		void printout( ostream&, bool ) ;

	public:
		StatStream( const char* fn, const char* g )
			: fn_(fn), g_(g), total_(0), mapped_(0), mapped_u_(0), different_(0)
		    , bases_(0), bases_gc_(0), bases_m_(0), bases_gc_m_(0)
		    , bases_squared_(0), bases_m_squared_(0) {}
		virtual ~StatStream() {}

		virtual void put_result( const Result& ) ;
		virtual void put_footer( const Footer& ) ;
} ;

void StatStream::put_result( const Result& r )
{
	unsigned bases = r.read().sequence().size() ;
	unsigned gc = count( r.read().sequence().begin(), r.read().sequence().end(), 'G' )
		        + count( r.read().sequence().begin(), r.read().sequence().end(), 'C' ) ;
	unsigned count = std::max( r.member_size(), 1 ) ;
	total_ += count ;
	bases_ += bases ;
	bases_gc_ += gc ;
	bases_squared_ += bases*bases ;
	if( has_hit_to( r, g_ ) )
	{
		mapped_ += count ;
		bases_m_ += bases ;
		bases_m_squared_ += bases*bases ;
		bases_gc_m_ += gc ;
		if( !hit_to( r, g_ ).has_diff_to_next() || hit_to( r, g_ ).diff_to_next() >= 60 )
		{
			mapped_u_ += count ;
			++different_ ;
		}
	}
	if( name_.empty() ) name_ = r.read().seqid() ;
}

void StatStream::put_footer( const Footer& f )
{
	Stream::put_footer( f ) ;
	if( !fn_ || 0 == strcmp( fn_, "-" ) ) printout( cout, true ) ;
	else if( 0 == strcmp( fn_, "+-" ) ) printout( cout, false ) ;
	else if( *fn_ == '+' ) 
	{
		ofstream s( fn_+1, ios_base::app ) ;
		printout( s, false ) ;
	}
	else 
	{
		ofstream s( fn_, ios_base::trunc ) ;
		printout( s, true ) ;
	}
}

static inline float std_dev( int n, uint64_t m1, uint64_t m2 )
{ return sqrt( (m2 - m1*m1 / (float)n) / (n-1) ) ; }

//! \brief prints statistics
//! The fields are:
//! -# some read name (as an audit trail)
//! -# total number of reads
//! -# number of mapped reads
//! -# number of uniquely mapped reads
//! -# percentage of mapped reads
//! -# number of unique sequences mapped uniquely
//! -# GC content in raw data
//! -# GC content in mapped data
//! -# average length of raw data
//! -# standard deviation of raw data
//! -# average length of mapped data
//! -# standard deviation of length of mapped data
void StatStream::printout( ostream& s, bool h )
{
	s << (h?"name\t#total\t#mapped\t#mapuniq\t%mapped\t#distinct\t"
		    "GCraw\tGCmapped\trawlen\tdevrawlen\tmaplen\tdevmaplen\n":"")
	  << name_ << '\t' << total_ << '\t' 							// arbitrary name, reads
	  << mapped_ << '\t' << mapped_u_ << '\t' 						// mapped reads, uniquely mapped reads
	  << 100*(float)mapped_/(float)total_ << '\t'					// percent hominid
	  << different_ << '\t'
      << 100*(float)bases_gc_ / (float)bases_ << '\t'				// GC content
	  << 100*(float)bases_gc_m_ / (float)bases_m_ << '\t'			// GC content, mapped only
	  << bases_ / (float)total_ << '\t'								// avg. length
	  << std_dev( total_, bases_, bases_squared_ ) << '\t'       	// std.dev. length
	  << bases_m_ / (float)mapped_ << '\t'  						// avg. mapped length
	  << std_dev( mapped_, bases_m_, bases_m_squared_ ) << endl ;   // std.dev. mapped length
} 

class IgnoreHit : public Filter
{
	private:
		const char* g_ ;
		const char* s_ ;

	public:
		IgnoreHit( const char* g, const char* s ) : g_(g), s_(s) {}
		virtual ~IgnoreHit() {}
		virtual bool xform( Result& r ) 
		{
			int j = 0 ;
			for( int i = 0 ; i != r.hit_size() ; ++i )
			{
				if( (g_ && *g_ && r.hit(i).genome_name() != g_)
						|| (s_ && *s_ && r.hit(i).sequence() != s_) )
				{
					swap( *r.mutable_hit(j), *r.mutable_hit(i) ) ;
					++j ;
				}
			}
			while( j != r.hit_size() ) r.mutable_hit()->RemoveLast() ;
			return true ;
		}
} ;

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

struct ParamBlock {
	float slope ;
	float intercept ;
	int context ;
	const char* genome ;
	const char* arg ;

	ParamBlock( float s, float i, int c, const char* g, const char* a )
		: slope(s), intercept(i), context(c), genome(g), arg(a) {}
} ;

typedef void (*G)( ostream&, const ParamBlock& ) ;

template< typename S > struct FilterParams : public ParamBlock {
	typedef S* (*F)( const ParamBlock& ) ;

	F maker ;
	G describe ;

	FilterParams( const ParamBlock& p, const char* a, F m, G d )
		: ParamBlock(p.slope, p.intercept, p.context, p.genome, a), maker(m), describe(d) {}
} ;

typedef std::vector< FilterParams< Stream > > FilterStack ;

float parse_float( const char* a )
{ 
	char *e ;
	float f = strtod( a, &e ) ;
	if( a && *a && !*e ) return f ;
	throw "expected real number, but found \"" + string(a) + "\"" ;
}
int parse_int( const char* a )
{ 
	char *e ;
	int i = strtol( a, &e, 10 ) ;
	if( a && *a && !*e ) return i ;
	throw "expected integer, but found \"" + string(a) + "\"" ;
}
int parse_int( const char* a, int d )
{ 
	char *e ;
	if( !a || !*a ) return d ;
	int i = strtol( a, &e, 10 ) ;
	if( !*e ) return i ;
	throw "expected integer, but found \"" + string(a) + "\"" ;
}

bool is_stdout( const char* a ) { return !a || 0 == strcmp( a, "-" ) ; }
const char* parse_fn( const char* a ) { return is_stdout( a ) ? "<stdout>" : a ; }

Stream* mk_sort_by_pos( const ParamBlock& p )
{ return new SortingStream<by_genome_coordinate>( parse_int( p.arg, 1024 ) * 1024 * 1024, 256, by_genome_coordinate(p.genome) ) ; }

void desc_sort_by_pos( ostream& ss, const ParamBlock& p )
{ 
	ss << "sort by position on " << ( p.genome ? p.genome : "any genome") 
	   << ", using " << + parse_int( p.arg, 1024 ) <<  " MB" ;
}

Stream* mk_sort_by_name( const ParamBlock& p )
{ return new SortingStream<by_seqid>( parse_int( p.arg, 1024 ) * 1024 * 1024 ) ; }

void desc_sort_by_name( ostream& ss, const ParamBlock& p )
{ ss << "sort by sequence id, using " << parse_int( p.arg, 1024 ) << " MB" ; }

Stream* mk_filter_by_length( const ParamBlock& p )
{ return new LengthFilter( parse_int( p.arg ) ) ; }

void desc_filter_by_length( ostream& ss, const ParamBlock& p )
{ ss << "remove alignments shorter than " << parse_int( p.arg ) ; }

Stream* mk_filter_by_score( const ParamBlock& p )
{ return new ScoreFilter( p.slope, p.intercept, p.genome ) ; }

void desc_filter_by_score( ostream& ss, const ParamBlock& p )
{
	ss << "remove alignments to " << (p.genome?p.genome:"any genome") 
		<< " scoring worse than ( " << p.slope << " * ( L - " << p.intercept << " ) )" ;
}

Stream* mk_filter_by_mapq( const ParamBlock& p )
{ return new MapqFilter( p.genome, parse_int( p.arg ) ) ; }

void desc_filter_by_mapq( ostream& ss, const ParamBlock& p )
{
	ss << "remove alignments where MAPQ" ;
	if( p.genome ) ss << " on genome " << p.genome ;
	ss << " is below " << parse_int( p.arg ) ;
}

Stream* mk_filter_by_hit( const ParamBlock& p )
{ return new HitFilter( p.genome, p.arg ) ; }

void desc_filter_by_hit( ostream& ss, const ParamBlock& p )
{
	ss << "remove sequences without hit" ;
	if( p.arg && *p.arg ) ss << " to sequence " << p.arg ;
	if( p.genome && *p.genome ) ss << " in genome " << p.genome ; 
}

Stream* mk_delete_hit( const ParamBlock& p )
{ return new IgnoreHit( p.genome, p.arg ) ; }

void desc_delete_hit( ostream& ss, const ParamBlock& p )
{
	ss << "delete hits" ;
	if( p.arg && *p.arg ) ss << " to sequence " << p.arg ;
	if( p.genome && *p.genome ) ss << " in genome " << p.genome ;
}

Stream* mk_filter_qual( const ParamBlock& p )
{ return new QualFilter( parse_int( p.arg ) ) ; }

void desc_filter_qual( ostream& ss, const ParamBlock& p )
{ ss << "mask bases with quality below " << parse_int( p.arg ) ; }

Stream* mk_filter_multi( const ParamBlock& p )
{ return new MultiFilter( parse_int( p.arg, 2 ) ) ; }

void desc_filter_multi( ostream& ss, const ParamBlock& p )
{ ss << "retain only sequences that were seen at least " << parse_int( p.arg, 2 ) << " times" ; }

Stream* mk_subsample( const ParamBlock& p )
{ return new Subsample( parse_float( p.arg ) ) ; }

void desc_subsample( ostream& ss, const ParamBlock& p )
{ ss << "subsample a " << parse_float(p.arg) << " fraction of sequences" ; }

Stream* mk_edit_header( const ParamBlock& p )
{ return new RepairHeaderStream( p.arg ) ; }

void desc_edit_header( ostream& ss, const ParamBlock& p )
{ ss << "invoke " << (p.arg?p.arg:" text editor ") << " on stream's header" ; }

Stream* mk_rmdup( const ParamBlock& p )
{ return new RmdupStream( p.slope, p.intercept ) ; }

void desc_rmdup( ostream& ss, const ParamBlock& p )
{
	ss << "coalesce duplicates as long as score is no worse than ( "
		<< p.slope << " * ( L - " << p.intercept << " ) )" ;
}

StreamBundle* mk_merge( const ParamBlock& )
{ return new MergeStream() ; }

void desc_merge( ostream& ss, const ParamBlock& )
{ ss << "merge sorted streams" ; }

StreamBundle* mk_join( const ParamBlock& )
{ return new BestHitStream() ; }

void desc_join( ostream& ss, const ParamBlock& )
{ ss << "join near-sorted streams and retain best hits to each genome" ; }

StreamBundle* mk_mega_merge( const ParamBlock& )
{ return new MegaMergeStream() ; }

void desc_mega_merge( ostream& ss, const ParamBlock& )
{ ss << "join fragments from grid jobs and retain best hits" ; }

StreamBundle* mk_concat( const ParamBlock& )
{ return new ConcatStream() ; }

void desc_concat( ostream& ss, const ParamBlock& )
{ ss << "concatenate streams" ; }

Stream* mk_output( const ParamBlock& p )
{ return is_stdout( p.arg ) ? new AnfoWriter( 1, "<stdout>", true ) : new AnfoWriter( p.arg, true ) ; } 

void desc_output( ostream& ss, const ParamBlock& p )
{ ss << "write native output to " << parse_fn( p.arg ) ; }

Stream* mk_output_text( const ParamBlock& p )
{ return is_stdout( p.arg ) ? new TextWriter( 1 ) : new TextWriter( p.arg ) ; } 

void desc_output_text( ostream& ss, const ParamBlock& p )
{ ss << "write in text format to " << parse_fn( p.arg ) ; }

Stream* mk_output_sam( const ParamBlock& p )
{ return is_stdout( p.arg ) ? new SamWriter( cout.rdbuf(), p.genome ) : new SamWriter( p.arg, p.genome ) ; } 

void desc_output_sam( ostream& ss, const ParamBlock& p )
{ 
	ss << "write alignments" ;
	if( p.genome ) ss << " to genome " << p.genome ;
	ss << " in SAM format to " << parse_fn( p.arg ) ;
}

Stream* mk_output_glz  ( const ParamBlock& p )
{ return is_stdout( p.arg ) ? new GlzWriter( 1 ) : new GlzWriter( p.arg ) ; } 

void desc_output_glz( ostream& ss, const ParamBlock& p )
{ ss << "write contigs in GLZ format to " << parse_fn( p.arg ) ; }

Stream* mk_output_3aln  ( const ParamBlock& p )
{ return new ThreeAlnWriter( p.arg ) ; } // XXX is_stdout( p.arg ) ? new GlzWriter( 1 ) : new GlzWriter( p.arg ) ; } 

void desc_output_3aln( ostream& ss, const ParamBlock& p )
{ ss << "write contigs in 3ALN format to " << p.arg ; } // XXX parse_fn( p.arg ) ; }

Stream* mk_output_fasta( const ParamBlock& p )
{
	return is_stdout( p.arg ) ? new FastaAlnWriter( cout.rdbuf(), p.genome, p.context )
	                          : new FastaAlnWriter( p.arg, p.genome, p.context ) ; 
} 

void desc_output_fasta( ostream& ss, const ParamBlock& p )
{ 
	ss << "write alignments(!) to " << (p.genome?p.genome:"any genome") 
	   << " in FASTA format to " << parse_fn( p.arg ) ;
	if( p.context ) ss << " with " << p.context << "nt of context" ;
}

Stream* mk_output_fastq( const ParamBlock& p )
{ return is_stdout( p.arg ) ? new FastqWriter( cout.rdbuf() ) : new FastqWriter( p.arg ) ; } 

void desc_output_fastq( ostream& ss, const ParamBlock& p )
{ ss << "write sequences(!) in FASTQ format to " << parse_fn( p.arg ) ; }

Stream* mk_output_table( const ParamBlock& p )
{ return is_stdout( p.arg ) ? new TableWriter( cout.rdbuf(), p.genome ) : new TableWriter( p.arg, p.genome ) ; }

void desc_output_table( ostream& ss, const ParamBlock& p )
{ 
	ss << "write useless table" ;
	if( p.genome ) ss << " about genome " << p.genome ;
	ss << " to " << parse_fn( p.arg ) ;
}

Stream* mk_duct_tape( const ParamBlock& p )
{ return new DuctTaper( p.genome, parse_int( p.arg ) ) ; }

void desc_duct_tape( ostream& ss, const ParamBlock& p )
{ 
	ss << "mock-assemble hits" ;
	if( p.genome ) ss << " to genome " << p.genome ;
	ss << " and clamp Q-score to no more than " << parse_int( p.arg ) ;
}

Stream* mk_stats( const ParamBlock& p )
{ return new StatStream( p.arg, p.genome ) ; }

void desc_stats( ostream& ss, const ParamBlock& p )
{ 
	if( p.arg && *p.arg == '+' )
		ss << "append statistics to " << parse_fn( p.arg+1 ) ;
	else
		ss << "write statistics to " << parse_fn( p.arg ) ;
}

const char *poptGetOptArg1( poptContext con )
{
	const char *p = poptGetOptArg( con ) ;
	if( !p || *p != '-' || !p[1] ) return p ;

	console.output( Console::warning, string("poptGetOptArg: ") + p + " treated as parameter" ) ;
	return p ;
}

int main_( int argc, const char **argv )
{
	GOOGLE_PROTOBUF_VERIFY_VERSION ;
	enum { opt_none, opt_sort_pos, opt_sort_name, opt_filter_length,
		opt_filter_score, opt_filter_mapq, opt_filter_hit, opt_delete_hit, opt_filter_qual, opt_subsample, opt_filter_multi, opt_edit_header, opt_merge, opt_join,
		opt_mega_merge, opt_concat, opt_rmdup, opt_output, opt_output_text, opt_output_sam, opt_output_glz, opt_output_3aln,
		opt_output_fasta, opt_output_fastq, opt_output_table, opt_duct_tape, opt_stats, opt_version, opt_MAX } ;

	FilterParams<Stream>::F filter_makers[opt_MAX] = {
		0, mk_sort_by_pos, mk_sort_by_name, mk_filter_by_length,
		mk_filter_by_score, mk_filter_by_mapq, mk_filter_by_hit, mk_delete_hit, mk_filter_qual, mk_subsample, mk_filter_multi, mk_edit_header, 0, 0,
		0, 0, mk_rmdup, 0, 0, 0, 0, 0,
		0, 0, 0, mk_duct_tape, 0, 0 } ;

	FilterParams<StreamBundle>::F merge_makers[opt_MAX] = {
		0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, mk_merge, mk_join,
		mk_mega_merge, mk_concat, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0 } ;

	FilterParams<Stream>::F output_makers[opt_MAX] = {
		0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, mk_output, mk_output_text, mk_output_sam, mk_output_glz, mk_output_3aln,
		mk_output_fasta, mk_output_fastq, mk_output_table, 0, mk_stats, 0 } ;

	G descriptions[opt_MAX] = {
		0, desc_sort_by_pos, desc_sort_by_name, desc_filter_by_length,
		desc_filter_by_score, desc_filter_by_mapq, desc_filter_by_hit, desc_delete_hit, desc_filter_qual, desc_subsample, desc_filter_multi, desc_edit_header, desc_merge, desc_join, 
		desc_mega_merge, desc_concat, desc_rmdup, desc_output, desc_output_text, desc_output_sam, desc_output_glz, desc_output_3aln,
		desc_output_fasta, desc_output_fastq, desc_output_table, desc_duct_tape, desc_stats, 0 } ;

	ParamBlock param( 7.5, 20.0, 0, 0, 0 ) ;
	int POPT_ARG_DFLT = POPT_ARG_FLOAT | POPT_ARGFLAG_SHOW_DEFAULT ;
	int POPT_ARG_DINT = POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT ;

	struct poptOption options[] = {
		{ "sort-pos",      's', POPT_ARG_INT,    0, opt_sort_pos,      "sort by alignment position [using <n MiB memory]", "n" },
		{ "sort-name",     'S', POPT_ARG_INT,    0, opt_sort_name,     "sort by read name [using <n MiB memory]", "n" },
		{ "filter-length", 'l', POPT_ARG_INT,    0, opt_filter_length, "filter for length of at least L", "L" },
		{ "filter-score",  'f', POPT_ARG_NONE,   0, opt_filter_score,  "filter for max score", 0 },
		{ "filter-mapq",    0 , POPT_ARG_INT,    0, opt_filter_mapq,   "remove alignments with MAPQ below Q", "Q" },
		{ "filter-hit",    'h', POPT_ARG_STRING, 0, opt_filter_hit,    "filter for hitting (in G) SEQ/anything", "SEQ" },
		{ "delete-hit",     0 , POPT_ARG_STRING, 0, opt_delete_hit,    "delete hits (in G/anywhere) to SEQ/anything", "SEQ" },
		{ "filter-qual",    0 , POPT_ARG_INT,    0, opt_filter_qual,   "delete bases with quality below Q", "Q" },
		{ "subsample",      0,  POPT_ARG_FLOAT,  0, opt_subsample,     "subsample a fraction F of the results", "F" },
		{ "multiplicity",   0 , POPT_ARG_INT,    0, opt_filter_multi,  "keep reads with multiplicity above N", "N" },
		{ "edit-header",    0 , POPT_ARG_STRING, 0, opt_edit_header,   "invoke editor ED on the stream's header", "ED" },
		{ "concat",        'c', POPT_ARG_NONE,   0, opt_concat,        "concatenate streams", 0 },
		{ "merge",         'm', POPT_ARG_NONE,   0, opt_merge,         "merge sorted streams", 0 },
		{ "join",          'j', POPT_ARG_NONE,   0, opt_join,          "join streams and retain best hits", 0 },
		{ "mega-merge",     0 , POPT_ARG_NONE,   0, opt_mega_merge,    "merge many streams, e.g. from grid jobs", 0 },
		{ "rmdup",         'd', POPT_ARG_NONE,   0, opt_rmdup,         "remove PCR duplicates", 0 },
		{ "output",        'o', POPT_ARG_STRING, 0, opt_output,        "write native stream to file FILE", "FILE" },
		{ "output-text",    0 , POPT_ARG_STRING, 0, opt_output_text,   "write protobuf text stream to FILE", "FILE" },
		{ "output-sam",     0 , POPT_ARG_STRING, 0, opt_output_sam,    "write alignments in sam format to FILE", "FILE" },
		{ "output-glz",     0 , POPT_ARG_STRING, 0, opt_output_glz,    "write contigs in GLZ format to FILE", "FILE" },
		{ "output-3aln",    0 , POPT_ARG_STRING, 0, opt_output_3aln,   "write contigs in 3ALN format to FILE", "FILE" },
		{ "output-fasta",   0 , POPT_ARG_STRING, 0, opt_output_fasta,  "write alignments(!) in fasta format to FILE", "FILE" },
		{ "output-fastq",   0 , POPT_ARG_STRING, 0, opt_output_fastq,  "write sequences(!) in fastq format to FILE", "FILE" },
		{ "output-table",   0 , POPT_ARG_STRING, 0, opt_output_table,  "write per-alignment stats to FILE", "FILE" },
		{ "duct-tape",      0 , POPT_ARG_STRING, 0, opt_duct_tape,     "mock-assemble while clamping Q-scores to Q", "Q" },
		{ "stats",          0,  POPT_ARG_STRING, 0, opt_stats,         "write simple statistics to FILE", "FILE" },

		{ "set-slope",      0 , POPT_ARG_DFLT,   &param.slope,      0, "set slope parameter to S", "S" },
		{ "set-intercept",  0 , POPT_ARG_DFLT,   &param.intercept,  0, "set length discount parameter to L", "L" },
		{ "set-context",    0 , POPT_ARG_DINT,   &param.context,    0, "set context parameter to C", "C" },
		{ "set-genome",     0 , POPT_ARG_STRING, &param.genome,     0, "set interesting genome parameter to G", "G" },
		{ "clear-genome",   0 , POPT_ARG_VAL,    &param.genome,     0, "clear interesting genome parameter", 0 },

		{ "quiet",         'q', POPT_ARG_VAL,    &console.loglevel, Console::error, "suppress most output", 0 },
		{ "verbose",       'v', POPT_ARG_VAL,    &console.loglevel, Console::info,  "produce more output", 0 },
		{ "debug",          0 , POPT_ARG_VAL,    &console.loglevel, Console::debug, "produce debugging output", 0 },
		{ "version",       'V', POPT_ARG_NONE,   0, opt_version,       "print version number and exit", 0 },
		POPT_AUTOHELP POPT_TABLEEND
	} ;

	FilterStack filters_initial ;
	FilterStack *filters_current = &filters_initial ;
	FilterParams< StreamBundle > merging_filter( ParamBlock(0,0,0,0,0), 0, mk_concat, desc_concat ) ;

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
			filters_current->push_back( FilterParams< Stream >(
						param, poptGetOptArg1( pc ), filter_makers[rc],
						descriptions[rc] ) );
		}
		else if( rc >= 0 && merge_makers[rc] )
		{
			// make sure we are still creating input filters
			if( filters_current != &filters_initial )
				throw "merge-like commands cannot not follow merge- or output-like commands" ;

			merging_filter = FilterParams< StreamBundle >(
				param, poptGetOptArg1( pc ), merge_makers[rc], descriptions[rc] ) ;

			// from now on we build output filters
			filters_current = &filters_terminal.back() ;
		}
		else if( rc >= 0 && output_makers[rc] )
		{
			// make sure we are creating output filters
			if( filters_current == &filters_initial )
				filters_current = &filters_terminal.back() ;

			// create filter
			filters_current->push_back( FilterParams< Stream >(
				param, poptGetOptArg( pc ), output_makers[rc], descriptions[rc] ) ) ;

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
	the_glob.gl_pathv = 0 ;
	the_glob.gl_pathc = 0 ;
	int glob_flag = GLOB_NOSORT ;

	while( const char* arg = poptGetArg( pc ) )
	{
		glob( arg, glob_flag, 0, &the_glob ) ;
		glob_flag |= GLOB_APPEND ;
	}

	// give a report of the filter stack (because it's so easy to get
	// wrong)
	console.output( Console::notice, "Input files:" ) ;
	if( the_glob.gl_pathc )
		for( char **arg = the_glob.gl_pathv ; arg != the_glob.gl_pathv + the_glob.gl_pathc ; ++arg )
			console.output( Console::notice, "  " + string(*arg) ) ;
	else console.output( Console::notice, "  <stdin>" ) ;

	if( !filters_initial.empty() ) 
	{
		console.output( Console::notice, "For every input file:" ) ;
		for( FilterStack::const_iterator i = filters_initial.begin() ; i != filters_initial.end() ; ++i )
		{
			stringstream s ;
			(i->describe)( s << "  ", *i ) ;
			console.output( Console::notice, s.str() ) ;
		}
	}
	{
		stringstream s ;
		merging_filter.describe( s, merging_filter ) ;
		console.output( Console::notice, s.str() ) ;
	}

	// last filter stack is not empty (== missing output filter) or only
	// one filter stack (== no output at all)?  --> add a writer for
	// stdout to last filter.  else remove the empty one
	if( !filters_terminal.back().empty() || filters_terminal.size() == 1 )
		filters_terminal.back().push_back( FilterParams< Stream >(
			param, 0, mk_output_text, desc_output_text ) ) ;
	else
		filters_terminal.pop_back() ;

	for( FilterStacks::const_iterator i = filters_terminal.begin() ; i != filters_terminal.end() ; ++i )
	{
		console.output( Console::notice, "Filter a copy of the result as follows:" ) ;
		for( FilterStack::const_iterator j = i->begin() ; j != i->end() ; ++j )
		{
			stringstream s ;
			(j->describe)( s << "  ", *j ) ;
			console.output( Console::notice, s.str() ) ;
		}
	}

	std::auto_ptr< StreamBundle > merging_stream( (merging_filter.maker)( merging_filter ) ) ;

	// iterate over glob results
	if( the_glob.gl_pathc )
	{
		vector< Compose* > cs ;
		for( char **arg = the_glob.gl_pathv ; arg != the_glob.gl_pathv + the_glob.gl_pathc ; ++arg )
		{
			Compose *c = new Compose ;
			c->add_stream( make_input_stream( *arg ) ) ; // new AnfoReader( *arg ) ) ;
			for( FilterStack::const_iterator i = filters_initial.begin() ; i != filters_initial.end() ; ++i )
				c->add_stream( (i->maker)( *i ) ) ;
			cs.push_back( c ) ;
		}
		for( vector< Compose* >::const_iterator i = cs.begin(), ie = cs.end() ; i != ie ; ++i )
			merging_stream->add_stream( *i ) ;
	}
	else merging_stream->add_stream( make_input_stream( /*new AnfoReader( 0,*/ dup( 0 ), "<stdin>" ) ) ;

	FanOut out ;
	for( FilterStacks::const_iterator i = filters_terminal.begin() ; i != filters_terminal.end() ; ++i )
	{
		auto_ptr< Compose > c( new Compose ) ;
		for( FilterStack::const_iterator j = i->begin() ; j != i->end() ; ++j )
			c->add_stream( (j->maker)( *j ) ) ;
		out.add_stream( c.release() ) ;
	}

	transfer( *merging_stream, out ) ;
	poptFreeContext( pc ) ;
	return 0 ;
}

