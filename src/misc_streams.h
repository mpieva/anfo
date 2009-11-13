#ifndef INCLUDED_MISC_STREAMS_H
#define INCLUDED_MISC_STREAMS_H

#include "stream.h"

#include <deque>

namespace streams {

	using namespace std ;

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

    bool is_sorted( const output::Header& h ) {
        return h.has_is_sorted_by_coordinate() 
            && (!g_ || h.is_sorted_by_coordinate() == g_) ;
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
    bool is_sorted( const output::Header& h ) {
        return h.is_sorted_by_name() ;
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

		//! \brief checks if enough inputs are known.
		bool enough_inputs() const ;

		//! \brief reports final results for one sequence.
		virtual Result fetch_result() ;
} ;


//! \brief presents a container as an input stream
//! The container must be of pointers to \c Result, the stream acts as
//! input stream.  A suitable header and footer are supplied at
//! construction time.
template< typename I > class ContainerStream : public Stream
{
	private:
		I cur_, end_ ;
		unsigned total_, done_ ;
		Chan chan_ ;

	public:
		ContainerStream( const Header &hdr, I begin, I end, const Footer &foot ) 
			: cur_( begin ), end_( end ), total_( std::distance( begin, end ) ), done_(0)
		{
			hdr_ = hdr ;
			foot_ = foot ;
			state_ = cur_ == end_ ? end_of_stream : have_output ;
		}

		virtual Result fetch_result() {
			if( ++done_ % 1024 == 0 )
			{
				stringstream s ;
				s << "(mem) " << done_ << "/" << total_
					<< " (" << (int)(100*done_/total_) << "%)" ;
				chan_( Console::info, s.str() ) ;
			}

			Result r = **cur_ ;
			++cur_ ;
			if( cur_ == end_ ) state_ = end_of_stream ;
			return r ;
		}
} ;

extern unsigned SortingStream__ninstances ;

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
	enqueue_stream( make_input_stream( fd, tempname.c_str() ), 1 ) ;

	for_each( scratch_space_.begin(), scratch_space_.end(), delete_ptr<Result>() ) ;
	scratch_space_.clear() ;
	total_scratch_size_ = 0 ;
}

template < typename Comp > void SortingStream<Comp>::enqueue_stream( streams::Stream* s, int level ) 
{
	Header h = s->fetch_header() ;
	assert( comp_.is_sorted( h ) ) ;
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
			enqueue_stream( make_input_stream( fd, fname.c_str() ), max_bin + 1 ) ;
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

class RepairHeaderStream : public Stream
{
	private:
		const char* editor_ ;

	public:
		RepairHeaderStream( const char* e ) : editor_( e ) {}
		virtual ~RepairHeaderStream() {}
		virtual void put_header( const Header& ) ;
} ;

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

class IgnoreHit : public Filter
{
	private:
		const char* g_ ;
		const char* s_ ;

	public:
		IgnoreHit( const char* g, const char* s ) : g_(g), s_(s) {}
		virtual ~IgnoreHit() {}
		virtual void put_header( const Header& h ) {
			Filter::put_header( h ) ;
			hdr_.clear_is_sorted_by_coordinate() ;
		}

		virtual bool xform( Result& r ) 
		{
			int j = 0 ;
			for( int i = 0 ; i != r.hit_size() ; ++i )
			{
				if( (g_ && *g_ && r.hit(i).genome_name() != g_)
						|| (s_ && *s_ && r.hit(i).sequence() != s_) )
				{
					if( i != j ) swap( *r.mutable_hit(j), *r.mutable_hit(i) ) ;
					++j ;
				}
			}
			while( j != r.hit_size() ) r.mutable_hit()->RemoveLast() ;
			return true ;
		}
} ;

class RegionFilter : public Filter
{
	private:
		typedef std::map< unsigned, unsigned > Regions ;		// end(!) & start
		typedef std::map< std::string, Regions > Regions2 ;		// chromosome
		typedef std::map< std::string, Regions2 > Regions3 ; 	// filename

		static Regions3 all_regions ;

		Regions2 *my_regions ;

	public:
		RegionFilter( const std::string &fn ) ;

		//! \brief looks for a region overlapping the current alignment
		//! This will only work correctly for non-overlapping
		//! annotations.  Deal with it.
		bool inside( const Result& res )
		{
			const Hit &h = hit_to( res, 0 ) ;
			unsigned x = h.start_pos() ;
			const Regions &r = (*my_regions)[ h.sequence() ] ;
			Regions::const_iterator i = r.lower_bound( x ) ;
			// we now got the leftmost region whose end is to the right
			// of our start.  if no such thing exists, nothing overlaps.
			if( i == r.end() ) return false ;
			// now check if what we got actually overlaps.  it does if
			// it starts before our alignment ends.
			return i->second <= x + abs(h.aln_length()) && i->first >= x ;
		}
} ;

class InsideRegion : public RegionFilter
{
	public:
		InsideRegion( const std::string& fn ) : RegionFilter( fn ) {}
		bool xform( Result& r ) { return inside( r ) ; }
} ;
class OutsideRegion : public RegionFilter
{
	public:
		OutsideRegion( const std::string& fn ) : RegionFilter( fn ) {}
		bool xform( Result& r ) { return !inside( r ) ; }
} ;


} // namespace

#endif
