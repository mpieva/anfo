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
	vector<string> gs_ ;

	by_genome_coordinate( const string &g ) : gs_() { gs_.push_back( g ) ; }
	by_genome_coordinate( const vector<string> &gs ) : gs_(gs) {}

	bool operator() ( const Result *a, const Result *b ) {
		if( gs_.empty() ) return compare( *a, *b ) ;
		for( vector<string>::const_iterator i = gs_.begin(), e = gs_.end() ; i != e ; ++i )
		{
			if( compare( *a, *b, *i ) ) return true ;
			if( compare( *b, *a, *i ) ) return false ;
		}
		return false ;
	}

	bool compare( const Result &a, const Result &b )
	{
		const Hit *u = hit_to( a ), *v = hit_to( b ) ;
		if( u && !v ) return true ;
		if( !u ) return false ;

		if( u->genome_name() < v->genome_name() ) return true ;
		if( v->genome_name() < u->genome_name() ) return false ;
		if( u->sequence() < v->sequence() ) return true ;
		if( v->sequence() < u->sequence() ) return false ;
		if( u->start_pos() < v->start_pos() ) return true ;
		if( v->start_pos() < u->start_pos() ) return false ;
		return u->aln_length() < v->aln_length() ;
	}

	bool compare( const Result &a, const Result &b, const string& g )
	{
		const Hit *u = hit_to( a, g ), *v = hit_to( b, g ) ;
		if( u && !v ) return true ;
		if( !u ) return false ;

		if( u->sequence() < v->sequence() ) return true ;
		if( v->sequence() < u->sequence() ) return false ;
		if( u->start_pos() < v->start_pos() ) return true ;
		if( v->start_pos() < u->start_pos() ) return false ;
		return u->aln_length() < v->aln_length() ;
	}

	void tag_header( output::Header& h ) {
		h.clear_is_sorted_by_name() ;
		h.clear_is_sorted_by_coordinate() ;
		if( gs_.empty() ) {
			h.set_is_sorted_by_all_genomes( true ) ;
		}
		else
		{
			h.clear_is_sorted_by_all_genomes() ;
			for( vector<string>::const_iterator i = gs_.begin(), e = gs_.end() ; i != e ; ++i )
				h.add_is_sorted_by_coordinate( *i ) ;
		}
	}

    bool is_sorted( const output::Header& h ) {
		if( gs_.empty() ) return h.is_sorted_by_all_genomes() ;
		
		if( (int)gs_.size() != h.is_sorted_by_coordinate_size() ) return false ;
		return equal( gs_.begin(), gs_.end(), h.is_sorted_by_coordinate().begin() ) ;
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
//!
//! \note All input streams must necessarily be opened at the same time.
//!       In principle, we could merge smaller chunks into temporary
//!       files to avoid hitting the file descriptor limit, but there
//!       hasn't been a practical need to do that yet.
class MergeStream : public StreamBundle
{
	private:
		deque< Result > rs_ ;
		enum { unknown, by_name, by_coordinate } mode_ ;
		vector<string> gs_ ;

	public:
		MergeStream() : mode_( unknown ) {}
		virtual Header fetch_header() ;
		virtual Result fetch_result() ;

		//! \todo totally broken, need to think about how to keep the
		//!       necessary information.
#if HAVE_ELK_SCHEME_H
		virtual Object get_summary() const { return False ; }
#endif
} ;

//! \brief merges multiple streams by taking the best hit
//! If asked for a result, this class takes results from each of the
//! streams in turn, until one sequence has received an entry from each.
//! This sequence is then delivered.  Everything works fine, as long as
//! the sequence names come in in roughly the same order and every name
//! is contained in every input.  Else it still works, but eats memory.
//! This is best used to merge chunks of work done by independent
//! processes where the order of records hasn't been disturbed too much
//! (multithreading and the implied slight shuffle is fine).
class NearSortedJoin : public StreamBundle
{
	private:
		typedef map< string, pair< size_t, Result > > Buffer ;
		Buffer buffer_ ;
		size_t nread_, nwritten_, nstreams_ ;
		deque< StreamHolder >::iterator cur_input_ ;
		Chan progress_ ;

	public:
		NearSortedJoin() : nread_(0), nwritten_(0), nstreams_(0) {}

		virtual Header fetch_header() ;
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
		typedef deque< Holder< streams::Stream > > MergeableQueue ;
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
		Holder< MergeStream > final_stream_ ;

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

		void enqueue_stream( streams::StreamHolder, int = 0 ) ;
		void flush_scratch() ;

		virtual ~SortingStream()
		{
			for_each( scratch_space_.begin(), scratch_space_.end(), delete_ptr<Result>() ) ;
			--SortingStream__ninstances ;
		}

	public:
		SortingStream( unsigned as = 256*1024*1024, unsigned qs = 256, Comp comp = Comp() )
			: final_stream_( new MergeStream ), total_scratch_size_(0), max_que_size_( qs ), max_arr_size_( as ), comp_( comp )
		{ foot_.set_exit_code( 0 ) ; ++SortingStream__ninstances ; }

		virtual void put_header( const Header& h ) { Stream::put_header( h ) ; comp_.tag_header( hdr_ ) ; }
		virtual void put_footer( const Footer& ) ;
		virtual void put_result( const Result& r ) {
			scratch_space_.push_back( new Result( r ) ) ;
			total_scratch_size_ += scratch_space_.back()->SpaceUsed() ;
			if( total_scratch_size_ >= max_arr_size_ ) flush_scratch() ;
		}

		virtual Result fetch_result()
		{
			Result r = final_stream_->fetch_result() ;
			state_ = final_stream_->get_state() ;
			return r ;
		}
		virtual Footer fetch_footer() { merge_sensibly( foot_, final_stream_->fetch_footer() ) ; return foot_ ; }
#if HAVE_ELK_SCHEME_H
		virtual Object get_summary() const { return final_stream_->get_summary() ; }
#endif
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
	google::protobuf::io::FileInputStream* is = new google::protobuf::io::FileInputStream( fd ) ;
	is->SetCloseOnDelete( true ) ;
	enqueue_stream( new UniversalReader( tempname, is ), 1 ) ;

	for_each( scratch_space_.begin(), scratch_space_.end(), delete_ptr<Result>() ) ;
	scratch_space_.clear() ;
	total_scratch_size_ = 0 ;
}

template < typename Comp > void SortingStream<Comp>::enqueue_stream( streams::StreamHolder s, int level ) 
{
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
			google::protobuf::io::FileInputStream* is = new google::protobuf::io::FileInputStream( fd ) ;
			is->SetCloseOnDelete( true ) ;
			enqueue_stream( new UniversalReader( fname, is ), max_bin + 1 ) ;
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
	// temporary storage.  (Also, if the temporay space is empty, we
	// *always* add the ContainerStream, otherwise we get strange
	// effects if the output turns out to be empty.)

	if( scratch_space_.begin() != scratch_space_.end() && SortingStream__ninstances > 1 )
		flush_scratch() ; 
	else {
		console.output( Console::notice, "SortingStream: final sort" ) ;
		sort_scratch() ;
		final_stream_->add_stream( new ContainerStream< deque< Result* >::const_iterator >(
					hdr_, scratch_space_.begin(), scratch_space_.end(), foot_ ) ) ;
	}

	// add any streams that have piled up
	for( MergeableQueues::const_iterator i = mergeable_queues_.begin() ; i != mergeable_queues_.end() ; ++i )
		for( MergeableQueue::const_iterator j = i->second.begin() ; j != i->second.end() ; ++j )
			final_stream_->add_stream( *j ) ;
	mergeable_queues_.clear() ;
	console.output( Console::notice, "SortingStream: merging everything to output" ) ;
	
	final_stream_->fetch_header() ;
	state_ = final_stream_->get_state() ;
}


class RepairHeaderStream : public Stream
{
	private:
		string editor_ ;

	public:
		RepairHeaderStream( const string &e ) : editor_( e ) {}
		virtual void put_header( const Header& ) ;
} ;

class FanOut : public StreamBundle
{
	public:
		virtual void put_header( const Header& ) ;
		virtual void put_result( const Result& ) ;
		virtual void put_footer( const Footer& ) ;
		virtual Footer fetch_footer() ;
} ;

class Compose : public StreamBundle
{
	public:
		virtual state get_state() ;

		virtual void put_header( const Header&   ) ;
		virtual void put_result( const Result& r ) { streams_.front()->put_result( r ) ; get_state() ; }
		virtual void put_footer( const Footer& f ) { streams_.front()->put_footer( f ) ; get_state() ; }

		virtual Header fetch_header() ;
		virtual Result fetch_result() { return streams_.back()->fetch_result() ; }
		virtual Footer fetch_footer() { return streams_.back()->fetch_footer() ; }
} ;

class StatStream : public Stream
{
	private:
		string fn_ ;
		string name_ ;

		unsigned total_, mapped_, mapped_u_, different_ ;
		uint64_t bases_, bases_gc_, bases_m_, bases_gc_m_ ;
		uint64_t bases_squared_, bases_m_squared_ ; 

		void printout( ostream&, bool ) ;

	public:
		StatStream( const string& fn )
			: fn_(fn), total_(0), mapped_(0), mapped_u_(0), different_(0)
		    , bases_(0), bases_gc_(0), bases_m_(0), bases_gc_m_(0)
		    , bases_squared_(0), bases_m_squared_(0) {}

		virtual void put_result( const Result& ) ;
		virtual void put_footer( const Footer& ) ;
#if HAVE_ELK_SCHEME_H
		virtual Object get_summary() const ;
#endif
} ;

//! \brief calculates divergence
//! Divergence is calculated by the Green Triangulation method.  Two
//! genomes are needed (primary and secondary) and alignments to both of
//! them.  Differences are counted (all equal, either sequence
//! different, all different), then error corrected (math stolen from
//! dropin_AHA.pl) and turned into divergence of the last common
//! ancestor, expressed as fraction of total divergence between primary
//! and secondary genome.

class DivergenceStream : public Stream
{
	private:
		string primary_genome_, secondary_genome_ ;
		int chop_ ;
        bool ancient_ ;
		int64_t b1, b2, b3, b4, b5 ;

	public:
		DivergenceStream( const string& primary, const string& secondary, int chop )
			: primary_genome_( primary ), secondary_genome_( secondary ), chop_( chop )
			, b1(0), b2(0), b3(0), b4(0), b5(0) {}

		virtual void put_header( const Header& ) ;
		virtual void put_result( const Result& ) ;
#if HAVE_ELK_SCHEME_H
		virtual Object get_summary() const ;
#endif
} ;

class MismatchStats : public Stream
{
	private:
		int mat_[4][4] ;

	public:
		MismatchStats() { memset( mat_, 0, sizeof( mat_ ) ) ; }
		virtual void put_result( const Result& ) ;
#if HAVE_ELK_SCHEME_H
		virtual Object get_summary() const ;
#endif
} ;

class RegionFilter : public HitFilter
{
	private:
		typedef std::map< unsigned, unsigned > Regions ;		// end(!) & start
		typedef std::map< std::string, Regions > Regions2 ;		// chromosome
		typedef std::map< std::string, Regions2 > Regions3 ; 	// filename

		static Regions3 all_regions ;

		Regions2 *my_regions ;

	public:
		RegionFilter( const pair< istream*, string >& ) ;

		//! \brief looks for a region overlapping the current alignment
		//! This will only work correctly for non-overlapping
		//! annotations.  Deal with it.
		bool inside( const Hit& h )
		{
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
		InsideRegion( const pair< istream*, string > &p ) : RegionFilter( p ) {}
		virtual bool keep( const Hit& h ) { return inside( h ) ; }
} ;
class OutsideRegion : public RegionFilter
{
	public:
		OutsideRegion( const pair< istream*, string > &p ) : RegionFilter( p ) {}
		virtual bool keep( const Hit& h ) { return !inside( h ) ; }
} ;

class Sanitizer : public Filter
{
	public:
		virtual void put_header( const Header& ) ;
		virtual bool xform( Result& ) ;
} ;

} // namespace

#endif
