#ifndef INCLUDED_MISC_STREAMS_H
#define INCLUDED_MISC_STREAMS_H

#include "stream.h"

#include <deque>
#include <map>

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
		const Hit *u = best_hit( a ), *v = best_hit( b ) ;
		if( u && !v ) return true ;
		if( !u ) return false ;

		if( u->genome_name() < v->genome_name() ) return true ;
		if( v->genome_name() < u->genome_name() ) return false ;
		if( u->sequence() < v->sequence() ) return true ;
		if( v->sequence() < u->sequence() ) return false ;
		if( u->start_pos() < v->start_pos() ) return true ;
		if( v->start_pos() < u->start_pos() ) return false ;
		return lexicographical_compare(
				u->cigar().begin(), u->cigar().end(), 
				v->cigar().begin(), v->cigar().end() ) ;
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
		return lexicographical_compare(
				u->cigar().begin(), u->cigar().end(), 
				v->cigar().begin(), v->cigar().end() ) ;
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
//! \note If records with the same name appear after one another, they
//!       get "merged sensibly".  This is the intended behaviour when
//!       sorting by name, otherwise it does no harm.
//! \todo Broken at the moment.
class MergeStream : public InputStream
{
	private:
		vector< pair< Result*, StreamHolder > > heap_ ;
		vector< StreamHolder > streams_ ;

		enum { unknown, by_name, by_coordinate } mode_ ;
		vector<string> gs_ ;
		auto_ptr< Footer > ftr_ ;

	public:
		MergeStream( const vector< StreamHolder >& ss ) :
			streams_( ss ), mode_( unknown ), ftr_( new Footer ) {}

		virtual state priv_get_state() {
			return !streams_.empty() ? have_header :
				   !heap_.empty() ? have_output :
				   end_of_stream ;
		}

		virtual auto_ptr< Header > priv_fetch_header() ;
		virtual auto_ptr< Result > priv_fetch_result() ;
		virtual auto_ptr< Footer > priv_fetch_footer() { return ftr_ ; }

		virtual string type_name() const { return "MergeStream" ; }

		bool compare( const Result* a, const Result* b ) const
		{
			return 
				mode_ == by_name ? by_seqid()( a, b ) :
				by_genome_coordinate( gs_ )( a, b ) ;
		}
} ;

//! \brief merges multiple streams by taking the best hit
//! If asked for a result, this class takes results from each of the
//! streams in turn, until one sequence has received an entry from each.
//! This sequence is then delivered.  Everything works fine, as long as
//! the sequence names come in in roughly the same order and every name
//! is contained in every input.  Else it still works, but eats memory.
//! You get a warning if that happens, and probably a crash later on.
//! This is best used to merge chunks of work done by independent
//! processes where the order of records hasn't been disturbed too much
//! (multithreading and the implied slight shuffle is fine).
class NearSortedJoin : public InputStream
{
	private:
		vector< StreamHolder > streams_ ;
		typedef map< string, pair< size_t, Result* > > Buffer ;
		Buffer buffer_ ;
		size_t nread_, nwritten_, nstreams_ ;
		vector< StreamHolder >::iterator cur_input_ ;
		Chan progress_ ;
		auto_ptr< Footer > ftr_ ;
		bool warned_ ;

	public:
		NearSortedJoin( const vector< StreamHolder >& ss ) :
			streams_( ss ), nread_(0), nwritten_(0), nstreams_( streams_.size() ),
			cur_input_( streams_.end() ), ftr_( new Footer ), warned_( false ) {}

		virtual state get_state() {
			return streams_.empty() && buffer_.empty() ? end_of_stream :
				cur_input_ == streams_.end() ? have_header :
				have_output ;
		}

		virtual auto_ptr< Header > priv_fetch_header() ;
		virtual auto_ptr< Result > priv_fetch_result() ;
		virtual auto_ptr< Footer > priv_fetch_footer() { return ftr_ ; }

		virtual string type_name() const { return "NearSortedJoin" ; }
} ;


//! \brief presents a container as an input stream
//! The container must be of pointers to \c Result, the stream acts as
//! input stream and takes ownership(!) of the contained pointers.  A
//! suitable header and footer are supplied at construction time.
template< typename I > class ContainerStream : public InputStream
{
	private:
		I cur_, end_ ;
		unsigned total_, done_ ;
		Chan chan_ ;

		auto_ptr< Header > hdr_ ;
		auto_ptr< Footer > ftr_ ;

	public:
		ContainerStream( const Header &hdr, I begin, I end, const Footer &foot ) :
			cur_( begin ), end_( end ), total_( std::distance( begin, end ) ), done_(0),
			hdr_( new Header( hdr ) ), ftr_( new Footer( foot ) )
		{}

		virtual state priv_get_state() {
			return hdr_.get() ? have_header :
				cur_ == end_ ? end_of_stream : 
				have_output ;
		}

		virtual auto_ptr< Header > priv_fetch_header() { return hdr_ ; }
		virtual auto_ptr< Footer > priv_fetch_footer() { return ftr_ ; }

		virtual auto_ptr< Result > priv_fetch_result() {
			if( ++done_ % 1024 == 0 )
			{
				stringstream s ;
				s << "(mem) " << done_ << "/" << total_
					<< " (" << (int)(100*done_/total_) << "%)" ;
				chan_( Console::info, s.str() ) ;
			}

			auto_ptr< Result > r( *cur_ ) ;
			*cur_ = 0 ;
			++cur_ ;
			return r ;
		}

		virtual string type_name() const { return "ContainerStream" ; }
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

template <class Comp> class SortStream : public Stream
{
	private:
		typedef deque< StreamHolder > MergeableQueue ;
		typedef map< unsigned, MergeableQueue > MergeableQueues ;
		typedef deque< Result* > ScratchSpace ;

		//! Storage to perform quicksort in.
		ScratchSpace scratch_space_ ;

		//! We keep multiple queues of streams ordered and separated by
		//! the number of times they have been merged.  This way we
		//! don't merge the stuff over and over, and instead tend to do
		//! a merge that gives the greatest reduction in number of open
		//! files with the least amount of IO.  The key is the number of
		//! times a file has already been merged with others, the value
		//! is just a set of streams.
		MergeableQueues mergeable_queues_ ;
		StreamHolder final_stream_ ;

		uint64_t total_scratch_size_ ;

		unsigned max_que_size_ ;
        unsigned num_open_files_ ;
		unsigned max_arr_size_ ;
		Comp comp_ ;

		auto_ptr< Header > hdr_ ;

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

		void enqueue_stream( StreamHolder, int = 0 ) ;
		void flush_scratch() ;

		virtual ~SortStream()
		{
			for_each( scratch_space_.begin(), scratch_space_.end(), delete_ptr<Result>() ) ;
			--SortingStream__ninstances ;
		}

	public:
		SortStream( unsigned as = 256, unsigned qs = 256, Comp comp = Comp() ) :
            final_stream_(), total_scratch_size_(0), max_que_size_(qs),
            num_open_files_(0), max_arr_size_( as), comp_( comp )
		{ ++SortingStream__ninstances ; }

		virtual state priv_get_state() 
		{
			return final_stream_ ? final_stream_->get_state() :
				hdr_.get() ? need_input : need_header ;
		}

		virtual void priv_put_header( auto_ptr< Header > h ) { hdr_ = h ; comp_.tag_header( *hdr_ ) ; }
		virtual void priv_put_footer( auto_ptr< Footer > ) ;

		virtual void priv_put_result( auto_ptr< Result > r ) {
			scratch_space_.push_back( r.release() ) ;
			total_scratch_size_ += scratch_space_.back()->SpaceUsed() ;
			if( (total_scratch_size_ >> 20) >= max_arr_size_ ) flush_scratch() ;
		}

		virtual auto_ptr< Header > priv_fetch_header() { return final_stream_->fetch_header() ; }
		virtual auto_ptr< Result > priv_fetch_result() { return final_stream_->fetch_result() ; }
		virtual auto_ptr< Footer > priv_fetch_footer() { return final_stream_->fetch_footer() ; }
		virtual Object get_summary() const { return final_stream_->get_summary() ; }
} ;

template < typename Comp > void SortStream<Comp>::flush_scratch()
{
	sort_scratch() ;
	string tempname ;
	int fd = mktempfile( &tempname ) ;
	{
		ContainerStream< deque< Result* >::const_iterator >
			sa( hdr_, scratch_space_.begin(), scratch_space_.end(), new Footer ) ;
		ChunkedWriter out( fd, 25, tempname.c_str() ) ;
		console.output( Console::notice, "SortStream: Writing to tempfile " + tempname ) ;
		transfer( sa, out ) ;
	}
	throw_errno_if_minus1( lseek( fd, 0, SEEK_SET ), "seeking in ", tempname.c_str() ) ;
	google::protobuf::io::FileInputStream* is = new google::protobuf::io::FileInputStream( fd ) ;
	is->SetCloseOnDelete( true ) ;
	enqueue_stream( new UniversalReader( tempname, is ), 0 ) ;

	for_each( scratch_space_.begin(), scratch_space_.end(), delete_ptr<Result>() ) ;
	scratch_space_.clear() ;
	total_scratch_size_ = 0 ;
}

template < typename Comp > void SortStream<Comp>::enqueue_stream( StreamHolder s, int level ) 
{
	mergeable_queues_[ level ].push_back( s ) ;

	if( ++num_open_files_ > max_que_size_ ) {
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
		// (only important in a weird corner case, in which we will
		// simply use one more file descriptor then requested)
		if( total_inputs > 2 ) {
			string fname ;
			int fd = mktempfile( &fname ) ;
			stringstream s ;
			s << "SortStream: Merging bins 0.." << max_bin << " to tempfile " << fname ;
			console.output( Console::notice, s.str() ) ;
			{
				vector< StreamHolder > ss ;
				for( MergeableQueues::iterator i = mergeable_queues_.begin() ; i->first <= max_bin ; ++i ) 
				{
                    num_open_files_ -= i->second.size() ;
					for( size_t j = 0 ; j != i->second.size() ; ++j )
						ss.push_back( i->second[j] ) ;
					i->second.clear() ;
				}
				MergeStream ms( ss ) ;
				ChunkedWriter out( fd, 25, fname.c_str() ) ;
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
template < typename Comp > void SortStream<Comp>::priv_put_footer( auto_ptr< Footer > f ) 
{
	// We have to be careful about buffering; if more than one
	// SortStream is active, we could run out of RAM.  Therefore, if
	// we're alone, we sort and add a a stream.  Else we flush to
	// temporary storage.  (Also, even if the temporay space is empty,
	// we *always* add the ContainerStream, otherwise we get strange
	// effects if the output turns out to be empty.)
	vector< StreamHolder > ss ;
	if( scratch_space_.begin() != scratch_space_.end() && SortingStream__ninstances > 1 )
		flush_scratch() ; 
	else {
		console.output( Console::notice, "SortingStream: final sort" ) ;
		sort_scratch() ;
	}
	ss.push_back( new ContainerStream< deque< Result* >::const_iterator >(
				*hdr_, scratch_space_.begin(), scratch_space_.end(), *f ) ) ;

	// add any streams that have piled up
	for( MergeableQueues::const_iterator i = mergeable_queues_.begin() ; i != mergeable_queues_.end() ; ++i )
		for( MergeableQueue::const_iterator j = i->second.begin() ; j != i->second.end() ; ++j )
			ss.push_back( *j ) ;
	mergeable_queues_.clear() ;

    StreamHolder m( new MergeStream( ss ) ) ;
    if( SortingStream__ninstances == 1 )
    {
        // we're alone.  delegate directly to MergeStream
        console.output( Console::notice, "SortStream: merging everything to output" ) ;
        final_stream_ = m ;
    }
    else
    {
        // not alone.  We'll conserve file handles by merging everything
        // into a temporary file and opening that.
        string fname ;
        int fd = mktempfile( &fname ) ;
        stringstream s ;
        s << "SortStream: Merging everything to tempfile " << fname ;
        console.output( Console::notice, s.str() ) ;
        ChunkedWriter out( fd, 25, fname.c_str() ) ;
        transfer( *m, out ) ;
        num_open_files_ = 1 ;
        throw_errno_if_minus1( lseek( fd, 0, SEEK_SET ), "seeking in ", fname.c_str() ) ;
        google::protobuf::io::FileInputStream* is = new google::protobuf::io::FileInputStream( fd ) ;
        is->SetCloseOnDelete( true ) ;
        final_stream_ = new UniversalReader( fname, is ) ; 
    }
}


class RepairHeaderStream : public Stream
{
	private:
		string editor_ ;

	public:
		RepairHeaderStream( const string &e ) : editor_( e ) {}
		virtual void put_header( const Header& ) ;
} ;

class FanOut : public OutputStream
{
	private:
		vector< StreamHolder > streams_ ;
		typedef vector< StreamHolder >::iterator iter ;

	public:
		FanOut( const vector< StreamHolder >& ss ) : streams_( ss ) {}

		virtual state priv_get_state() ;

		virtual void priv_put_header( auto_ptr< Header > ) ;
		virtual void priv_put_result( auto_ptr< Result > ) ;
		virtual void priv_put_footer( auto_ptr< Footer > ) ;

		virtual string type_name() const { return "FanOut" ; }
} ;

class Compose : public Stream
{
	private:
		vector< StreamHolder > streams_ ;
		typedef vector< StreamHolder >::const_reverse_iterator criter ;

	public:
		Compose( const vector< StreamHolder >& ss ) : streams_( ss ) {}
		virtual state priv_get_state() ;

		virtual void priv_put_header( auto_ptr< Header > h ) { streams_.front()->put_header( h ) ; }
		virtual void priv_put_result( auto_ptr< Result > r ) { streams_.front()->put_result( r ) ; }
		virtual void priv_put_footer( auto_ptr< Footer > f ) { streams_.front()->put_footer( f ) ; }

		virtual auto_ptr< Header > priv_fetch_header() { return streams_.back()->fetch_header() ; }
		virtual auto_ptr< Result > priv_fetch_result() { return streams_.back()->fetch_result() ; }
		virtual auto_ptr< Footer > priv_fetch_footer() { return streams_.back()->fetch_footer() ; }

		virtual string type_name() const { return "Compose" ; }
} ;

class StatStream : public OutputStream
{
	private:
		unsigned total_, mapped_, mapped_u_, different_ ;
		uint64_t bases_, bases_gc_, bases_m_, bases_gc_m_ ;
		uint64_t bases_squared_, bases_m_squared_ ; 

	public:
		StatStream()
			: total_(0), mapped_(0), mapped_u_(0), different_(0)
		    , bases_(0), bases_gc_(0), bases_m_(0), bases_gc_m_(0)
		    , bases_squared_(0), bases_m_squared_(0) {}

		virtual state priv_get_state() { return need_input ; }
		virtual void priv_put_result( auto_ptr< Result > ) ;
		virtual Object get_summary() const ;
} ;

//! \brief calculates divergence
//! Divergence is calculated by the Green Triangulation method.  Two
//! genomes are needed (primary and secondary) and alignments to both of
//! them.  Differences are counted (all equal, either sequence
//! different, all different), then error corrected (math stolen from
//! dropin_AHA.pl) and turned into divergence of the last common
//! ancestor, expressed as fraction of total divergence between primary
//! and secondary genome.

class DivergenceStream : public OutputStream
{
	private:
		string primary_genome_, secondary_genome_ ;
		int chop_ ;
        bool ancient_, saw_header_ ;
		int64_t b1, b2, b3, b4, b5 ;

	public:
		DivergenceStream( const string& primary, const string& secondary, int chop )
			: primary_genome_( primary ), secondary_genome_( secondary ), chop_( chop )
			, saw_header_(false), b1(0), b2(0), b3(0), b4(0), b5(0) {}

		virtual state priv_get_state() { return saw_header_ ? need_header : need_input ; }
		virtual void priv_put_header( auto_ptr< Header > ) ;
		virtual void priv_put_result( auto_ptr< Result > ) ;
		virtual Object get_summary() const ;
} ;

class MismatchStats : public OutputStream
{
	private:
		int mat_[4][4] ;

	public:
		MismatchStats() { memset( mat_, 0, sizeof( mat_ ) ) ; }
		virtual state priv_get_state() { return need_input ; }
		virtual void priv_put_result( auto_ptr< Result > ) ;
		virtual Object get_summary() const ;
} ;

//! \brief checks for hits to homologous regions
//! This filter reads a UCSC Chain file and stores the "homologous"
//! ranges.  A record passes the filter iff it has hits to both genomes
//! that make up the chain and both hits are inside the two regions of
//! the same chain.  This should be equivalent to the "traditional"
//! sequence of two liftovers and check for agreement.
//!
//! \todo I swear, one of those days I'll implement a symbol table for
//!       those repeated chromosome names.
class AgreesWithChain : public Filter
{
	private:
		string left_genome_, right_genome_ ;

		struct Entry ;
		typedef map< unsigned, Entry > Chains ;	// =^= left_start
		typedef map< string, Chains > Map1 ;	// =^= left_chr

		// Chains have a hierarchical structure: below any chain, there
		// can be a collection of more.  We have one such top-level
		// collection per chromosome.
		struct Entry {
			unsigned left_end ;
			unsigned right_start ;
			unsigned right_end : 31 ;
			unsigned strand : 1 ;
			string right_chr ;
			Chains nested ;
		} ;

		Map1 map_ ;

		static Chains::iterator find_any_most_specific_overlap( 
				unsigned start, unsigned end, Chains *chains ) ;
	public:
		AgreesWithChain( const string& l, const string& r, const pair<istream*,string>& s ) ;
		virtual bool xform( Result& ) ;
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
