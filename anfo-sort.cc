#include "config.h"

#include "compress_stream.h"
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
//! files is usefult, but we take that from the environment.

//! \brief compares hits by smallest genome coordinate
//! Comparison is first done lexically on the subject name, then on the
//! smallest coordinate that's part of the alignment, then on the
//! alignment length.
struct by_genome_coordinate {
	bool operator() ( const Result *a, const Result *b ) {
		if( a->has_best_to_genome() && !b->has_best_to_genome() ) return true ;
		if( !a->has_best_to_genome() ) return false ;
		const Hit& u = a->best_to_genome(), v = b->best_to_genome() ;
		if( u.genome() < v.genome() ) return true ;
		if( v.genome() < u.genome() ) return false ;
		if( u.sequence() < v.sequence() ) return true ;
		if( v.sequence() < u.sequence() ) return false ;
		if( u.start_pos() + max(u.aln_length(),0) < v.start_pos() + max(v.aln_length(),0) ) return true ;
		if( v.start_pos() + max(v.aln_length(),0) < u.start_pos() + max(u.aln_length(),0) ) return false ;
		return u.aln_length() < v.aln_length() ;
	}
} ;

//! \brief stream of result messages
class Stream
{
	public:
		virtual ~Stream() {}
		virtual output::Header read_header() = 0 ;
		virtual output::Result read_result() = 0 ;
		virtual output::Footer read_footer() = 0 ;
} ;

class FileStream : public Stream
{
	private:
		AnfoFile file_ ;
		static int num_files_ ; // tracked to avoid bumping into the file descriptor limit

	public:
		FileStream( const std::string& name, bool unlink_on_delete = false )
			: file_( name, unlink_on_delete )
		{
			clog << "Reading from " << name << endl ;
			++num_files_ ;
		}
		virtual ~FileStream() { --num_files_ ; }
		virtual output::Header read_header() { return file_.read_header() ; }
		virtual output::Result read_result() { return file_.read_result() ; }
		virtual output::Footer read_footer() { return file_.read_footer() ; }

		static int num_open_files() { return num_files_ ; }
} ;

int FileStream::num_files_ = 0 ;

template< typename T > struct delete_ptr { void operator()( T* p ) const { delete p ; } } ;

class MergeStream : public Stream
{
	private:
		std::vector< Stream* > streams_ ;
		deque< Result > rs_ ;
		output::Header hdr_ ;
		output::Footer foot_ ;

	public:
		MergeStream() {}
		virtual ~MergeStream() { std::for_each( streams_.begin(), streams_.end(), delete_ptr<Stream>() ) ; }

		void add_stream( const Header& hdr, std::auto_ptr< Stream > s )
		{
			merge_sensibly( hdr_, hdr ) ;
			rs_.push_back( s->read_result() ) ;
			if( rs_.back().has_seqid() )
			{
				streams_.push_back( s.release() ) ;
			}
			else
			{
				rs_.pop_back() ;
				merge_sensibly( foot_, s->read_footer() ) ;
			}
		}

		virtual output::Header read_header() { return hdr_ ; }
		virtual output::Result read_result() ;
		virtual output::Footer read_footer() { return foot_ ; }
} ;

Result MergeStream::read_result() 
{
	for( size_t i = 0 ; i != streams_.size() ; ++i )
	{
		Result r = streams_[i]->read_result() ;
		if( r.has_seqid() ) rs_.push_back( r ) ;
	}

	if( !rs_.size() ) return Result() ;
	int min_idx = 0 ;
	for( size_t i = 1 ; i != rs_.size() ; ++i ) 
		if( by_genome_coordinate()( &rs_[i], &rs_[ min_idx ] ) )
			min_idx = i ;

	Result &rm = rs_[ min_idx ] ;
	Stream *sm = streams_[ min_idx ] ;

	Result r = rm ;
	rm = sm->read_result() ;
	if( !rm.has_seqid() ) 
	{
		merge_sensibly( foot_, sm->read_footer() ) ;
		delete sm ;
		streams_.erase( streams_.begin() + min_idx ) ;
		rs_.erase( rs_.begin() + min_idx ) ;
	}
	return r ;
}

int write_stream_to_file( int fd, Stream &s, bool expensive = false )
{
	FileOutputStream fos( fd ) ;
	std::auto_ptr< ZeroCopyOutputStream > zos(
			expensive ? compress_fast( &fos ) : compress_small( &fos ) ) ;
	CodedOutputStream o( zos.get() ) ;
	
	o.WriteRaw( "ANFO", 4 ) ;
	write_delimited_message( o, 1, s.read_header() ) ;
	for( Result r ; ( r = s.read_result() ).has_seqid() ; )
		write_delimited_message( o, 2, r ) ;
	Footer f = s.read_footer() ;
	write_delimited_message( o, 3, f ) ;
	return f.exit_code() ;
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
		std::vector< Stream* > streams_ ;
		output::Header hdr_ ;

		typedef map< string, pair< size_t, output::Result > > Buffer ;
		Buffer buffer_ ;
		size_t cur_input_, nread_, nwritten_ ;

	public:
		BestHitStream() : streams_(), cur_input_(0), nread_(0), nwritten_(0) {}
		virtual ~BestHitStream() { std::for_each( streams_.begin(), streams_.end(), delete_ptr<Stream>() ) ; }

		//! \brief reads a stream's header and adds the stream as input
		void add_stream( Header hdr, std::auto_ptr< Stream > s ) {
			merge_sensibly( hdr_, hdr ) ;
			streams_.push_back( s.release() ) ;
		}

		//! \brief checks if enough inputs are known
		//! Checks if at least as many inputs are known as requested.
		//! If that number isn't known, it assumed only one index was
		//! used and takes the number of slices declared for it.  
		bool enough_inputs( int n ) const { return streams_.size() >= n ; }
		bool enough_inputs() const
		{
			if( streams_.empty() ) return false ;
			if(  hdr_.config().policy_size() == 0 ) return true ;
			if(  hdr_.config().policy(0).use_compact_index_size() == 0 ) return true ;
			if( !hdr_.config().policy(0).use_compact_index(1).has_number_of_slices() ) return true ;
			return enough_inputs( hdr_.config().policy(0).use_compact_index(1).number_of_slices() ) ;
		}

		//! \brief returns the merged headers
		//! Redundant information is removed from the headers (e.g.
		//! repeated paths), the rest is merged as best as possible.
		virtual output::Header read_header() { return hdr_ ; }

		//! \brief reports on one sequence
		//! Enough information is read until one sequence has been
		//! described by each input stream.  Everything else is buffered
		//! if necessary.
		virtual output::Result read_result()
		{
			bool cont = true ;
			for(;;)
			{
				if( !cur_input_ ) cont = false ;
				Result r = streams_[cur_input_]->read_result() ;
				if( ++cur_input_ == streams_.size() ) cur_input_ = 0 ;
				if( r.has_seqid() ) 
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
						Result q = p.second ;
						buffer_.erase( r.seqid() ) ;

						if( ((nread_ + nwritten_) & 0xFFFF) == 0 ) 
							clog << "\033[KRead " << nread_ << ", wrote "
								 << nwritten_ << ", buffering "
								 << buffer_.size() << '\r' << flush ;

						return q ;
					}
				}
			} while( cur_input_ || cont ) ;

			Result q = buffer_.begin()->second.second ;
			buffer_.erase( buffer_.begin()->first ) ;
			return q ;
		}

		//! \brief reads a footer
		//! All footers are merged, the exit codes are logically OR'ed,
		//! and the LSB is set if something wen't wrong internally.
		virtual output::Footer read_footer()
		{
			output::Footer foot ;
			for( size_t i = 0 ; i != streams_.size() ; ++i )
				merge_sensibly( foot, streams_[i]->read_footer() ) ;
			return foot ;
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

		virtual output::Header read_header() { return hdr_ ; }
		virtual output::Footer read_footer() { return foot_ ; }
		virtual output::Result read_result() {
			if( cur_ == end_ ) return Result() ;
			else return *(*cur_++) ;
		}
} ;


// keep queues of streams ordered and separate, so we don't merge the
// same stuff repeatedly
typedef std::deque< std::pair< Header, Stream* > > MergeableQueue ;
typedef std::map< int, MergeableQueue > MergeableQueues ;
MergeableQueues mergeable_queues ;

typedef std::deque< Result* > ScratchSpace ;
ScratchSpace scratch_space ;
Header scratch_header ;
Footer scratch_footer ;
int64_t total_scratch_size = 0 ;

void sort_scratch() {
    if( scratch_space.size() > 1 )
		clog << "Sorting " << scratch_space.size() << " results in memory" << endl ;
    sort( scratch_space.begin(), scratch_space.end(), by_genome_coordinate() ) ;
}

int mktempfile( std::string &name )
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
	name = n1 ;
	return fd ;
}

void enqueue_stream( auto_ptr<Stream>, int = 0 ) ;

void flush_scratch() {
	sort_scratch() ;
	std::string tempname ;
	int fd = mktempfile( tempname ) ;
	scratch_header.set_is_sorted_by_coordinate( true ) ;
	StreamAdapter< std::deque< Result* >::const_iterator >
		sa( scratch_header, scratch_space.begin(), scratch_space.end(), scratch_footer ) ;
	clog << "Writing to tempfile" << endl ;
	write_stream_to_file( fd, sa ) ;
	auto_ptr<Stream> fs( new FileStream( tempname.c_str() ) ) ;
	enqueue_stream( fs ) ;
	unlink( tempname.c_str() ) ;

	std::for_each( scratch_space.begin(), scratch_space.end(), delete_ptr<Result>() ) ;
	scratch_space.clear() ;
	total_scratch_size = 0 ;
}

void enqueue_stream( std::auto_ptr<Stream> s, int level ) 
{
	Header h = s->read_header() ;
	if( h.is_sorted_by_coordinate() ) {
		mergeable_queues[ level ].push_back( make_pair( h, s.release() ) ) ;

		if( FileStream::num_open_files() > max_que_size ) {
			// get the biggest bin, we'll merge everything below that
			int max_bin = 0 ;
			for( MergeableQueues::const_iterator i = mergeable_queues.begin() ;
					i != mergeable_queues.end() ; ++i ) 
				if( i->second.size() > max_bin )
					max_bin = i->first ;

			MergeStream ms ;
			for( MergeableQueues::iterator i = mergeable_queues.begin() ; i->first <= max_bin ; ++i ) 
			{
				for( size_t j = 0 ; j != i->second.size() ; ++j )
					ms.add_stream( i->second[j].first, auto_ptr<Stream>( i->second[j].second ) ) ;
				i->second.clear() ;
			}

			std::string fname ;
			int fd = mktempfile( fname ) ;
			write_stream_to_file( fd, ms ) ;
			auto_ptr<Stream> fs( new FileStream( fname.c_str() ) ) ;
			enqueue_stream( fs, max_bin + 1 ) ;
			unlink( fname.c_str() ) ;
		}
	}
	else
	{
		scratch_header.MergeFrom( h ) ;
		for(;;) {
			Result r = s->read_result() ;
			if( !r.has_seqid() ) break ;

			scratch_space.push_back( new Result( r ) ) ;
			total_scratch_size += r.SpaceUsed() ;
			if( total_scratch_size >= max_arr_size ) flush_scratch() ;
	    }
		scratch_footer.MergeFrom( s->read_footer() ) ;
	}
}

int main_( int argc, const char **argv )
{
	if( argc == 1 ) return 0 ;

	// iterate over command line, glob everything
	glob_t the_glob ;
	glob( argv[1], GLOB_NOSORT, 0, &the_glob ) ;
    for( const char **arg = argv+2 ; arg != argv+argc ; ++arg )
		glob( *arg, GLOB_NOSORT | GLOB_APPEND, 0, &the_glob ) ;

	// iterate over glob results, collect into BestHitStreams
	std::map< int, BestHitStream* > stream_per_slice ;
	for( char **arg = the_glob.gl_pathv ; 
			arg != the_glob.gl_pathv + the_glob.gl_pathc ; ++arg )
	{
		std::auto_ptr<Stream> s (new FileStream( *arg )) ;
		Header h = s->read_header() ;
		if( h.has_sge_slicing_stride() ) {
			BestHitStream* &bhs = stream_per_slice[ h.sge_slicing_index(0) ] ;
			if( !bhs ) bhs = new BestHitStream ;
			bhs->add_stream( h, s ) ;

			// if a BestHitStream is ready, add it to the queue of
			// mergeable streams
			if( bhs->enough_inputs() ) {
				enqueue_stream( auto_ptr<Stream>(bhs) ) ;
				stream_per_slice.erase( h.sge_slicing_index(0) ) ;
			}
		}
		// if no slicing was done, add to queue of mergeable streams
		else enqueue_stream( s ) ;
	}

	// at the end, merge everything
	MergeStream final_stream ;
	auto_ptr<Stream> as( new StreamAdapter< std::deque< Result* >::const_iterator >(
				scratch_header, scratch_space.begin(), scratch_space.end(), scratch_footer ) ) ;
	final_stream.add_stream( as->read_header(), as ) ;

	for( MergeableQueues::const_iterator i = mergeable_queues.begin() ;
			i != mergeable_queues.end() ; ++i )
	{
		for( MergeableQueue::const_iterator j = i->second.begin() ;
				j != i->second.end() ; ++j )
		{
			final_stream.add_stream( j->first, auto_ptr<Stream>( j->second ) ) ;
		}
	}
    clog << "Merging everything to output." << endl ;
	return write_stream_to_file( 0, final_stream, true ) ;
}


