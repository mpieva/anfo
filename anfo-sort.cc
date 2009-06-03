#include "config.h"

#include "compress_stream.h"
#include "stream.h"
#include "util.h"

#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/text_format.h>

#include <algorithm>
#include <cstdlib>
#include <deque>
#include <iostream>

#include <glob.h>
#include <popt.h>

#if HAVE_ALLOCA_H
#  include <alloca.h>
#endif

#if HAVE_UNISTD_H
#  include <unistd.h>
#endif

using namespace google::protobuf ;
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


//! \brief merges sorted streams into a sorted stream
//! \todo needs to be parameterized with a comparison function (we
//!       already have two).
class MergeStream : public FanInStream
{
	private:
		vector< Stream* > streams_ ;
		deque< Result > rs_ ;
		output::Header hdr_ ;
		output::Footer foot_ ;

	public:
		MergeStream() {}
		virtual ~MergeStream() { for_each( streams_.begin(), streams_.end(), delete_ptr<Stream>() ) ; }

		virtual void add_stream( Stream* s )
		{
			assert( s->get_header().is_sorted_by_coordinate() ) ;
			merge_sensibly( hdr_, s->get_header() ) ;

			Result r ;
			if( s->read_result( r ) ) 
			{
				rs_.push_back( r ) ;
				streams_.push_back( s ) ;
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
class BestHitStream : public FanInStream
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
		virtual void add_stream( Stream* s ) {
			merge_sensibly( hdr_, s->get_header() ) ;
			streams_.push_back( s ) ;
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

void enqueue_stream( Stream*, int = 0 ) ;

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
	
	enqueue_stream( new AnfoFile( fd, tempname.c_str() ) ) ;

	for_each( scratch_space.begin(), scratch_space.end(), delete_ptr<Result>() ) ;
	scratch_space.clear() ;
	total_scratch_size = 0 ;
}

void enqueue_stream( Stream* s, int level ) 
{
	Header h = s->get_header() ;
	if( h.is_sorted_by_coordinate() ) {
		mergeable_queues[ level ].push_back( s ) ;

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
							ms.add_stream( i->second[j] ) ;
						i->second.clear() ;
					}
					write_stream_to_file( fd, ms ) ;
				}

				throw_errno_if_minus1( lseek( fd, 0, SEEK_SET ), "seeking in ", fname.c_str() ) ;
				enqueue_stream( new AnfoFile( fd, fname.c_str() ), max_bin + 1 ) ;
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
		Stream* s = new AnfoFile( *arg ) ;
		Header h = s->get_header() ;
		if( h.has_sge_slicing_stride() ) {
			BestHitStream* &bhs = stream_per_slice[ h.sge_slicing_index(0) ] ;
			if( !bhs ) bhs = new BestHitStream ;
			bhs->add_stream( s ) ;

			// if a BestHitStream is ready, add it to the queue of
			// mergeable streams
			if( bhs->enough_inputs() ) {
				clog << "\033[KGot everything for slice " << h.sge_slicing_index(0) << endl ;
				enqueue_stream( bhs ) ;
				stream_per_slice.erase( h.sge_slicing_index(0) ) ;
			}
		}
		// if no slicing was done, add to queue of mergeable streams
		else enqueue_stream( s ) ;
	}
	if( !stream_per_slice.empty() ) 
		cerr << "\033[KWARNING: input appears to be incomplete" << endl ;

	for( map< int, BestHitStream* >::iterator l = stream_per_slice.begin(), r = stream_per_slice.end() ;
			l != r ; ++l ) enqueue_stream( l->second ) ;

	MergeStream final_stream ;

	// at the end, merge everything
	if( scratch_space.begin() != scratch_space.end() ) {
		clog << "\033[Kfinal sort" << endl ;
		sort_scratch() ;
		clog << "\033[Kdone sorting" << endl ;
		final_stream.add_stream( new StreamAdapter< deque< Result* >::const_iterator >(
					scratch_header, scratch_space.begin(), scratch_space.end(), scratch_footer ) ) ;
	}

	for( MergeableQueues::const_iterator i = mergeable_queues.begin() ; i != mergeable_queues.end() ; ++i )
		for( MergeableQueue::const_iterator j = i->second.begin() ; j != i->second.end() ; ++j )
			final_stream.add_stream( *j ) ;
	
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

class RepairHeaderStream : public Stream
{
	private:
		std::auto_ptr< Stream > s_ ;
		std::auto_ptr< output::Header > h_ ;

	public:
		RepairHeaderStream( Stream* s ) : s_(s) {}
		virtual ~RepairHeaderStream() {}

		virtual const output::Header& get_header() ;
		virtual bool read_result( output::Result& r ) { return s_->read_result( r ) ; }
		virtual const output::Footer& get_footer() { return s_->get_footer() ; }
} ;

const output::Header& RepairHeaderStream::get_header() 
{
	if( !h_.get() ) {
		char tmpname[] = "/tmp/anfo_header_XXXXXX" ;
		int fd = mkstemp( tmpname ) ;
		{
			FileOutputStream fos( fd ) ;
			TextFormat::Print( s_->get_header(), &fos ) ;
		}

		string cmd = string( getenv("EDITOR") ) + " " + tmpname ;
		for(;;) {
			system( cmd.c_str() ) ;
			h_.reset( new Header ) ;
			lseek( fd, 0, SEEK_SET ) ;
			FileInputStream fis( fd ) ;
			if( TextFormat::Parse( &fis, h_.get() ) ) break ;
		} 
	}
	return *h_ ;
}

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

typedef Stream* (*FilterMaker)( float, float, const char*, const char*, Stream* ) ;

struct FilterParams {
	float slope, intercept ;
	const char *genome, *arg ;
	FilterMaker maker ;
} ;

typedef std::vector< FilterParams > FilterStack ;
Stream* run_filter_stack( const FilterStack& st, Stream* in )
{
	for( FilterStack::const_iterator l = st.begin(), r = st.end() ; l != r ; ++l )
		in = (l->maker)( l->slope, l->intercept, l->genome, l->arg, in ) ;
	return in ;
}

Stream* mk_sort_by_pos( float, float, const char* genome, const char*, Stream* s ) {}
Stream* mk_sort_by_name( float, float, const char*, const char*, Stream* s ) {}

Stream* mk_filter_by_length( float, float, const char*, const char* arg, Stream* in )
{ return new LengthFilterStream( in, atoi(arg) ) ; }

Stream* mk_filter_by_score( float s, float i, const char* g, const char*, Stream* in )
{ return new ScoreFilterStream( in, s, i, g ) ; }

Stream* mk_filter_by_hit( float, float, const char*, const char* arg, Stream* s )
{ return new HitFilterStream( s, arg ) ; }

Stream* mk_edit_header( float, float, const char*, const char* arg, Stream* s )
{ return new RepairHeaderStream( s ) ; }

Stream* mk_rmdup( float, float, const char*, const char*, Stream* s )
{ return new RmdupStream( s ) ; }

Stream* mk_output( float, float, const char*, const char* arg, Stream* s ) {}
Stream* mk_output_sam( float, float, const char*, const char* arg, Stream* s ) {}
Stream* mk_output_fasta( float, float, const char*, const char* arg, Stream* s ) {}
Stream* mk_output_glz( float, float, const char*, const char* arg, Stream* s ) {}

int main_( int argc, const char **argv )
{
	GOOGLE_PROTOBUF_VERIFY_VERSION ;
	enum { opt_none, opt_sort_pos, opt_sort_name, opt_filter_length,
		opt_filter_score, opt_filter_hit, opt_edit_header, opt_merge,
		opt_mega_merge, opt_concat, opt_rmdup, opt_output, opt_output_sam,
		opt_output_fasta, opt_output_glz, opt_version } ;

	FilterParams params = { 0, 0, 0, 0 } ;
	int param_noise = 1 ;

	struct poptOption options[] = {
		{ "sort-pos",      's', POPT_ARG_NONE,   0, opt_sort_pos,      "sort by alignment position", 0 },
		{ "sort-name",     'S', POPT_ARG_NONE,   0, opt_sort_name,     "sort by read name", 0 },
		{ "filter-length", 'l', POPT_ARG_INT,    0, opt_filter_length, "filter for length of at least L", "L" },
		{ "filter-score",  'f', POPT_ARG_NONE,   0, opt_filter_score,  "filter for max score", 0 },
		{ "filter-hit",    'h', POPT_ARG_NONE,   0, opt_filter_hit,    "filter for having a hit", 0 },
		{ "edit-header",    0 , POPT_ARG_NONE,   0, opt_edit_header,   "invoke the text editor on the stream's header", 0 },
		{ "concat",        'c', POPT_ARG_NONE,   0, opt_concat,        "concatenate streams", 0 },
		{ "merge",         'm', POPT_ARG_NONE,   0, opt_merge,         "merge streams to retain best hits", 0 },
		{ "mega-merge",     0 , POPT_ARG_NONE,   0, opt_mega_merge,    "merge many streams, e.g. from grid", 0 },
		{ "rmdup",         'd', POPT_ARG_NONE,   0, opt_rmdup,         "remove dups", 0 },
		{ "output",        'o', POPT_ARG_STRING, 0, opt_output,        "write native stream to file FILE", "FILE" },
		{ "output-sam",     0 , POPT_ARG_STRING, 0, opt_output_sam,    "write alns in sam format to FILE", "FILE" },
		{ "output-fasta",   0 , POPT_ARG_STRING, 0, opt_output_fasta,  "write alignments in fasta format to FILE", "FILE" },
		{ "output-glz",     0 , POPT_ARG_STRING, 0, opt_output_glz,    "write consensus in glz format to FILE", "FILE" },

		{ "set-slope",      0 , POPT_ARG_FLOAT,  &params.slope,     0, "set slope for subsequent filters to S", "S" },
		{ "set-intercept",  0 , POPT_ARG_FLOAT,  &params.intercept, 0, "set length intercept for filters to L", "L" },
		{ "set-genome",     0 , POPT_ARG_STRING, &params.genome,    0, "set interesting genome for most operations to G", "G" },
		{ "clear-genome",   0 , POPT_ARG_VAL,    &params.genome,    0, "clear interesting genome", 0 },

		{ "quiet",         'q', POPT_ARG_VAL,    &param_noise,      0, "suppress most output", 0 },
		{ "verbose",       'v', POPT_ARG_VAL,    &param_noise,      2, "produce more output", 0 },
		{ "version",       'V', POPT_ARG_NONE,   0, opt_version,       "print version number and exit", 0 },
		POPT_AUTOHELP POPT_TABLEEND
	} ;

	FilterMaker filter_makers[] = { 0, mk_sort_by_pos, mk_sort_by_name,
		mk_filter_by_length, mk_filter_by_score, mk_filter_by_hit,
		mk_edit_header, 0, 0, 0, mk_rmdup, mk_output, mk_output_sam,
		mk_output_fasta, mk_output_glz } ;

	FilterStack filters_initial, filters_terminal ;
	FilterStack *filters_current = &filters_initial ;
	std::auto_ptr< FanInStream > merging_stream ;

	poptContext pc = poptGetContext( "anfo", argc, argv, options, 0 ) ;
	poptSetOtherOptionHelp( pc, "[OPTION...] [sequence-file...]" ) ;
	poptReadDefaultConfig( pc, 0 ) ;

	if( argc <= 1 ) { poptPrintHelp( pc, stderr, 0 ) ; return 1 ; }
	for( int rc = poptGetNextOpt( pc ) ; rc > 0 ; rc = poptGetNextOpt(pc) )
	{
		if( rc == opt_version ) {
			std::cout << poptGetInvocationName(pc) << ", revision " << PACKAGE_VERSION << std::endl ;
			return 0 ;
		}
		else if( rc >= 0 && rc * sizeof(filter_makers[0]) < sizeof(filter_makers) && filter_makers[rc] )
		{
			filters_current->push_back( params ) ;
			filters_current->back().arg = poptGetOptArg( pc ) ;
			filters_current->back().maker = filter_makers[rc] ;
		}
		else if( rc == opt_merge ) {
			merging_stream.reset( new BestHitStream ) ;
			filters_current = &filters_terminal ;
		}
		else if( rc == opt_mega_merge ) {
			merging_stream.reset( new MegaMergeStream ) ;
			filters_current = &filters_terminal ;
		}
		else if( rc == opt_concat ) {
			merging_stream.reset( new ConcatStream ) ;
			filters_current = &filters_terminal ;
		}
		else {
			std::clog << poptGetInvocationName(pc) << ": " << poptStrerror( rc ) 
				<< ' ' << poptBadOption( pc, 0 ) << std::endl ;
			return 1 ; 
		}
	}

	if( !merging_stream.get() ) merging_stream.reset( new ConcatStream ) ;

	while( const char* arg = poptGetArg( pc ) ) {
		merging_stream->add_stream( run_filter_stack( filters_initial, new AnfoFile( arg ) ) ) ;
	}
	std::auto_ptr< Stream > final_stream( run_filter_stack( filters_terminal, merging_stream.release() ) ) ;
	for( Result r ; final_stream->read_result( r ) ; ) {
		// XXX build some stats?
	}
	poptFreeContext( pc ) ;
	return 0 ;
}

