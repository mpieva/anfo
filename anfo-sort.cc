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

#if HAVE_UNISTD_H
#  include <unistd.h>
#endif

#if HAVE_ALLOCA_H
#  include <alloca.h>
#endif

using namespace google::protobuf::io ;
using namespace std ;
using namespace output ;

static unsigned max_que_size = 256 ;
static unsigned max_arr_size = 1024*1024*1024 ;

// read ANFO files, sort and merge them
// Outline: keep a queue of opened files and a deque of results.  A new
// file (e.g. from command line) is added to the queue if it is sorted,
// else it is read and added to the deque.  If the deque becomes large,
// it is sorted and written to tempory storage, which is then added to
// the queue.  If the queue becomes long, it is merged into temporary
// storage, the new file is added to the queue.  Finally, all open files
// and the contents of the deque are merged and written out.

struct by_genome_coordinate {
    bool operator() ( const Result *a, const Result *b ) {
	if( a->has_best_to_genome() && !b->has_best_to_genome() ) return true ;
	if( !a->has_best_to_genome() ) return false ;
	const Hit& u = a->best_to_genome(), v = b->best_to_genome() ;
	if( u.sequence() < v.sequence() ) return true ;
	if( v.sequence() < u.sequence() ) return false ;
	if( u.start_pos() + max(u.aln_length(),0) < v.start_pos() + max(v.aln_length(),0) ) return true ;
	if( v.start_pos() + max(v.aln_length(),0) < u.start_pos() + max(u.aln_length(),0) ) return false ;
	return false ;
    }
} ;

deque< const Result* > arr ;
deque< AnfoFile* > que ;
unsigned total_arr_size = 0 ;

Header hdr ;
Footer ftr ;

void sort_arr() {
    if( arr.size() > 1 ) clog << "sorting " << arr.size() << " results in memory" << endl ;
    sort( arr.begin(), arr.end(), by_genome_coordinate() ) ;
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

void dump_arr() {
    sort_arr() ;
	std::string name ;
    int fd = mktempfile( name ) ;
	{
		FileOutputStream fos( fd ) ;
		fos.SetCloseOnDelete( true ) ;
		std::auto_ptr< ZeroCopyOutputStream > zos( compress_fast( &fos ) ) ;
		CodedOutputStream cos( zos.get() ) ;
		hdr.set_is_sorted_by_coordinate( true ) ;
		cos.WriteRaw( "ANFO", 4 ) ;
		write_delimited_message( cos, 1, hdr ) ;
		clog << "writing temp file" << endl ;
		for( deque< const Result* >::iterator i = arr.begin(), j = arr.end() ; i!=j ; ++i )
		{
			write_delimited_message( cos, 2, **i ) ;
			delete *i ;
		}
		arr.clear() ;
		total_arr_size = 0 ;
		write_delimited_message( cos, 3, ftr ) ;
		que.push_back( new AnfoFile( name, true ) ) ;
		que.back()->read_header() ;
	}
}

void merge_footer( AnfoFile* f ) 
{
    ftr.set_exit_code( ftr.exit_code() | f->read_footer().exit_code() ) ;
}

void merge_all( int fd, bool final ) 
{
	FileOutputStream fos( fd ) ;
	std::auto_ptr< ZeroCopyOutputStream > zos(
                    final ? compress_small( &fos ) : compress_fast( &fos ) ) ;
	CodedOutputStream cos( zos.get() ) ;
	cos.WriteRaw( "ANFO", 4 ) ;
	write_delimited_message( cos, 1 , hdr ) ;

	deque< Result > rs ;
	for( size_t i = 0 ; i != que.size() ; ++i )
	{
		Result r = que[i]->read_result() ;
		if( r.has_seqid() ) rs.push_back( r ) ;
		else
		{
			merge_footer( que[i] ) ;
			delete que[i] ;
			que.erase( que.begin()+i ) ;
		}
	}

	while( rs.size() ) 
	{
		int m = 0 ;
		Result rm = rs[0] ;
		for( size_t i = 1 ; i != rs.size() ; ++i ) 
		{
			if( by_genome_coordinate()( &rs[i], &rm ) )
			{
				m = i ;
				rm = rs[i] ;
			}
		}
		if( final && !arr.empty() && by_genome_coordinate()( arr[0], &rm ) ) {
			m = -1 ;
			rm = *arr[0] ;
		}

		write_delimited_message( cos, 2, rm ) ;
		if( m == -1 ) {
			delete arr.front() ;
			arr.pop_front() ;
		} else {
			rs[m] = que[m]->read_result() ;
			if( !rs[m].has_seqid() ) {
				merge_footer( que[m] ) ;
				delete que[m] ;
				que.erase( que.begin() + m ) ;
				rs.erase( rs.begin() + m ) ;
			}
		}
	}

	for( deque< const Result* >::iterator i = arr.begin(), j = arr.end() ; i!=j ; ++i )
	{
		write_delimited_message( cos, 2, **i ) ;
		delete *i ;
	}
	arr.clear() ;
	total_arr_size = 0 ;
	write_delimited_message( cos, 3, ftr ) ;
}

void flush_queue() {
	std::string name ;
    int fd = mktempfile( name ) ;
    FileOutputStream fos( fd ) ;
	fos.SetCloseOnDelete( true ) ;
	std::auto_ptr< ZeroCopyOutputStream > zos( compress_fast( &fos ) ) ;
    CodedOutputStream cos( zos.get() ) ;
    hdr.set_is_sorted_by_coordinate( true ) ;
    cos.WriteRaw( "ANFO", 4 ) ;
    write_delimited_message( cos, 1, hdr ) ;
    clog << "merging to temp file" << endl ;
    merge_all( fd, false ) ;
    write_delimited_message( cos, 3, ftr ) ;
    que.push_back( new AnfoFile( name, true ) ) ;
    que.back()->read_header() ;
}

int main_( int argc, const char **argv )
{
    for( const char **arg = argv+1 ; arg != argv+argc ; ++arg )
    {
	AnfoFile *f = new AnfoFile( *arg ) ;
	Header h = f->read_header() ;
	hdr.MergeFrom( h ) ;
	if( h.is_sorted_by_coordinate() ) {
	    que.push_back( f ) ;
	    if( que.size() == max_que_size ) flush_queue() ;
	} else {
	    clog << "reading " << *arg << endl ;
	    for(;;) {
		Result r = f->read_result() ;
		if( !r.has_seqid() ) break ;
		arr.push_back( new Result( r ) ) ;
		total_arr_size += r.SpaceUsed() ;
		if( total_arr_size >= max_arr_size ) {
		    dump_arr() ;
		    if( que.size() == max_que_size ) flush_queue() ;
		}
	    }
	    merge_footer( f ) ;
	    delete f ;
	}
    }
			
    sort_arr() ;
    clog << "merging everything to output" << endl ;
    merge_all( 1, true ) ;
    return ftr.exit_code() ;
}


