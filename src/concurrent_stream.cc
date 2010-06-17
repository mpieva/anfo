#include "concurrent_stream.h"

#include "util.h"

#include <cassert>

namespace {
	void mutex_init( pthread_mutex_t& m )
	{
		throw_strerror_if_not_null(
				pthread_mutex_init( &m, 0 ), "pthread_mutex_init" ) ;
	}

	void mutex_lock( pthread_mutex_t& m )
	{
		throw_strerror_if_not_null( 
				pthread_mutex_lock( &m ), "pthread_mutex_lock" ) ;
	}

	void mutex_unlock( pthread_mutex_t& m )
	{
		throw_strerror_if_not_null( 
				pthread_mutex_unlock( &m ), "pthread_mutex_unlock" ) ;
	}

	void cond_init( pthread_cond_t& c )
	{
		throw_strerror_if_not_null(
				pthread_cond_init( &c, 0 ), "pthread_cond_init" ) ;
	}

	void cond_wait( pthread_cond_t& c, pthread_mutex_t& m )
	{
		throw_strerror_if_not_null(
				pthread_cond_wait( &c, &m ), "pthread_cond_wait" ) ;
	}

	void cond_signal( pthread_cond_t& c )
	{
		throw_strerror_if_not_null(
				pthread_cond_signal( &c ), "pthread_cond_signal" ) ;
	}
} ;

namespace streams {

void ConcurrentStream::init()
{
	cond_init( input_available_ ) ;
	cond_init( room_for_output_ ) ;
	cond_init( action_needed_ ) ;
	mutex_init( mutex_ ) ;
}

ConcurrentStream::~ConcurrentStream()
{
	pthread_mutex_destroy( &mutex_ ) ;
	pthread_cond_destroy( &action_needed_ ) ;
	pthread_cond_destroy( &room_for_output_ ) ;
	pthread_cond_destroy( &input_available_ ) ;
	assert( outgoing_empty() ) ;
	assert( incoming_terminated() ) ;
	assert( (next_in_+1) % size_ == next_in_free_ ) ;
}

//! \brief determines stream state
//! State is entirely determined by the queues, so that it can
//! change due to background activity.  We determine the state of the
//! queues, and if neither input is required nor output available, we
//! wait.  Note that we need to dequeue termination markers from the
//! outgoing queue and we need to signal if we did so.
//! \todo We may need a way to indicate "both input need and output
//!       available".
Stream::state ConcurrentStream::get_state() 
{
	mutex_lock( mutex_ ) ;
	state s = invalid ;
	for(;;)
	{
		if( !incoming_full() && !incoming_terminated() ) { s = need_input ; break ; }
		else if( !eos_outstanding_ ) { s = end_of_stream ; break ; }
		else {
			while( !outgoing_empty() && !outgoing_[ next_out_ ] )
			{
				--eos_outstanding_ ;
				++next_out_ ;
				next_out_ %= size_ ;
				cond_signal( room_for_output_ ) ;
			}
			if( !eos_outstanding_ ) { s = end_of_stream ; break ; }
			else if( !outgoing_empty() ) { s = have_output ; break ; }
		}
		cond_wait( action_needed_, mutex_ ) ;
	}
	mutex_unlock( mutex_ ) ;
	return s ;
}

//! Putting a header means we put the header to each stream, then
//! fire off separate threads for each one.  Threads are started in a
//! detached state, we get results from the stream objects themselves.
void ConcurrentStream::put_header( const Header& h )
{
	if( !eos_outstanding_ )
	{
		for( size_t i = 0 ; i != streams_.size() ; ++i )
		{
			streams_[i].first->put_header( h ) ;

			pthread_t tid ;
			throw_strerror_if_not_null(
					pthread_create( &tid, 0, &start_routine, &streams_[i] ),
					"pthread_create" ) ;
			throw_strerror_if_not_null( pthread_detach( tid ), "pthread_detach" ) ;
		}
		eos_outstanding_ = streams_.size() ;
	}
}

//! we know we have room for input, so thre's no waiting, we just put it
//! in incoming and signal
void ConcurrentStream::put_result( const Result& r )
{
	mutex_lock( mutex_ ) ;
	assert( !incoming_full() ) ;
	incoming_[ next_in_free_++ ] = new Result( r ) ;
	next_in_free_ %= size_ ;
	cond_signal( input_available_ ) ;
	mutex_unlock( mutex_ ) ;
}

//! Putting a footer means we terminate the input queue and store the
//! footer for later.  Here we don't need to wait for a slot in the
//! incoming queue, because that condition must be checked by the caller
//! and the queue cannot fill up by itself.
void ConcurrentStream::put_footer( const Footer& f )
{
	mutex_lock( mutex_ ) ;
	assert( !incoming_full() ) ;
	foot_ = f ;
	incoming_[ next_in_free_ ] = 0 ;
	next_in_free_ = (next_in_free_+1) % size_ ;
	cond_signal( input_available_ ) ;
	mutex_unlock( mutex_ ) ;
}

//! We fetch each header, then merge them.
//! \note This assumes that reading the header from any stream is safe
//!       at any time.  While this is unlikely to be a wrong assumption,
//!       it still feels brittle.
Header ConcurrentStream::fetch_header()
{
	Header h ;
	for( size_t i = 0 ; i != streams_.size() ; ++i )
		h.MergeFrom( streams_[i].first->fetch_header() ) ;
	return h ;
}

Result ConcurrentStream::fetch_result()
{
	mutex_lock( mutex_ ) ;
	assert( !outgoing_empty() ) ;
	auto_ptr< Result > r( outgoing_[ next_out_++ ] ) ;
	next_out_ %= size_ ;
	cond_signal( room_for_output_ ) ;
	mutex_unlock( mutex_ ) ;
	return *r ;
}

//! To fetch the footer, output must already have been terminated.  That
//! means threads aren't running anymore and we can simply get the
//! footers and merge them, without regard for concurrent writes.
Footer ConcurrentStream::fetch_footer()
{
	int exit_code = 0 ;
	Footer f ;
	for( size_t i = 0 ; i != streams_.size() ; ++i )
	{
		Footer f1 = streams_[i].first->fetch_footer() ;
		exit_code |= f1.exit_code() ;
		f.MergeFrom( f1 ) ;
	}
	f.set_exit_code( exit_code ) ;
	return f ;
}

void *ConcurrentStream::start_routine( void* param ) {
	try {
		pair< StreamHolder, ConcurrentStream* > *p = 
			(pair< StreamHolder, ConcurrentStream* >*) param ;
		return p->second->run_thread( &*p->first ) ;
	}
	catch( const std::string& e ) { perr( e ) ; }
	catch( const char *e ) { perr( e ) ; }
	catch( char *e ) { perr( e ) ; }
	catch( const Exception& e ) { perr( e ) ; }
	catch( const std::exception& e ) { perr( e.what() ) ; }
	catch( ... ) { perr( "Oh noes!" ) ; }
	return 0 ;
}

//! \brief dequeues a record from incoming and passes it on
//! We dequeue a record, waiting for it if necessary.  If we got a
//! record, we unlock and pass it on.  Else we got a terminator, which
//! we leave in, then we copy the footer to change to the end-of-stream
//! state.
void ConcurrentStream::dequeue( Stream *s )
{
	mutex_lock( mutex_ ) ;
	while( incoming_empty() ) cond_wait( input_available_, mutex_ ) ;

	auto_ptr< Result > r( incoming_[ next_in_ ] ) ;
	if( r.get() ) {
		++next_in_ ;
		next_in_ %= size_ ;
		cond_signal( action_needed_ ) ;
		mutex_unlock( mutex_ ) ;
		s->put_result( *r ) ;
	}
	else
	{
		// note ordering: putting the footer can cause output to
		// appear, which will fiddle with the queue.  best not
		// nest this code inside our critical section.
		cond_signal( input_available_ ) ;
		mutex_unlock( mutex_ ) ;
		s->put_footer( foot_ ) ;
	}
}

void ConcurrentStream::enqueue( Result* r )
{
	mutex_lock( mutex_ ) ;
	while( outgoing_full() ) cond_wait( room_for_output_, mutex_ ) ;

	outgoing_[ next_out_free_++ ] = r ;
	next_out_free_ %= size_ ;
	cond_signal( action_needed_ ) ;
	mutex_unlock( mutex_ ) ;
}

// driven by state of stream.  if input is requested, wait for it and
// supply it; if eos comes up, fetch the footer, put eos back and supply
// the footer.  if output is available, wait for space, then deliver it.
// if the stream ends, wait for space, deliver eos, and terminate.
void* ConcurrentStream::run_thread( Stream* s ) {
	for(;;) {
		switch( s->get_state() ) {
			case invalid:
				throw "stream in invalid state: not supposed to happen!" ;

			case need_input:
				dequeue( s ) ;
				break ;

			case have_output:
				enqueue( new Result( s->fetch_result() ) ) ;
				break ;

			case end_of_stream:
				enqueue( 0 ) ;
				return 0 ;
		}
	}
}

Object ConcurrentStream::get_summary()
{
	GC_Node ;
	Object l = Null ;
	GC_Link( l ) ;
	for( size_t i = streams_.size() ; i != 0 ; --i )
		l = Cons( streams_[i-1].first->get_summary(), l ) ;
	GC_Unlink ;
	return l ;
}

} ;

