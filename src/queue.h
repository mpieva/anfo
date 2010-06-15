/*
inspired by
c-pthread-queue - c implementation of a bounded buffer queue using posix threads
Copyright (C) 2008  Matthew Dickinson

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef INCLUDED_QUEUE_H
#define INCLUDED_QUEUE_H

#include "stream.h"

#include <pthread.h>

namespace streams {


//! \brief multithreaded adapter for streams
//! This is our primitive for multithreading: something that looks like
//! an ordinary stream, but in reality drives stream processors in
//! separate threads.
//!
//! \todo we need to know if we need a header.  It's awkward otherwise.
//!
//! Here's how it works: we initialize with a list of streams.  They
//! should probably all be clones of a single stream processor,
//! otherwise weird stuff might happen, but we don't enforce that.
//!
//! We now check the state of every stream.  If at least one needs a
//! header, we don't do anything until we get a header.  When we get the
//! header, it is distributed to every stream, then each stream is
//! started in its own thread.  If no header is wanted, we start the
//! threads immediately.
//!
//! We maintain two queues of 'Result' records, incoming and outgoing.
//! Each thread checks the stream status:  if output is available, it is
//! put into outgoing, the thread sleeps if there's no room.  If input is needed,
//! it's taken from incoming, the thread sleeps if nothing's there.  If the stream
//! ends, a marker is put into outgoing and the thread is terminated.  If
//! the input ends, the end marker is put back(!), and the footer is
//! processed (which means we must first put a footer in, then terminate
//! the input queue).
//!
//! The main thread decides what to do depending on the queues.  If
//! there's room in incoming and we haven't received a footer, we ask for
//! input.  Else if output is available, we offer it.  Else we wait (on
//! the combined condition).  If input ends, we store the footer and enter a zero into
//! incoming.  If output ends often enough (it does so once per thread),
//! we wait on all threads, collect the footers, combine them with error
//! codes and signal end of stream.  (XXX We might need a state for "have
//! output *and* need input".)
//! 
//! If a result is asked for (from ELK), we collect all of them into a
//! list.  No interlock is necessary, as long is the result is only
//! demanded after processing finishes.
//!
//! XXX result collection is missing
//! XXX error handling is basically absent

class ConcurrentStream : public Stream
{
	private:
		enum { size_ = 16 } ;
		vector< pair< StreamHolder, ConcurrentStream* > > streams_ ;

		Result *incoming_[ size_ ] ;
		Result *outgoing_[ size_ ] ;

		pthread_cond_t input_available_ ;
		pthread_cond_t room_for_output_ ;
		pthread_cond_t action_needed_ ;

		int next_in_, next_in_free_ ; 		// equal iff incoming is empty
		int next_out_, next_out_free_ ;		// equal iff outgoing is empty

		pthread_mutex_t mutex_ ;
		unsigned eos_outstanding_ ;

		bool incoming_empty() const { return next_in_ == next_in_free_ ; }
		bool outgoing_empty() const { return next_out_ == next_out_free_ ; }
		bool incoming_full() const { return (next_in_+1) % size_ == next_in_free_ ; }
		bool outgoing_full() const { return (next_out_+1) % size_ == next_out_free_ ; }
		bool incoming_terminated() const { return !incoming_empty() && !incoming_[next_in_] ; }

	public:
		//! \brief starts streams in parallel in the background
		//! This takes a list of streams (need not be identical, but the
		//! results will be weird if they aren't) and wires them in
		//! parallel.  Every stream gets its own thread.
		template< typename Iter > ConcurrentStream( Iter begin, Iter end ) :
			next_in_(0), next_in_free_(0),
			next_out_(0), next_out_free_(0),
			eos_outstanding_(0)
		{
			pthread_cond_init( &input_available_, 0 ) ;
			pthread_cond_init( &room_for_output_, 0 ) ;
			pthread_cond_init( &action_needed_, 0 ) ;
			pthread_mutex_init( &mutex_, 0 ) ;

			for( ; begin != end ; ++begin )
				streams_.push_back( make_pair( *begin, this ) ) ;
		}

		~ConcurrentStream() 
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
		//! change due to background activity.
		//! \todo We may need a way to indicate "both input need and
		//!       output available".
		virtual state get_state() 
		{
			pthread_mutex_lock( &mutex_ ) ;
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
					}
					if( !eos_outstanding_ ) { s = end_of_stream ; break ; }
					else if( !outgoing_empty() ) { s = have_output ; break ; }
				}
				pthread_cond_wait( &action_needed_, &mutex_ ) ;
			}
			pthread_mutex_unlock( &mutex_ ) ;
			return s ;
		}

		// Putting a header means we put the header to each stream, then
		// fire off separate threads for each one.
		virtual void put_header( const Header& h )
		{
			if( !eos_outstanding_ )
			{
				for( size_t i = 0 ; i != streams_.size() ; ++i )
				{
					streams_[i].first->put_header( h ) ;

					pthread_t tid ;
					pthread_create( &tid, 0, &start_routine, &streams_[i] ) ;
					pthread_detach( tid ) ;
				}
				eos_outstanding_ = streams_.size() ;
			}
		}

		// we know we have room for input, so no waiting; just put it
		// there and signal
		virtual void put_result( const Result& r )
		{
			pthread_mutex_lock( &mutex_ ) ;
			incoming_[ next_in_free_++ ] = new Result( r ) ;
			next_in_free_ %= size_ ;
			pthread_mutex_unlock( &mutex_ ) ;
			pthread_cond_signal( &input_available_ ) ;
		}

		// Putting a footer means we terminate the input queue and store
		// the footer for later.  Here we don't need to wait for a slot
		// in the incoming queue, because that condition *must* be
		// checked by the caller and the queue cannot fill up by itself.
		virtual void put_footer( const Footer& f )
		{
			assert( !incoming_full() ) ;
			pthread_mutex_lock( &mutex_ ) ;
			foot_ = f ;
			incoming_[ next_in_free_ ] = 0 ;
			next_in_free_ = (next_in_free_+1) % size_ ;
			pthread_mutex_unlock( &mutex_ ) ;
			pthread_cond_signal( &input_available_ ) ;
		}

		// we fetch each header, then merge them.
		// (XXX if there is anything anywhere that needs locking around
		// reading the header, it needs to be implemented there; but
		// that's unlikely to be needed).
		virtual Header fetch_header()
		{
			Header h ;
			for( size_t i = 0 ; i != streams_.size() ; ++i )
				h.MergeFrom( streams_[i].first->fetch_header() ) ;
			return h ;
		}

		// we know output is available, so no waiting, we just take it
		// off
		virtual Result fetch_result()
		{
			pthread_mutex_lock( &mutex_ ) ;
			auto_ptr< Result > r( outgoing_[ next_out_++ ] ) ;
			next_out_ %= size_ ;
			pthread_mutex_unlock( &mutex_ ) ;
			pthread_cond_signal( &room_for_output_ ) ;
			return *r ;
		}

		// to fetch the footer, output must already have been
		// terminated.  that means threads aren't running anymore and we
		// can simply get the footers and merge them.
		virtual Footer fetch_footer()
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


		static void *start_routine( void* param ) {
			pair< StreamHolder, ConcurrentStream* > *p = 
				(pair< StreamHolder, ConcurrentStream* >*) param ;
			return p->second->run_thread( &*p->first ) ;
		}

		//! \brief dequeues a record from incoming and passes it on
		//! We dequeue a record (waiting for it if necessary).  If we
		//! got a record, we unlock and pass it on.  Else we leave the
		//! terminator in and fetch the footer before releasing the
		//! mutex.
		void dequeue( Stream *s )
		{
			pthread_mutex_lock( &mutex_ ) ;
			while( incoming_empty() )
				pthread_cond_wait( &input_available_, &mutex_ ) ;

			auto_ptr< Result > r( incoming_[ next_in_ ] ) ;
			if( r.get() ) {
				++next_in_ ;
				next_in_ %= size_ ;
				pthread_mutex_unlock( &mutex_ ) ;
				pthread_cond_signal( &action_needed_ ) ;
				s->put_result( *r ) ;
			}
			else
			{
				// note ordering: putting the footer can cause output to
				// appear, which will fiddle with the queue.  best not
				// nest this code inside out critical section.
				pthread_mutex_unlock( &mutex_ ) ;
				pthread_cond_signal( &input_available_ ) ;
				s->put_footer( foot_ ) ;
			}
		}

		void enqueue( Result* r )
		{
			pthread_mutex_lock( &mutex_ ) ;
			while( outgoing_full() )
				pthread_cond_wait( &room_for_output_, &mutex_ ) ;

			outgoing_[ next_out_free_++ ] = r ;
			next_out_free_ %= size_ ;
			pthread_mutex_unlock( &mutex_ ) ;
			pthread_cond_signal( &action_needed_ ) ;
		}

		// driven by state of stream.  if input is requested, wait for
		// it and supply it; if eos comes up, fetch the footer, put eos
		// back and supply the footer.
		// if output is available, wait for space, then deliver it.  if
		// the stream ends, wait for space, deliver eos, and terminate.
		void* run_thread( Stream* s ) {
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
} ;

} ; // namespace streams

#endif
