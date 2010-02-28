/*
stolen and adapted from
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

#include <pthread.h>

//! \brief interlocked queue for thread communication
//! This is a generic queue with statically limited capacity that blocks
//! on attempts to dequeue from an empty queue or enqueue into a full
//! one.
template< typename T, size_t capacity > class Queue
{
	private:
		T buffer[capacity] ;
		size_t size;
		size_t in;
		size_t out;
		pthread_mutex_t mutex;
		pthread_cond_t cond_full;
		pthread_cond_t cond_empty;

	public:
		//! \brief initializes empty queue
		Queue()
			: size(0), in(0), out(0) 
		{
			pthread_mutex_init( &mutex, 0 ) ;
			pthread_cond_init( &cond_full, 0 ) ;
			pthread_cond_init( &cond_empty, 0 ) ;
		}

		//! \brief tries to enqueue a value
		//! If the queue is full, the element is not enqueued and false
		//! is returned.
		//! \param value element to be enqueued
		//! \return true iff value could be enqueued
		bool try_enqueue( const T& value )
		{
			pthread_mutex_lock(&mutex);
			bool have_space = size != capacity ;
			if( have_space ) {
				buffer[in] = value;
				++ size;
				++ in;
				in %= capacity;
			}
			pthread_mutex_unlock(&mutex);
			pthread_cond_broadcast(&cond_empty);
			return have_space ;
		}

		//! \brief enqueues a value
		//! If the queue is full, the enqueuing thread is blocked until
		//! at least one element is dequeued.
		//! \param value element to be enqueued
		void enqueue( const T& value )
		{
			pthread_mutex_lock(&mutex);
			while (size == capacity)
				pthread_cond_wait(&cond_full, &mutex);
			buffer[in] = value;
			++ size;
			++ in;
			in %= capacity;
			pthread_mutex_unlock(&mutex);
			pthread_cond_broadcast(&cond_empty);
		}

		//! \brief dequeus a value
		//! If the queue is empty, the calling thread is blocked until
		//! at least one element is enqueued.
		//! \return the dequeued element
		T dequeue()
		{
			pthread_mutex_lock(&mutex);
			while (size == 0)
				pthread_cond_wait(&cond_empty, &mutex);
			T value = buffer[out];
			-- size;
			++ out;
			out %= capacity;
			pthread_mutex_unlock(&mutex);
			pthread_cond_broadcast(&cond_full);
			return value;
		}

		//! \brief tries to deque an element
		//! If no element is available, false is returned.
		//! \param value space to store the dequed element in
		//! \return true iff an element was dequeued
		bool try_dequeue( T& value )
		{
			pthread_mutex_lock(&mutex);
			bool have_element = size != 0 ;
			if( have_element ) {
				value = buffer[out];
				-- size;
				++ out;
				out %= capacity;
			}
			pthread_mutex_unlock(&mutex);
			pthread_cond_broadcast(&cond_full);
			return have_element ;
		}

		size_t get_size()
		{
			pthread_mutex_lock(&mutex);
			size_t size_ = size;
			pthread_mutex_unlock(&mutex);
			return size_;
		}
} ;

/* XXX
 * A clever multithreaded primitive is missing.  That would be a stream
 * that in reality just pushes work off to a number of identical stream
 * processors in their own threads.  We need essentially two queues for
 * incoming and outgoing data, but combined conditions need to be
 * signalled to create a facade of an ordinary stream.  Hmmm...
 *
 * XXX we need to know if we need a header.  It's awkward otherwise.
 */

// Outline:  We check the state of every stream.  If at least one needs
// a header, we'll simply wait.  When we get the header, it is
// distributed to every stream, then each stream is started in its own
// thread.  If no header is wanted, we  (or immediately, if the header wasn't even necessary).  
//
// We maintain two queues of 'Result' records, incoming and outgoing.
// Each thread checks the stream status:  if output is available, it is
// put into outgoing, sleeping if there's no room.  If input is needed,
// it's taken from incoming, sleeping if nothing's there.  If the stream
// ends, a marker is put into outgoing and the thread is terminated.  If
// the input ends, the end marker is put back(!), and the footer is
// processed (must be available).
//
// The main thread decides what to do depending on the queues.  If
// there's room in incoming (XXX what about the footer?) we ask for
// input.  Else if output is available, we offer it.  Else we wait (on
// the combined condition).  If input ends, we enter a zero into
// incoming.  If output ends often enough (it does so once per thread),
// we wait on all threads, collect the footers, combine them with error
// codes and signal end of stream.
// 
// If a result is asked for (from ELK), we collect all of them into a
// list.
class ConcurrentStream : public Stream
{
	private:
		static const int size_ = 10 ;
		vector< StreamHolder > streams_ ;

		Result *incoming_[10] ;
		Result *outgoing_[10] ;

		pthread_cond_t input_available_ ;
		pthread_cond_t room_for_output_ ;
		pthread_cond_t action_needed_ ;

		int next_in_, next_in_free_ ; 		// equal iff incoming is empty
		int next_out_, next_out_free_ ;		// equal iff outgoing is empty

		pthread_mutex_t mutex_ ;
		vector< pthread_t > threads_ ;
		bool accepting_input_ ;

		bool incoming_empty() const { return next_in_ == next_in_free_ ; }
		bool outgoing_empty() const { return next_out_ == next_out_free_ ; }
		bool incoming_full() const { return (next_in_+1) % size_ == next_in_free_ ; }
		bool outgoing_full() const { return (next_out_+1) % size_ == next_out_free_ ; }

	public:
		// Takes a list of streams (need not be identical, but the
		// results will be weird if they aren't) and wires them in
		// parallel.  Every stream gets its own thread.
		template< typename Iter > ConcurrentStream( Iter begin, Iter end ) :
			streams_( begin, end ), 
			next_in_(0), next_in_free_(0),
			next_out_(0), next_out_free_(0),
			accepting_input_(true)
		{
			pthread_cond_init( input_available_ ) ;
			pthread_cond_init( room_for_output_ ) ;
			pthread_cond_init( action_needed_ ) ;
			pthread_mutex_init( mutex_ ) ;
		}

		virtual state get_state() const
		{
			state s = invalid ;
			pthread_mutex_lock(&mutex);
			if( accepting_input_ && !incoming_full() ) s = need_input ;
			else if( !outgoing_empty() ) s = have_output ;
			else s = end_of_stream ;
			pthread_mutex_unlock(&mutex);
			return s ;
		}

























} ;


#endif
