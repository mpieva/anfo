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

#endif
