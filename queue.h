/*
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

template< typename T, int capacity > class Queue
{
	private:
		T buffer[capacity] ;
		int size;
		int in;
		int out;
		pthread_mutex_t mutex;
		pthread_cond_t cond_full;
		pthread_cond_t cond_empty;

	public:
		Queue()
			: size(0), in(0), out(0) 
		{
			pthread_mutex_init( &mutex, 0 ) ;
			pthread_cond_init( &cond_full, 0 ) ;
			pthread_cond_init( &cond_empty, 0 ) ;
		}

		void enqueue( const T& value )
		{
			pthread_mutex_lock(&mutex);
			while (size == capacity)
				pthread_cond_wait(&cond_full, &mutex);
			// printf("enqueue %d\n", *(int *)value);
			buffer[in] = value;
			++ size;
			++ in;
			in %= capacity;
			pthread_mutex_unlock(&mutex);
			pthread_cond_broadcast(&cond_empty);
		}

		T dequeue()
		{
			pthread_mutex_lock(&mutex);
			while (size == 0)
				pthread_cond_wait(&cond_empty, &mutex);
			T value = buffer[out];
			// printf("dequeue %d\n", *(int *)value);
			-- size;
			++ out;
			out %= capacity;
			pthread_mutex_unlock(&mutex);
			pthread_cond_broadcast(&cond_full);
			return value;
		}

		int get_size()
		{
			pthread_mutex_lock(&mutex);
			int size_ = size;
			pthread_mutex_unlock(&mutex);
			return size_;
		}
} ;

#endif
