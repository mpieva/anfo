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

#ifndef INCLUDED_CONCURRENT_STREAM_H
#define INCLUDED_CONCURRENT_STREAM_H

#include "stream.h"

#include <utility>
#include <vector>

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
//! codes and signal end of stream.  
//
//! \todo We might need a state for "have output *and* need input".
//! 
//! If a result is asked for (from ELK), we collect all of them into a
//! list.  No interlock is necessary, as long is the result is only
//! demanded after processing finishes.
//!
//! \todo result collection is missing

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

		void init() ;
		void dequeue( Stream* ) ;
		void enqueue( Result* ) ;
		void* run_thread( Stream* ) ;
		static void *start_routine( void* ) ;

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
			init() ;
			for( ; begin != end ; ++begin )
				streams_.push_back( make_pair( *begin, this ) ) ;
		}

		~ConcurrentStream() ;

		virtual state get_state() ;
		virtual void put_header( const Header& h ) ;
		virtual void put_result( const Result& r ) ;
		virtual void put_footer( const Footer& f ) ;
		virtual Header fetch_header() ;
		virtual Result fetch_result() ;
		virtual Footer fetch_footer() ;
} ;

} ; // namespace streams

#endif
