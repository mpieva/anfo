#ifndef INCLUDED_ALIGN_H
#define INCLUDED_ALIGN_H

#include "index.h"
#include "judy++.h"

#include <deque>

#if INT_BITS <= 32 
#warning "I will compile, but this code won't work on a 32 bit machine"
#endif

/* Alignment by Dijkstra's Algorithm
 *
 * The idea, in part already already realized in TWA, is to treat the
 * alignment problem as finding the shortest path through a graph, in
 * this case a regular grid.  In our case, we have many alignments and
 * are only interested in the best one (or maybe a few exceptionally
 * good ones).  So what we're looking for is the shortest path from a
 * number of starting places (the seeds) to a number of goals (any
 * completed semi-global alignment).  This can all be done in parallel,
 * and the single loop around the algorithm terminates as soon as at
 * least one alignment is done.  There's no need to consider all the
 * marginal alignments.  (Mega-)blast has to do the latter, so there's
 * hope we can get a lot faster by avoiding all this.
 *
 * To guarantee that we find the best alignment first, all costs (aka
 * scores) have to be positive.  To put all this on a solid base, the
 * scores are the negative natural logarithms of probabilities.
 * Naturally, all scores are now positive and have a sensible
 * interpretation, though a different one from the blast scores.  The
 * algorithm will run a bit quicker if we subtract a small constant from
 * the scores so the minimum cost (for perfect matches) goes to zero.
 *
 *
 * Potential alignments are kept in a priority queue, obviously sorted
 * by score.  We will sort alignments of equal score additionally by the
 * number of remaining nucleotides.  That way, advanced alignments will
 * be considered first.  The idea here is that we quickly find the best
 * alignment from the best seed and are done, maybe without even
 * /touching/ the other seeds.
 */


/* Dijkstra's Algorithm.
 *
 * We need a sparse bitmap of "closed" states, a priority queue of
 * "open" states, and a representation of states themselves.  We choose
 * to encode the complete automaton into the state.  As much as possible
 * of the automaton should be static information, which is quite
 * possible with only one type of automaton running at the same time.
 */

/* This automaton realizes two alignments (to the left and to the
 * right of a seed) in sequence, allowing for affine gap penalties and
 * modelling ancient DNA by separate states for double stranded and
 * single stranded parts.
 * XXX: To be implemented.
 */
struct simple_adna {} ;

/* Extremely simple automaton for testing.  There are no internal
 * states, gap costs are linear, aDNA is not modelled, but the alignment
 * consists of two parts in series.
 */

struct flat_alignment {
	// What to store for an alignment?  We'll ultimately need a place to
	// set up on when doing backtracing; this will be within the seed.
	// We can define that by just two DnaPs.  Only every second position in the
	// genome can be represented by a plain pointer, but we can always
	// choose such a position within the seed.  The same is true for the
	// sequence itself, but we cannot choose freely here.  Instead,
	// we'll store the sequence twice, if necessary.  A
	// reverse-complemented copy of the sequence is made for reversed
	// alignments; other wise they work the same way.  Sequences are
	// self-delimiting, so the two pointers already define the alignment
	// task---good!.  For the ongoing state, we need two offsets (one on
	// each sequence) and a bit to know whether we're already in the
	// second part.  We'll use signed shorts for the offsets; ±32k
	// should be plenty for the intended application of a mapper.
	// Furthermore, the score accumulated so far must be stored, and it
	// is here where we bite off a bit.  That limits the score (cost,
	// actually) to 2G (not really a limit), and we can store the whole
	// state in 4 double words (16 byte in a 32-bit memory model, 32
	// byte in a 64 bit memory model.

	DnaP reference, query ;					// entry point
	signed short ref_offs, query_offs ; 	// current positions
	uint32_t skew  :  1 ;					// initial offset for query
	uint32_t state :  1 ;					// and state (left/right)
	uint32_t score : 30 ;					// cost so far

	// associated types: 
	// - set of closed nodes (ClosedSet :: flat_alignment -> Bool)
	// - map of closed nodes to ancestor states
	//   (ClosedMap :: flat_alignment -> Maybe flat_alignment)
	typedef JudyL< JudyL< JudyL< Judy1 > > > ClosedSet ;
	typedef JudyL< JudyL< JudyL< JudyL< const flat_alignment* > > > > ClosedMap ;

	static int subst_mat( Ambicode a, Ambicode b ) {
		// this is a trivial one: any overlap gives a point, everything
		// else costs.
		return (a & b) != 0 ? 1 : -2 ;
	}
	static int gap_penalty() {
		return -2 ;
	}
} ;

// basic operations (overloaded for different automata)
bool is_empty( const std::deque< flat_alignment >& ol ) { return ol.empty() ; }
flat_alignment get_top( const std::deque< flat_alignment >& ol ) { return ol.front() ; }

void pop_top( std::deque< flat_alignment >& ol ) { std::pop_heap( ol.begin(), ol.end() ) ; ol.pop_back() ; }

bool is_set_bit( flat_alignment::ClosedSet& cl, const flat_alignment& s )
{
	return cl.insert( (Word_t)s.reference.unsafe_ptr() )
		     ->insert( (Word_t)s.query.unsafe_ptr() )
			 ->insert( s.ref_offs )
			 ->test( s.query_offs | s.state << 16 ) ;
}

bool is_present( flat_alignment::ClosedMap& cl, const flat_alignment& s )
{
	return cl.insert( (Word_t)s.reference.unsafe_ptr() )
		     ->insert( (Word_t)s.query.unsafe_ptr() )
			 ->insert( s.ref_offs )
			 ->insert( s.query_offs | s.state << 16 ) ;
}

void set_bit( flat_alignment::ClosedSet& cl, const flat_alignment& s )
{
	cl.insert( (Word_t)s.reference.unsafe_ptr() )
	  ->insert( (Word_t)s.query.unsafe_ptr() )
	  ->insert( s.ref_offs )
	  ->set( s.query_offs | s.state << 16 ) ;
}

void insert( flat_alignment::ClosedMap& cl, const flat_alignment& s, const flat_alignment *p )
{ 
	*cl.insert( (Word_t)s.reference.unsafe_ptr() )
	   ->insert( (Word_t)s.query.unsafe_ptr() )
	   ->insert( s.ref_offs )
	   ->insert( s.query_offs | s.state << 16 ) = p ;
}

bool finished( const flat_alignment& s ) {
	return s.state == 1 
		&& s.reference[ s.ref_offs ] == 0
		&& s.query[ s.query_offs ] == 0 ;
}

template< typename F > void greedy( flat_alignment& s )
{
	for(;;) 
	{
		if( s.reference[s.ref_offs] == 0 || s.query[s.query_offs] == 0 ) 
		{
			if( s.state == 1 ) return ;
			s.state = 0 ;
			s.ref_offs = 0 ;
			s.query_offs = 0 ;
		}
		else if( s.reference[s.ref_offs] == s.query[s.query_offs] )
		{
			if( s.state == 0 ) 
			{
				++s.ref_offs ;
				++s.query_offs ;
			}
			else
			{
				--s.ref_offs ;
				--s.query_offs ;
			}
		}
	}
}

// XXX this is all very ugly; need to abstract the gunk out of here.
template< typename F > void forward( const flat_alignment& s, F f )
{
	// what to do?  
	// in state 0: mismatch, gap ref, gap query
	// in state 1: mismatch, gap ref, gap query
	// In any case, we can rely on the fact that greedy() was called
	// before.
	if( s.state == 0 )
	{
		{
			flat_alignment s1 = s ;
			s1.score += flat_alignment::subst_mat( s1.reference[s1.ref_offs], s1.query[s1.query_offs] ) ;
			++s1.ref_offs ;
			++s1.query_offs ;
			f( s1 ) ;
		}{
			flat_alignment s2 = s ;
			s2.score += flat_alignment::gap_penalty() ;
			++s2.ref_offs ;
			f( s2 ) ;
		}{
			flat_alignment s3 = s ;
			s3.score += flat_alignment::gap_penalty() ;
			++s3.query_offs ;
			f( s3 ) ;
		}
	}
	else
	{
		{
			flat_alignment s1 = s ;
			s1.score += flat_alignment::subst_mat( s1.reference[s1.ref_offs], s1.query[s1.query_offs] ) ;
			--s1.ref_offs ;
			--s1.query_offs ;
			f( s1 ) ;
		}{
			flat_alignment s2 = s ;
			s2.score += flat_alignment::gap_penalty() ;
			--s2.ref_offs ;
			f( s2 ) ;
		}{
			flat_alignment s3 = s ;
			s3.score += flat_alignment::gap_penalty() ;
			--s3.query_offs ;
			f( s3 ) ;
		}
	}
}

// a heap puts the largest(!) element at the top, so we want an
// alignment with lower penalty to compare greater.
bool operator < ( const flat_alignment& a, const flat_alignment& b )
{ return a.score > b.score ; }

template< typename State > struct enter {
	private: 
		std::deque< State > &o_ ;
	public:
		enter( std::deque< State > &o ) : o_(o) {}
		void operator()( State &s ) {
			greedy( s ) ;
			o_.push_back( s ) ;
			std::push_heap( o_.begin(), o_.end() ) ;
		}
} ;

template< typename State > struct enter_bt {
	private: 
		std::deque< std::pair< State, const State* > > &o_ ;
		const State *s_ ;
	public:
		enter_bt( const State *s, std::deque< std::pair< State, const State* > > &o ) : o_(o), s_(s) {}
		void operator()( State &s ) {
			greedy( s ) ;
			o_.push_back( std::make_pair( s, s_ ) ) ;
			std::push_heap( o_.begin(), o_.end() ) ;
		}
} ;


// Dijkstra's without backtracing
template< typename State >
State find_cheapest( std::deque< State > &open_list )
{
	typename State::ClosedSet closed_list ;
	while( !open_list.empty() )
	{
		State s = open_list.front() ;
		std::pop_heap( open_list.begin(), open_list.end() ) ;
		open_list.pop_back() ;

		if( !is_set_bit( closed_list, s ) ) 
		{
			set_bit( closed_list, s ) ;
			if( finished( s ) ) return s ;
			forward( s, enter<State>( open_list ) ) ;
		}
	}
	throw "ran into a dead end" ;
}

// Dijkstra's with backtracing
template< typename State >
void find_cheapest( std::deque< std::pair< State, const State *> > &open_list )
{
	typename State::ClosedMap closed_list ;
	std::deque< State > used_states ;
	while( !open_list.empty() )
	{
		std::pair< State, const State *> p = open_list.front() ;
		std::pop_heap( open_list.begin(), open_list.end() ) ;
		open_list.pop_back() ;

		if( !is_present( closed_list, p.first ) ) 
		{
			insert( closed_list, p.first, p.second ) ;
			if( finished( p.first ) ) return ; // XXX do the actual backtracing... somehow

			used_states.push_back( p.first ) ;
			forward( p.first, enter_bt<State>( &used_states.back(), open_list ) ) ;
		}
	}
	throw "ran into a dead end" ;
}



#endif
