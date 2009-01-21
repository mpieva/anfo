#ifndef INCLUDED_ALIGN_H
#define INCLUDED_ALIGN_H

#include "index.h"
#include "judy++.h"

#include <deque>

/*!
\page alignment_algorithm Alignment by Dijkstra's Algorithm

\section idea Basic Idea

The idea, in part already realized in \c TWA, is to treat the alignment
problem as finding the shortest path through a generic graph, which in
our case happens to be a regular grid.  Here we have many alignments
starting from different seeds, but we're only interested in the best one
(or maybe a few exceptionally good ones and/or runner-ups on different
chromosomes, species or what have you).  So what we're looking for is
the shortest path from a number of starting places (the seeds) to a
number of goals (any completed semi-global alignment).  This can all be
done in parallel, and the single loop around the algorithm terminates as
soon as at least one alignment is done.  There's no need to ever
consider all the marginal alignments.  (Mega-)blast has to do the
latter, so there's hope we can get a lot faster by avoiding all this.

We choose to encode the complete automaton into the state.  As much as
possible of the automaton should be static information, which is quite
possible with only one type of automaton running at the same time.

To guarantee that we find the best alignment first, all costs (aka
scores) have to be positive (a requirement for anything derived from
Dijkstra's Algorithm).  To put all this on a solid base, the scores are
the negative natural logarithms of probabilities.  Naturally, all scores
are now positive and have a sensible interpretation, though a different
one from the blast scores.  The algorithm will run a bit quicker if we
subtract a small constant from the scores so the minimum cost (for
perfect matches) becomes zero.

\section data_structures Data Structures 

The algorithm keeps a list of \em open \em states (incomplete alignments
that need to be extended) and \em closed \em states (incomplete
alignments that have been extended).  The cheapest \em open \em state is
removed from the queue, made a \em closed \em state, advanced by one
step, and the results become new \em open \em states.  This way,
whenever we close a state, we know that we'll never again visit it at
lower cost.  

The list of potential alignments is a priority queue, obviously sorted
by score.  We will sort alignments of equal score additionally by the
number of remaining nucleotides.  That way, advanced alignments will be
considered first.  The idea here is that we quickly find the best
alignment from the best seed and are done, maybe without even \em
touching the other seeds.  The list of closed states is simply a set.

\todo Consider number of unaligned nucleotides when ordering
      alignments.

\subsection closed_list Closed List

This list may get quite large, so we want a compact representation with
fast lookup.  The "natural" ::std::set uses way too much memory.
Another option would be a simple array, however, since generic pointers
are part of our state (could be changed, but that's cumbersome), an
array would be infeasible.  Also, our state space has a lot of
structure, so instead of O(n<sup>2</sup>) bits, we can reasonably assume to get
away with O(n) bits.  That's why the closed list is a Judy array (nested
Judy arrays, actually).

When supporting backtracing, every state that was closed needs a link to
its parent.  We store that by replacing the Judy1 array with a JudyL
array containing a pointer.  Space usage is bigger now, but we backtrace
only a single alignment.
 
\subsection open_list Open List

We need to be able to extract the cheapest potential alignment, need to
add new ones, and actually we also want to \em overwrite states as soon
as we reach the state on a cheaper route.  The correct data structure,
therefore, is a \em Priority \em Search \em Queue.  Such a beast apparently
doesn't exist in a nice package, so for the time being, we'll just use a
heap and live with the fact that it fill up with additional states that
are never removed because they are simply too bad to ever be touched.

\todo Find or create a priority search queue.  

\todo Include penalty for unaligned tails of the query in score.
*/


/*! \brief Automaton to align according to simple aDNA model
 *
 * This automaton realizes two alignments (to the left and to the
 * right of a seed) in sequence, allowing for affine gap penalties and
 * modelling ancient DNA by separate states for double stranded and
 * single stranded parts.
 * 
 * \todo To be implemented for real.
 */
struct simple_adna {} ;

/*! \brief Extremely simple automaton for testing.
 *
 * This automaton realizes the most trivial alignment possible.  It has
 * no internal states, gap costs are linear, aDNA is not modelled, but
 * the alignment consists of two parts in series.  Matching is also
 * greedy.  The resulting alignment is suitable for testing and also to
 * map sequences from the 454 instrument (due to the very common gaps).
 * \see alignment_rep
 */

//! \page alignment_rep Representation Of Alignments
//!
//! What to store for a (partial) alignment and how to do it?  We'll
//! ultimately need a place to set up on when doing backtracing; this
//! will be within the seed.  We can define that by just two pointers
//! (one reversible, high-resolution DnaP for the reference, one QDnaP
//! with quality scores for the query).  Sequences are self-delimiting,
//! so the two pointers completely define the alignment task---good!.
//! For the ongoing state, we need two offsets (one on each sequence)
//! and a field to store internal state (to keep track of which half
//! we're in, open gaps, aDNA state and possibly more).  We'll use
//! signed shorts for the offsets; �32k should be plenty for the
//! intended applications of a mapper.  Furthermore, the score
//! accumulated so far must be stored, and it is here where we bite off
//! some bits for the state.  That limits the score (cost, actually) to
//! 16M (not really a limit, even if scaled up by ~100, a limit of 160k
//! is far away), and we can store the whole state in 20 bytes in a
//! 32-bit memory model or 24 bytes in a 64 bit memory model.
//!
//! The most simple alignment is ::flat_alignment, see ::simple_adna for
//! a more complex one.

template< typename T > struct gen_alignment {
	DnaP reference ; 						// entry point
	QDnaP query ;							// entry point
	uint32_t state     :  8 ;				// state
	uint32_t penalty   : 24 ;				// cost so far
	int32_t ref_offs   : 16 ;				// cur. position
	int32_t query_offs : 16 ; 				// cur. position

	//! \brief default constructor
	//! Constructs an alignment easily recognizable as invalid.
	gen_alignment() : reference(0), query(0), state(0), penalty(0), ref_offs(0), query_offs(0) {}

	//! \brief prepares an alignment from a seed
	//! Prepares an alignment from a genomic sequence, a sample sequence
	//! and an appropriate seed.  The alignment starts in the middle of
	//! the seed in state 0 at no penalty.
	//!
	//! \param g genome the seeds were prepared from
	//! \param ps sample sequence in compact format
	//! \param s the seed
	gen_alignment( const CompactGenome& g, const QSequence& ps, const Seed& s )
		: reference( (g.get_base() + s.diagonal + s.offset + s.size / 2).reverse_if( s.offset < 0 ) )
		, query( ps.start() + ( s.offset >= 0 ? s.offset + s.size / 2 : -s.offset - s.size/2 ) )
		, state(0), penalty(0), ref_offs(0), query_offs(0) {}


	// associated types: 
	// - set of closed nodes (ClosedSet :: flat_alignment -> Bool)
	// - map of closed nodes to ancestor states
	//   (ClosedMap :: flat_alignment -> Maybe flat_alignment)
#if SMALL_SYS
	typedef JudyL< JudyL< JudyL< JudyL< Judy1 > > > > ClosedSet ;
	typedef JudyL< JudyL< JudyL< JudyL< JudyL< const T* > > > > > ClosedMap ;
#else
	typedef JudyL< JudyL< JudyL< Judy1 > > > ClosedSet ;
	typedef JudyL< JudyL< JudyL< JudyL< const T* > > > > ClosedMap ;
#endif

	//! \brief checks whether this alignment is valid
	//! \return a null pointer iff this is an invalid alignment
	operator const void * () const { return (const void *)reference ; }

	//! \brief returns the current nucleotide on the reference
	Ambicode get_ref() const { return reference[ref_offs] ; }

	//! \brief returns the current nucleotide on the query
	Ambicode get_qry() const { return query[query_offs] ; }

	//! \brief returns the current quality score on the query
	uint8_t  get_qlt() const { return query.qual( query_offs ) ; }
} ;

struct flat_alignment : public gen_alignment<flat_alignment> {
	//! \brief Trivial substitution matrix.
	//! Aligning two codes that have any overlap costs nothing, else it
	//! costs 1.
	static uint32_t subst_mat( Ambicode a, Ambicode b ) { return (a & b) != 0 ? 0 : 1 ; }

	//! \brief trivial gap costs
	//!A gap costs one, whether it is opened or extended.
	static uint32_t gap_penalty() { return 1 ; }

	flat_alignment() : gen_alignment<flat_alignment>() {}
	flat_alignment( const CompactGenome& g, const QSequence& ps, const Seed& s ) : gen_alignment<flat_alignment>(g,ps,s) {}
} ;

//! formats an intermediate alignment state to a stream.
//! This is useful mostly for debugging.
//! \internal 
template< typename A > std::ostream& operator << ( std::ostream& s, const gen_alignment<A>& fa )
{
	return s << fa.reference + fa.ref_offs << " x "
		     << fa.query + fa.query_offs << " @"
			 << fa.state << ": " << fa.penalty ;
}

//! \brief removes cheapest alignment from list
//! Behaviour on an empty list is undefined.
//! \param ol open list
template< typename A > void pop_top( std::deque< A >& ol ) { std::pop_heap( ol.begin(), ol.end() ) ; ol.pop_back() ; }

//! \page generic_alignment_ops Basic Operations over Alignment States
//! The functions in here are still reasonably general and might be
//! useful for a variety of concrete alignments.  Anything not sensibly
//! derivable from ::gen_alignment should specialze them.
//!
//! \page trivial_alignments Trivial Alignments
//! \subpage basic_alignment_ops
//! @{

//! \brief cheks whether an alignment is already closed
//! \param cl closed set
//! \param s alignment state
//! \return true iff \c s in contained in \c cl
template< typename A > bool is_present( const typename gen_alignment<A>::ClosedSet& cl, const gen_alignment<A>& s )
{
	return cl, s.reference.get(), s.query.get(),
#if SMALL_SYS
			   s.reference.high() << 16 | s.query.high(),
#endif
			   s.ref_offs, s.query_offs | s.state << 16 ;
}

template< typename A > const_ref<const A*> lookup(
		const typename gen_alignment<A>::ClosedMap& cl, const gen_alignment<A>& s )
{
	return cl, s.reference.get(), s.query.get(),
#if SMALL_SYS
			   s.reference.high() << 16 | s.query.high(),
#endif
			   s.ref_offs, s.query_offs | s.state << 16 ;
}

template< typename A > void set_bit( typename gen_alignment<A>::ClosedSet& cl, const gen_alignment<A>& s )
{
	cl.insert( s.reference.get() )
	  ->insert( s.query.get() )
#if SMALL_SYS
	  ->insert( s.reference.high() << 16 | s.query.high() )
#endif
	  ->insert( s.ref_offs )
	  ->set( s.query_offs | s.state << 16 ) ;
}

template< typename A > void insert(
		typename gen_alignment<A>::ClosedMap& cl, const gen_alignment<A>& s, const A *p )
{ 
	*cl.insert( s.reference.get() )
	   ->insert( s.query.get() )
#if SMALL_SYS
       ->insert( s.reference.high() << 16 | s.query.high() )
#endif
	   ->insert( s.ref_offs )
	   ->insert( s.query_offs | s.state << 16 ) = p ;
}

//! \brief resets alignment to initial state.
//! The internal state is reset, the penalty set to zero and the
//! reference and query positions to their initial values.  The alignment
//! can be repeated now.
template< typename A > void reset( gen_alignment<A>& fa )
{
	fa.ref_offs = 0 ;
	fa.query_offs = 0 ;
	fa.state = 0 ;
	fa.penalty = 0 ;
}

/*! \brief compares two simple alignments.
 *
 * Since we are building a heap and the STL heap puts the \em largest
 * element at the top, we want an alignment with lower penalty to
 * compare greater.
 */
template< typename A >
bool operator < ( const gen_alignment<A>& a, const gen_alignment<A>& b )
{ return b.penalty < a.penalty ; }

//! }@

//! \page special_align_ops Specific Alignment Operations
//! These functions should be specialized for every concrete type of
//! alignment.
//!
//! {@

//! \brief checks whether an alignment is finished
//! A \c flat_alignment is finished iff either sequence hit a gap while
//! doing the second half of an alignment.
//! \param s alignment state
//! \return true iff the alignment is finished
bool finished( const flat_alignment& s ) {
	return s.state == 1 &&
		( s.reference[ s.ref_offs ] == 0 || s.query[ s.query_offs ] == 0 ) ;
}

//! \brief greedily extends an alignment.
//! This function simply extends an alignment as long at it finds
//! matches.  It is important that this is done in a way that leaves
//! enough traces if backtracing is desired.  Here we do this by never
//! changing state in greedy(), but only in forward().  greedy() just
//! advances two pointers synchronously.
//! \param s alignment state
void greedy( flat_alignment& s )
{
	while( s.query[s.query_offs] && s.reference[s.ref_offs] == s.query[s.query_offs] )
	{
		s.penalty += flat_alignment::subst_mat( s.reference[s.ref_offs], s.query[s.query_offs] ) ;
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

//! \todo this is all very ugly; need to abstract the gunk out of here.
template< typename F > void forward( const flat_alignment& s, F f )
{
	// what to do?  
	// in state 0: mismatch, gap ref, gap query
	// in state 1: mismatch, gap ref, gap query
	// In any case, we can rely on the fact that greedy() was called
	// before.
	if( s.state == 0 )
	{
		if( s.reference[s.ref_offs] == 0 || s.query[s.query_offs] == 0 )
		{
			flat_alignment s1 = s ;
			s1.state = 1 ;
			s1.ref_offs = -1 ;
			s1.query_offs = -1 ;
			f( s1 ) ;
		}
		else
		{
			flat_alignment s1 = s ;
			s1.penalty += flat_alignment::subst_mat( s1.reference[s1.ref_offs], s1.query[s1.query_offs] ) ;
			++s1.ref_offs ;
			++s1.query_offs ;
			f( s1 ) ;
		
			flat_alignment s2 = s ;
			s2.penalty += flat_alignment::gap_penalty() ;
			++s2.ref_offs ;
			f( s2 ) ;
	
			flat_alignment s3 = s ;
			s3.penalty += flat_alignment::gap_penalty() ;
			++s3.query_offs ;
			f( s3 ) ;
		}
	}
	else
	{
		{
			flat_alignment s1 = s ;
			s1.penalty += flat_alignment::subst_mat( s1.reference[s1.ref_offs], s1.query[s1.query_offs] ) ;
			--s1.ref_offs ;
			--s1.query_offs ;
			f( s1 ) ;
		}{
			flat_alignment s2 = s ;
			s2.penalty += flat_alignment::gap_penalty() ;
			--s2.ref_offs ;
			f( s2 ) ;
		}{
			flat_alignment s3 = s ;
			s3.penalty += flat_alignment::gap_penalty() ;
			--s3.query_offs ;
			f( s3 ) ;
		}
	}
}

//! @}

template< typename State > struct enter {
	private: 
		std::deque< State > &o_ ;
		typename State::ClosedSet *c_ ;
		uint32_t mp_ ;
	public:
		enter(
				std::deque< State > &o,
				uint32_t mp = std::numeric_limits<uint32_t>::max(),
				typename State::ClosedSet *c = 0
			 ) : o_(o), c_(c), mp_(mp) {}
		void operator()( State s ) {
			greedy( s ) ;
			if( (!c_ || !is_present( *c_, s )) && s.penalty <= mp_ )
			{
				o_.push_back( s ) ;
				std::push_heap( o_.begin(), o_.end() ) ;
			}
		}
} ;

template< typename State > struct enter_bt {
	private: 
		std::deque< std::pair< State, const State* > > &o_ ;
		typename State::ClosedMap *c_ ;
		const State *s_ ;
		uint32_t mp_ ;
	public:
		enter_bt(
				std::deque< std::pair< State, const State* > > &o,
				uint32_t mp = std::numeric_limits<uint32_t>::max(),
				typename State::ClosedMap *c = 0,
				const State *s = 0
				) : o_(o), c_(c), s_(s), mp_(mp) {}
		void operator()( State s ) {
			greedy( s ) ;
			if( (!c_ || !lookup( *c_, s )) && s.penalty <= mp_ )
			{
				o_.push_back( std::make_pair( s, s_ ) ) ;
				std::push_heap( o_.begin(), o_.end() ) ;
			}
		}
} ;


/*! \brief Dijkstra's without backtracing.
 * If called with an appropriately prepared \c open_list, this function
 * runs Disjkstra's algorithm by repeatedly calling \c forward and \c
 * greedy on the states in the \c open_list.  It returns the first state
 * to reach its goal according to finished().  The algorithm is
 * configured by using different types in place of \c State.
 *
 * \param open_list Queue of open states, must form a heap and \c greedy
 *                  must have been called on each of its elements.
 * \param max_penalty Penalty at which alignments are no longer
 *                    interesting.
 * \param report_success if set, a statistic of the finished alignment
 *                       is sent to std::clog (only useful for
 *                       debugging)
 * \return First state to be detected as finished().
 */
template< typename State >
State find_cheapest( 
		std::deque< State > &open_list, 
		uint32_t max_penalty = std::numeric_limits<uint32_t>::max(),
		bool report_success = false )
{
	typename State::ClosedSet closed_list ;
	while( !open_list.empty() )
	{
		State s = open_list.front() ;
		std::pop_heap( open_list.begin(), open_list.end() ) ;
		open_list.pop_back() ;

		if( !is_present( closed_list, s ) ) 
		{
			set_bit( closed_list, s ) ;
			if( finished( s ) )
			{
				int dups = 0 ;
				for( size_t i = 0 ; i != open_list.size() ; ++i )
					if( is_present( closed_list, open_list[i] ) ) ++dups ;

				if( report_success ) std::clog
						<< "Finished alignment, open list contains "
						<< open_list.size() << " nodes, " << deep_count( closed_list )
						<< " nodes are closed, " << dups << " of which are still tracked."
						<< std::endl ;
				
				return s ;
			}
			forward( s, enter<State>( open_list, max_penalty, &closed_list ) ) ;
		}
	}
	return State() ;
}

typedef std::deque< std::pair< Ambicode, Ambicode > > Trace_ ;
struct Trace
{
	Trace_ trace ;
	DnaP minpos, maxpos ;
} ;

//! prints a backtrace in three-line format.
//! This is intended for debugging, it prints a backtraced alignment in
//! two lines of sequence and one "conservation" line.
//! \internal
inline std::ostream& operator << ( std::ostream& s, const Trace_& t )
{
	for( Trace_::const_iterator i = t.begin(), e = t.end() ; i != e ; ++i )
		s << from_ambicode( i->first ) ; 
	s << '\n' ;
	for( Trace_::const_iterator i = t.begin(), e = t.end() ; i != e ; ++i )
		s << from_ambicode( i->second ) ; 
	s << '\n' ;
	for( Trace_::const_iterator i = t.begin(), e = t.end() ; i != e ; ++i )
		s << ( i->first == i->second ? '*' : i->first & i->second ? '.' : ' ') ;
	s << '\n' ;
	return s ;
}

//! \brief backtraces a simple alignment returning sequences
//!
//! This backtraces an alignment after it has been created by
//! find_cheapest() called with a backtracing structure.  Backtracing
//! works by looking at two intermediate states, the \em current one and
//! its \em predecessor.  In between those two states, exactly one
//! invocation of forward() and one of greedy() have happened.  To
//! backtrace, we first check in which direction we moved (depends on
//! which half of the alignment we were in).  Next, if both the reference
//! and the query offsets differ between states, we copy symbols from
//! both (this corresponds to the greedy extension or a mismatch).  If
//! only one differs, we copy one symbol and fill it up with a gap.
//! Depending on the direction we moved in, the pair is added at the
//! front or the end of the trace.  If the internal state changed, we
//! don't trace at all (we just jump).
//!
//! \param cl the ClosedMap that was used in find_cheapest()
//! \param a pointer to the alignment that needs to be traced
//! \return a trace, that is a sequence of pairs of Ambicodes
//! \internal

Trace backtrace( const flat_alignment::ClosedMap &cl, const flat_alignment *a )
{
	Trace r ;

	// When the alignment finished, it pointed to the minimum
	// coordinate (to a gap actually).  Just store it.
	r.minpos = a->reference + a->ref_offs ;

	// Only trace back second state here, this ends up at the front of
	// the alignment, but we add stuff to the back as we generate it in
	// the wrong order.
	while( const flat_alignment *b = *lookup( cl, *a ) ) 
	{
		if( b->state == 0 ) break ;
		for( signed short ro = a->ref_offs, qo = a->query_offs ; 
				ro != b->ref_offs || qo != b->query_offs ; )
		{
			Ambicode x = ro == b->ref_offs   ? 0 : a->reference[++ro] ;
			Ambicode y = qo == b->query_offs ? 0 : a->query[++qo] ;
			r.trace.push_back( std::make_pair( x,y ) ) ;
		}
		a = b ;
	}

	// Trace further to initiation of second state
	for( signed short ro = a->ref_offs, qo = a->query_offs ; ro != -1 || qo != -1 ; )
	{
		Ambicode x = ro ? a->reference[++ro] : 0 ;
		Ambicode y = qo ? a->query[++qo] : 0 ;
		r.trace.push_back( std::make_pair( x,y ) ) ;
	}

	// Skip one state, this is the one aligning the terminal gap.
	a = *lookup( cl, *a ) ; 
		
	// Now at the end of the first phase, we got a pointer to the
	// maximum coordinate (again a gap).  Store it.
	r.maxpos = a->reference + a->ref_offs ;

	// Trace back first state now, generating a new trace which needs to
	// be reversed in the end.
	Trace_ t2 ;
	while( const flat_alignment *b = *lookup( cl, *a ) ) 
	{
		for( signed short ro = a->ref_offs, qo = a->query_offs ; 
				ro != b->ref_offs || qo != b->query_offs ; )
		{
			Ambicode x = ro == b->ref_offs   ? 0 : a->reference[--ro] ;
			Ambicode y = qo == b->query_offs ? 0 : a->query[--qo] ;
			t2.push_back( std::make_pair( x,y ) ) ;
		}
		a = b ;
	}

	// Now (*b) is null, (*a) is the last state ever generated.  Now
	// trace further until we hit the initial state (both offsets
	// vanish).  We can again add this to t2, as we're going in the same
	// direction.
	for( signed short ro = a->ref_offs, qo = a->query_offs ; ro || qo ; )
	{
		Ambicode x = ro ? a->reference[--ro] : 0 ;
		Ambicode y = qo ? a->query[--qo] : 0 ;
		t2.push_back( std::make_pair( x,y ) ) ;
	}

	// To see the crack between the two halves, use this:
	// t2.push_back( std::make_pair( 0,0 ) ) ;

	r.trace.insert( r.trace.end(), t2.rbegin(), t2.rend() ) ;
	return r ;
}

/*! \brief Dijkstra's with backtracing.
 *
 * See \c find_cheapest, but this implementation also does backtracing
 * (at higher memory cost, naturally).  
 */
template< typename State >
Trace find_cheapest(
		std::deque< std::pair< State, const State *> > &open_list,
		uint32_t max_penalty = std::numeric_limits<uint32_t>::max() ) 
{
	typename State::ClosedMap closed_list ;
	std::deque< State > used_states ;
	while( !open_list.empty() )
	{
		std::pair< State, const State *> p = open_list.front() ;
		pop_top( open_list ) ;

		if( !lookup( closed_list, p.first ) ) 
		{
			insert( closed_list, p.first, p.second ) ;
			if( finished( p.first ) ) return backtrace( closed_list, &p.first ) ;
			used_states.push_back( p.first ) ;
			forward( p.first, enter_bt<State>( open_list, max_penalty, &closed_list, &used_states.back() ) ) ;
		}
	}
	return Trace() ;
}

//! \brief initializes alignments from a list of seeds
//!
//! Each seed, no matter its quality, gives rise to a potential
//! alignment.  A heap of alignments is updated; on return it can
//! immediately be fed to find_cheapest().
// 
//! \param g Genome the seeds came from.
//! \param ps Sequence to map in compact format.
//! \param begin Iterator to first seed.
//! \param end Iterator to last seed.
//! \param ol \c Open \c List to be updated.  Must form a heap (an empty
//!           heap is fine).

template< typename Aln, typename Iter >
void setup_alignments( const CompactGenome& g, const QSequence& ps,
		Iter begin, const Iter& end, std::deque< Aln >& ol )
{
	for( ; begin != end ; ++begin ) (enter<Aln>( ol ))( Aln( g, ps, *begin ) ) ;
}


#endif
