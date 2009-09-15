#ifndef INCLUDED_ALIGN_H
#define INCLUDED_ALIGN_H

#include "index.h"
#include "judy++.h"
#include "logdom.h"
#include "stream.h"

#include <output.pb.h>

#include <cmath>
#include <deque>
#include <sstream>

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
//! signed shorts for the offsets; ±32k should be plenty for the
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
	const QSequence::Base *query ;						// entry point
	uint32_t state     :  8 ;				// state
	uint32_t penalty   : 24 ;				// cost so far
	int32_t ref_offs   : 16 ;				// cur. position
	int32_t query_offs : 16 ; 				// cur. position

	enum { mask_dir = 1 } ;

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
	typedef JudyL< JudyL< JudyL< Judy1 > > > ClosedSet ;
	typedef JudyL< JudyL< JudyL< JudyL< const T* > > > > ClosedMap ;

	//! \brief checks whether this alignment is valid
	//! \return a null pointer iff this is an invalid alignment
	operator const void * () const { return (const void *)reference ; }

	//! \brief returns the current nucleotide on the reference
	Ambicode get_ref() const { return reference[ref_offs] ; }

	//! \brief returns the current nucleotide on the query
	QSequence::Base get_qry() const { return query[query_offs] ; }

	//! \brief returns the current quality score on the query
	// uint8_t  get_qlt() const { return query[query_offs].qualities[.qual( query_offs ) ; }

	void adv_ref() { if( state & mask_dir ) --ref_offs ; else ++ref_offs ; }
	void adv_qry() { if( state & mask_dir ) --query_offs ; else ++query_offs ; }
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
//! \subpage generic_alignment_ops
//! @{

//! \brief cheks whether an alignment is already closed
//! \param cl closed set
//! \param s alignment state
//! \return true iff \c s in contained in \c cl
template< typename A > bool is_present( const typename gen_alignment<A>::ClosedSet& cl, const gen_alignment<A>& s )
{
	return cl, s.reference.get(), (Word_t)s.query,
			   s.reference.high() << 16 | s.ref_offs,
			   s.query_offs | s.state << 16 ;
}

template< typename A > const_ref<const A*> lookup(
		const typename gen_alignment<A>::ClosedMap& cl, const gen_alignment<A>& s )
{
	return cl, s.reference.get(), (Word_t)s.query,
			   s.reference.high() << 16 | s.ref_offs,
			   s.query_offs | s.state << 16 ;
}

template< typename A > void set_bit( typename gen_alignment<A>::ClosedSet& cl, const gen_alignment<A>& s )
{
	cl.insert( s.reference.get() )
	  ->insert( (Word_t)s.query )
	  ->insert( s.reference.high() << 16 | s.ref_offs )
	  ->set( s.query_offs | s.state << 16 ) ;
}

template< typename A > void insert(
		typename gen_alignment<A>::ClosedMap& cl, const gen_alignment<A>& s, const A *p )
{ 
	*cl.insert( s.reference.get() )
	   ->insert( (Word_t)s.query )
       ->insert( s.reference.high() << 16 | s.ref_offs )
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


//! \page trivial_alignment Operations For Trivial Alignments
//! These functions should be specialized for every concrete type of
//! alignment.
//!
//! @{

//! \brief Extremely simple automaton for testing.
//! This automaton realizes the most trivial alignment possible.  It has
//! no internal states, gap costs are linear, aDNA is not modelled, but
//! the alignment consists of two parts in series.  Matching is also
//! greedy.  The resulting alignment is suitable for testing and also to
//! map sequences from the 454 instrument (due to the very common gaps).
//! \see alignment_rep

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

//! \brief checks whether an alignment is finished
//! The typical alignment is finished iff the query sequence hit a gap
//! while doing the second half of an alignment.
//! \param s alignment state
//! \return true iff the alignment is finished
template< typename A > bool finished( const gen_alignment<A>& s )
{
	return (s.state & gen_alignment<A>::mask_dir) &&
		s.get_qry().ambicode == 0 ;
}

//! \brief greedily extends an alignment.
//! This function simply extends an alignment as long at it finds
//! matches.  It is important that this is done in a way that leaves
//! enough traces if backtracing is desired.  Here we do this by never
//! changing state in greedy(), but only in forward().  greedy() just
//! advances two pointers synchronously.
//! \param s alignment state
inline void greedy( flat_alignment& s )
{
	while( s.get_qry().ambicode && s.get_ref() == s.get_qry().ambicode )
	{
		s.penalty += flat_alignment::subst_mat( s.get_ref(), s.get_qry().ambicode ) ;
		s.adv_ref() ;
		s.adv_qry() ;
	}
}

template< typename F > void forward( const flat_alignment& s, F f )
{
	// what to do?  
	// in state 0: mismatch, gap ref, gap query
	// in state 1: mismatch, gap ref, gap query
	// In any case, we can rely on the fact that greedy() was called
	// before.
	if( (s.state & flat_alignment::mask_dir) == 0 &&
			(s.get_ref() == 0 || s.get_qry().ambicode == 0) )
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
		s1.penalty += flat_alignment::subst_mat( s1.get_ref(), s1.get_qry().ambicode ) ;
		s1.adv_ref() ;
		s1.adv_qry() ;
		f( s1 ) ;

		flat_alignment s2 = s ;
		s2.penalty += flat_alignment::gap_penalty() ;
		s2.adv_ref() ;
		f( s2 ) ;

		flat_alignment s3 = s ;
		s3.penalty += flat_alignment::gap_penalty() ;
		s3.adv_qry() ;
		f( s3 ) ;
	}
}

//! @}

//! \page adna_alignment Operations on Ancient DNA Alignments
//! This type of alignment models aDNA as two additional states, one
//! with higher deamination rate, with an asymmetric substitution
//! matrix, and with affine gap costs.  Additional parameters are all
//! static.
//!
//! @{

//! \brief substitution matrix
//! We prepare a full matrix of 16 ambiguity codes vs. 16 ambiguity
//! codes, this avoid the need for expensive additions in the log-domain
//! should the need to align ambiguity codes arise.  First index is
//! "from" (reference code), second index is "to" (query code).
typedef Logdom subst_mat[16][16] ;

//! \brief aDNA alignment automaton
//! We encode the state as follows: Bit 0 is set if we're in the second
//! half (5' end), bit 1 is set if we're threating DNA as single
//! stranded, bit 2 is set if we're gapping the reference and bit 3 is
//! set if we're gapping the query.
//! \see alignment_rep

struct simple_adna : public gen_alignment<simple_adna> {
	//! \brief DS substitution matrix, forward direction
	static subst_mat ds_mat ;

	//! \brief SS substitution matrix, forward direction
	//! Deamination shows up as C->T as it is best understood this way.
	//! To process reverse-complemented deamination, we have to do
	//! rev-complemented lookups while actually moving in the forward
	//! (5'->3') direction.  \see simple_adna::subst_penalty()
	static subst_mat ss_mat ;

	//! \brief does a lookup in the appropriate subst matrix
	//! This retrieves the current codes in reference and query, if
	//! necessary complements them, then retrieves the appropriate score
	//! from the the currently active scoring matrix.  Note that
	//! complementing is necessary as long as the mask_dir bit is \e not
	//! set, since in the forward direction, aDNA damage actually
	//! operates on the reverse strand.
	Logdom subst_penalty() const {
		Ambicode r = state & mask_dir ? get_ref() : complement( get_ref() ) ;
		if( !r ) r = 15 ; // if reference has a gap, pretend it was an N

		Logdom prob = Logdom::null() ;
		for( uint8_t p = 0 ; p != 4 ; ++p )
		{
			Ambicode q = state & mask_dir ? (1<<p) : complement(1<<p) ;
			prob += ( state & mask_ss ? ss_mat[r][q] : ds_mat[r][q] )
				* Logdom::from_phred( get_qry().qscores[p] ) ;
		}
		return prob ;
	}

	//! \brief Penalty for extending an overhang.
	//! Having a constant penalty for the overhang length models its
	//! length distribution as geometric.
	static Logdom overhang_ext_penalty ;

	//! \brief penalty for entering SS state
	//! This is essentially the probability of having an overhang at
	//! all.
	static Logdom overhang_enter_penalty ;

	//! \brief gap open penalty
	static Logdom gap_open_penalty ;

	//! \brief gap extension penalty
	static Logdom gap_ext_penalty ;

	//! \brief sets up parameters from configuration block
	static void configure( const config::Aligner&, std::ostream* = 0 ) ;

	simple_adna() : gen_alignment<simple_adna>() {}
	simple_adna( const CompactGenome& g, const QSequence& ps, const Seed& s ) : gen_alignment<simple_adna>(g,ps,s) {}

	enum {
		mask_gap_ref = 2,
		mask_gap_qry = 4,
		mask_ss      = 8,

		mask_gaps = mask_gap_ref | mask_gap_ref
	} ;
} ;

//! \brief greedily extends an alignment.
//! We extend greedily if no gap is open and the two sequences match.
//! \see ::greedy( flat_alignment& )
//! \param s alignment state
inline void greedy( simple_adna& s )
{
	if( (s.state & simple_adna::mask_gaps) == 0 )
	{
		while( s.get_qry().ambicode && s.get_ref() == s.get_qry().ambicode )
		{
			s.penalty += s.subst_penalty().to_phred() ;
			if( s.state & simple_adna::mask_ss ) s.penalty += simple_adna::overhang_ext_penalty.to_phred() ;
			s.adv_ref() ;
			s.adv_qry() ;
		}
	}
}

template< typename F > void forward( const simple_adna& s, F f )
{
	// what to do?  
	// If in matching state, we know there's no immediate match, so we can...
	// - mismatch
	// - detect deamination and change to SS state (while matching)
	// - open ref gap
	// - open query gap
	//
	// If a gap is open, we can...
	// - extend it
	// - close it
	//
	// If we hit a gap symbol, we must...
	// - start over at second half in initial state

	// Note the penalties: The appropriate substitution penalty is
	// applied whenever we (mis-)match two codes, the gap open penalties
	// are applied when opening/extending a gap, the
	// overhang_enter_penalty is applied when changing to SS mode and
	// the overhang_ext_penalty is applied whenever moving along the
	// query while single stranded, even when a gap is open!  This gives
	// correct scores for a geometric distribution of overhang lengths.
	if( !s.get_ref() && s.get_qry().ambicode )
	{
		// We hit a gap in the reference, whatever is left of the query
		// must be penalized.  To do this, we virtually extend the
		// reference with Ns and align to those.  This is a white lie
		// in that it wil overestimate the real penalty, but that's
		// okay, because such an alignment isn't all that interesting
		// in reality anyway.  It woould feel more natural to do this
		// in greedy(), but we do it here to leave sufficient traces
		// for correct backtracing.
		simple_adna s1 = s ;
		while( s1.get_qry().ambicode )
		{
			s1.penalty += s1.subst_penalty().to_phred() ;
			if( s1.state & simple_adna::mask_ss ) s1.penalty += simple_adna::overhang_ext_penalty.to_phred() ;
			s1.adv_qry() ;
			f( s1 ) ;
		}
	}
	else if( (s.state & simple_adna::mask_dir) == 0 && !s.get_qry().ambicode )
	{
		// forward dir, hit gap in query --> start over in reverse dir
		simple_adna s1 = s ;
		s1.state = 1 ;
		s1.ref_offs = -1 ;
		s1.query_offs = -1 ;
		f( s1 ) ;
	}
	else if( (s.state & simple_adna::mask_gaps) == 0 )
	{
		// no gaps open --> mismatch, open either gap, enter SS
		{
			simple_adna s1 = s ;
			s1.penalty += s1.subst_penalty().to_phred() ;
			if( s.state & simple_adna::mask_ss ) s1.penalty += simple_adna::overhang_ext_penalty.to_phred() ;
			s1.adv_ref() ;
			s1.adv_qry() ;
			f( s1 ) ;
		}{
			simple_adna s2 = s ;
			s2.penalty += simple_adna::gap_open_penalty.to_phred() ;
			s2.state |= simple_adna::mask_gap_qry ;
			s2.adv_ref() ;
			f( s2 ) ;
		}{
			simple_adna s3 = s ;
			s3.penalty += simple_adna::gap_open_penalty.to_phred() ;
			if( s.state & simple_adna::mask_ss ) s3.penalty += simple_adna::overhang_ext_penalty.to_phred() ;
			s3.state |= simple_adna::mask_gap_ref ;
			s3.adv_qry() ;
			f( s3 ) ;
		}
		if( simple_adna::overhang_enter_penalty.is_finite() && (s.state & simple_adna::mask_ss) == 0 )
		{
			// To enter single stranded we require that the penalty for
			// doing so is immediately recovered by the better match.
			// This is easily the case for the observed deamination
			// rates in aDNA.
			simple_adna s4 = s ;
			s4.state |= simple_adna::mask_ss ;
			uint32_t p4 = ( s4.subst_penalty() + simple_adna::overhang_enter_penalty 
			                                   + simple_adna::overhang_ext_penalty ).to_phred() ;
			uint32_t p0 = s.subst_penalty().to_phred() ;
			if( p4 < p0 ) {
				s4.penalty += p4 ;
				s4.adv_ref() ;
				s4.adv_qry() ;
				f( s4 ) ;
			}
		}
	}
	else 
	{
		// already gapping (ref or qry) --> continue or close
		{
			simple_adna s1 = s ;
			s1.penalty += simple_adna::gap_ext_penalty.to_phred() ;
			bool which = (s.state & simple_adna::mask_gaps) == simple_adna::mask_gap_ref ;
			if( which && (s.state & simple_adna::mask_ss) ) s1.penalty += simple_adna::overhang_ext_penalty.to_phred() ;
			if( which ) s1.adv_qry() ; else s1.adv_ref() ;
			f( s1 ) ;
		}{
			simple_adna s2 = s ;
			s2.state &= ~simple_adna::mask_gaps ;
			s2.penalty += s2.subst_penalty().to_phred() ;
			if( s.state & simple_adna::mask_ss ) s2.penalty += simple_adna::overhang_ext_penalty.to_phred() ;
			s2.adv_ref() ;
			s2.adv_qry() ;
			f( s2 ) ;
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
 * \param closed_list Set of closed nodes, initially empty.
 * \param max_penalty Penalty at which alignments are no longer
 *                    interesting.
 * \param open_nodes_after_alignment
 *      if not 0, will contain the number of nodes left in the open list
 *      after aligning (useful for debugging and tuning only)
 * \param closed_nodes_after_alignment
 *      if not 0, will contain the number of nodes in the closed list
 *      aftert alignment (useful for debugging and tuning only)
 * \param tracked_closed_nodes_after_alignment
 *      if not 0, will contain the number of nodes in the open list that
 *      are already closed (useful for debugging and tuning only)
 * \return First state to be detected as finished() or an invalid state
 *         if no good alignment could be found.
 */
template< typename State >
State find_cheapest( 
		std::deque< State > &open_list, 
		typename State::ClosedSet &closed_list,
		uint32_t max_penalty = std::numeric_limits<uint32_t>::max(),
		uint32_t *open_nodes_after_alignment = 0,
		uint32_t *closed_nodes_after_alignment = 0,
		uint32_t *tracked_closed_nodes_after_alignment = 0 )
{
	while( !open_list.empty() && !exit_with )
	{
		State s = open_list.front() ;
		std::pop_heap( open_list.begin(), open_list.end() ) ;
		open_list.pop_back() ;

		if( !is_present( closed_list, s ) ) 
		{
			set_bit( closed_list, s ) ;
			if( finished( s ) )
			{
				if( open_nodes_after_alignment ) *open_nodes_after_alignment = open_list.size() ;
				if( closed_nodes_after_alignment ) *closed_nodes_after_alignment = deep_count( closed_list ) ;
				if( tracked_closed_nodes_after_alignment ) {
					int dups = 0 ;
					for( size_t i = 0 ; i != open_list.size() ; ++i )
						if( is_present( closed_list, open_list[i] ) ) ++dups ;
					*tracked_closed_nodes_after_alignment = dups ;
				}
				return s ;
			}
			forward( s, enter<State>( open_list, max_penalty, &closed_list ) ) ;
		}
	}

    if( open_nodes_after_alignment ) *open_nodes_after_alignment = 0 ;
    if( closed_nodes_after_alignment ) *closed_nodes_after_alignment = deep_count( closed_list ) ;
    if( tracked_closed_nodes_after_alignment ) *tracked_closed_nodes_after_alignment = 0 ;
	return State() ;
}


//! \brief backtraces an alignment and return a CIGAR line
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
//! don't trace at all (we just jump).  See output.proto for the
//! encoding of these CIGAR lines.
//!
//! \param cl the ClosedMap that was used in find_cheapest()
//! \param a final state to start backtracing from
//! \param minpos will be filled by position of smaller end of alignment
//! \param maxpos will be filled by position of greater end of alignment
//! \return binary CIGAR string
//! \internal
//! \todo Write mismatches with different code (to allow various
//!       calculations without the genome being available).

template< typename State > std::vector<unsigned>
backtrace( const typename State::ClosedMap &cl, const State *a, DnaP &minpos, DnaP &maxpos )
{
	std::vector<unsigned> fwd, rev ;

	// When the alignment finished, it pointed to the minimum
	// coordinate (to a gap actually).  Just store it.
	minpos = a->reference + a->ref_offs ;

	// Only trace back second state here, this ends up at the front of
	// the alignment, but we add stuff to the back as we generate it in
	// the wrong order.
	while( const State *b = *lookup( cl, *a ) ) 
	{
		if( (b->state & simple_adna::mask_dir) == 0 ) break ;
		int dr = b->ref_offs - a->ref_offs ;
		int dq = b->query_offs - a-> query_offs ;
		if( int m = std::min( dr, dq ) ) {
			dr -= m ;
			dq -= m ;
			streams::push_m( fwd, m ) ;
		}
		streams::push_i( fwd, dq ) ;
		streams::push_d( fwd, dr ) ;
		a = b ;
	}

	assert( a->ref_offs == a->query_offs ) ;
	streams::push_m( fwd, -a->ref_offs-1 ) ;
	fwd.push_back( 0 ) ;

	// Skip one state, this is the one aligning the terminal gap.
	a = *lookup( cl, *a ) ; 
		
	// Now at the end of the first phase, we got a pointer to the
	// maximum coordinate (again a gap).  Store it.
	maxpos = a->reference + a->ref_offs ;

	// Trace back first state now, generating a new trace which needs to
	// be reversed in the end.
	while( const State *b = *lookup( cl, *a ) ) 
	{
		int dr = a->ref_offs - b->ref_offs ;
		int dq = a->query_offs - b-> query_offs ;
		if( int m = std::min( dr, dq ) ) {
			dr -= m ;
			dq -= m ;
			streams::push_m( rev, m ) ; 
		}
		streams::push_i( rev, dq ) ;
		streams::push_d( rev, dr ) ;
		a = b ;
	}

	// Now (*b) is null, (*a) is the last state ever generated.  Now
	// trace further until we hit the initial state (both offsets
	// vanish).  We can again add this to t2, as we're going in the same
	// direction.
	assert( a->ref_offs == a->query_offs ) ;
	streams::push_m( rev, a->ref_offs ) ;

	fwd.insert( fwd.end(), rev.rbegin(), rev.rend() ) ;
	return fwd ;
}

/*! \brief Dijkstra's with backtracing.
 *
 * See \c find_cheapest, but this implementation also does backtracing
 * (at higher memory cost, naturally).  
 */
template< typename State >
std::vector<unsigned> find_cheapest(
		std::deque< std::pair< State, const State *> > &open_list,
		DnaP &minpos, DnaP &maxpos,
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
			if( finished( p.first ) ) return backtrace( closed_list, &p.first, minpos, maxpos ) ;

			used_states.push_back( p.first ) ;
			forward( p.first, enter_bt<State>( open_list, max_penalty, &closed_list, &used_states.back() ) ) ;
		}
	}
	return std::vector<unsigned>() ;
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
