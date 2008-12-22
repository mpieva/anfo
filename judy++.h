#ifndef JUDYPP_H
#define JUDYPP_H

/*! \page judypp Wrapper around Judy

This is a wrapper around the Judy library.  Judy pretends to be useable
from C++, but its macros actually aren't.  So this wrapper avoids them,
instead making direct use of the ordinary functions.

Due to the way Judy plays games with the representation of pointers,
this class may not be copied.  Wrap it in a smart pointer if you want
to pass it around.

\todo There is no sensible error handling in here, only error codes.
	  Once time permits, this should be enhanced by proper exception
	  handling.

@{ */

#include <Judy.h>
#include <string>

//! \brief Judy Set
//! This class keeps a set of \c Word_t.  Straight-forward wrapper around
//! Judy1 arrays with some cleanup.
class Judy1
{
	private:
		void *arr ;

		Judy1( const Judy1& ) ;            // not implemented
		void operator = ( const Judy1& ) ; // not implemented

	public:
		//! \brief creates a new, empty set.
		Judy1() : arr(0) {}

		//! \brief deallocates before getting rid of a set.
		~Judy1() { clear() ; }

		//! \brief replaces this array with an empty one.
		void clear() { Judy1FreeArray( &arr, 0 ) ; arr = 0 ; }

		//! \brief sets  a bit in this array
		int set( Word_t index ) { return Judy1Set( &arr, index, 0 ) ; }

		//! \brief clears a bit in the array
		int unset( Word_t index ) { return Judy1Unset( &arr, index, 0 ) ; }

		//! \brief finds the first set bit
		//! Finds the first set bit with a key greater or equal to \c
		//! index.  Sets \c index to the found key, else returns an error
		//! code.
		int first( Word_t &index ) const { return Judy1First( arr, &index, 0 ) ; }

		//! \brief finds the next set bit
		//! Finds the first set bit with a key greater than \c index.
		//! Sets \c index to the found key, else returns an error code.
		int next( Word_t &index ) const { return Judy1Next( arr, &index, 0 ) ; }

		//! \brief tells how many bytes are used by this array
		int mem_used() const { return Judy1MemUsed( arr ) ; }

		//! \brief counts set bits in a range
		//! Counts how many bits are set with a key greater or equal to
		//! \c i1 and less than or equal to \c i2.  Default is to count
		//! everything.
		int count( Word_t i1 = 0, Word_t i2 = (Word_t)(-1) ) const { return Judy1Count( arr, i1, i2, 0 ) ; }

		//! \brief tests whether a bit is set.
		int test( Word_t index ) const { return Judy1Test( arr, index, 0 ) ; }
} ;

//! \brief Judy Map
//! This class maps a \c Word_t to an object that fits into a word.  In
//! practice, that means a pointer, an integer or another Judy array.
//! Objects are cleanly constructed and destructed.
template< typename value > class JudyL
{
	private:
		void *arr ;

		JudyL( const JudyL& ) ;            // not implemented
		void operator = ( const JudyL& ) ; // not implemented

	public:
		JudyL() : arr(0) {} 
		~JudyL() { clear() ; }

		//! \brief clears the array
		//! All elements are deleted before deleting the array itself.
		void clear() {
			Word_t w = 0 ;
			for( value *v = first( w ) ; v ; v = next( w ) ) v->~value() ;
			JudyLFreeArray( &arr, 0 ) ;
			arr = 0 ;
		}

		value *insert( Word_t  index ) {
			if( void *p = JudyLGet( arr, index, 0 ) )
				return (value*)p ;
			else return new (JudyLIns( &arr, index, 0 )) value() ;
		}

		value *first(  Word_t& index ) { return (value*)JudyLFirst( arr, &index, 0 ) ; }
		value *next(   Word_t& index ) { return (value*)JudyLNext( arr, &index, 0 ) ; }
		value *last(   Word_t& index ) { return (value*)JudyLLast( arr, &index, 0 ) ; }
		value *get(    Word_t  index ) { return (value*)JudyLGet( arr, index, 0 ) ; }

		int remove( Word_t index ) { if( value *v = get( index ) ) v->~value() ; return JudyLDel( &arr, index, 0 ) ; }
		int mem_used() const { return JudyLMemUsed( arr ) ; }
		int count( Word_t i1 = 0, Word_t i2 = (Word_t)(-1) ) const { return JudyLCount( arr, i1, i2, 0 ) ; }

		//! tests whether the array is empty
		bool empty() { Word_t i = 0 ; return !JudyLFirst( arr, &i, 0 ) ; }
} ;

//! \brief Judy String Map
//! Maps a string to an object that fits into a word, quite similar to
//! JudyL.
template< typename value > struct JudyS
{
	void *arr ;
	size_t maxlen ;

	JudyS() : arr(0), maxlen(0) {} 
	~JudyS() { clear() ; }
	
	void clear() {
		uint8_t s[ maxlen+1 ] ;
		for( value *v = (value*)JudySLFirst( arr, s, 0 ) ; v ;
				v = (value*)JudySLNext( arr, s, 0 ) )
			v->~value() ;

		JudySLFreeArray( &arr, 0 ) ;
		arr = 0 ;
	}

	value *get( const std::string& index ) {
		return (value*)JudySLGet( arr, (uint8_t*)index.c_str(), 0 ) ;
	}

	int remove( const std::string& index ) { 
		uint8_t *i = (uint8_t*)index.c_str() ;
		if( value *v = (value*)JudySLGet( arr, i, 0 ) ) v->~value() ;
		return JudySLDel( &arr, i, 0 ) ;
	}

	value *first( std::string& index )
	{
		char s[ maxlen+1 ] ;
		std::copy( index.begin(), index.end(), s ) ;
		s[ index.size() ] = 0 ;

		value *p = (value*)JudySLFirst( arr, (uint8_t*)s, 0 ) ;
		index = s ;
		return p ;
	}

	value *next( std::string& index )
	{
		char s[ maxlen+1 ] ;
		std::copy( index.begin(), index.end(), s ) ;
		s[ index.size() ] = 0 ;

		value *p = (value*)JudySLNext( arr, (uint8_t*)s, 0 ) ;
		index = s ;
		return p ;
	}

	value *insert( const std::string& index )
	{
		if( maxlen < index.size() ) maxlen = index.size() ;
		uint8_t *s = (uint8_t*)index.c_str() ;
		if( void *p = JudySLGet( arr, s, 0 ) )
			return (value*)p ;
		else
			return new (JudySLIns( &arr, (uint8_t*)index.c_str(), 0 )) value() ; 
	}
} ;

//! @}
#endif
