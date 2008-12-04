#ifndef JUDYPP_H
#define JUDYPP_H

// Wrapper around Judy.h.  This is only useable from C++ because it
// doesn't use the macros.

#include <Judy.h>
#include <string>

class Judy1
{
	private:
		void *arr ;

		Judy1( const Judy1& ) ;            // not implemented
		void operator = ( const Judy1& ) ; // not implemented

	public:
		Judy1() : arr(0) {}
		~Judy1() { clear() ; }

		void clear() { Judy1FreeArray( &arr, 0 ) ; arr = 0 ; }
		int set( Word_t index ) { return Judy1Set( &arr, index, 0 ) ; }
		int unset( Word_t index ) { return Judy1Unset( &arr, index, 0 ) ; }
		int first( Word_t &index ) const { return Judy1First( arr, &index, 0 ) ; }
		int next( Word_t &index ) const { return Judy1Next( arr, &index, 0 ) ; }
		int mem_used() const { return Judy1MemUsed( arr ) ; }
		int count( Word_t i1 = 0, Word_t i2 = (Word_t)(-1) ) const { return Judy1Count( arr, i1, i2, 0 ) ; }
		int test( Word_t index ) const { return Judy1Test( arr, index, 0 ) ; }
} ;

template< typename value > class JudyL
{
	private:
		void *arr ;

		JudyL( const JudyL& ) ;            // not implemented
		void operator = ( const JudyL& ) ; // not implemented

	public:
		JudyL() : arr(0) {} 
		~JudyL() { clear() ; }

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
		bool empty() { Word_t i = 0 ; return !JudyLFirst( arr, &i, 0 ) ; }
} ;

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
#endif
