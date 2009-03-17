#ifndef INCLUDED_COMPRESS_STREAM_H
#define INCLUDED_COMPRESS_STREAM_H

#include "config.h"
#include "util.h"

#include <google/protobuf/io/zero_copy_stream.h>
#include <memory>
#include <ostream>

#include <iostream>

//! \brief decompression filter that doesn't actually decompress
//! A placeholder for cases where a filter is needed that does nothing.
//! Forwards all calls to another stream
class IdInputStream : public google::protobuf::io::ZeroCopyInputStream
{
	private:
		google::protobuf::io::ZeroCopyInputStream *is_ ;

	public:
		IdInputStream( google::protobuf::io::ZeroCopyInputStream *is ) : is_( is ) {}
		~IdInputStream() {}

		bool Next( const void **data, int *size ) { return is_->Next( data, size ) ; }
		void BackUp( int count ) { is_->BackUp( count ) ; }
		bool Skip( int count ) { return is_->Skip( count ) ; }
		int64_t ByteCount() const { return is_->ByteCount() ; }
} ;

class IdOutputStream : public google::protobuf::io::ZeroCopyOutputStream
{
	private:
		google::protobuf::io::ZeroCopyOutputStream *is_ ;

	public:
		IdOutputStream( google::protobuf::io::ZeroCopyOutputStream *is ) : is_( is ) {}
		~IdOutputStream() {}

		bool Next( void **data, int *size ) { return is_->Next( data, size ) ; }
		void BackUp( int count ) { is_->BackUp( count ) ; }
		int64_t ByteCount() const { return is_->ByteCount() ; }
} ;

#if HAVE_LIBZ
#include <zlib.h>

//! \brief decompression filter that uses the zlib
class InflateStream : public google::protobuf::io::ZeroCopyInputStream
{
	private:
		google::protobuf::io::ZeroCopyInputStream *is_ ;

		z_stream zs_ ;
		const Bytef *next_undelivered_ ;
		int64_t total_ ;
		Bytef obuf_[15500] ;

		bool next( const void **data, int *size )
		{
			// nothing in output buffer?  inflate more
			while( zs_.next_out == next_undelivered_ )
			{
				zs_.avail_out = sizeof( obuf_ ) ;
				zs_.next_out = obuf_ ;
				next_undelivered_ = obuf_ ;

				// input buffer empty?  try and get more
				if( !zs_.avail_in ) is_->Next( (const void**)&zs_.next_in, (int*)&zs_.avail_in ) ;

				// inflate some; after that we either have output or
				// need to get more input (or have an error)
				int r = inflate( &zs_, Z_NO_FLUSH ) ;
				// when are we finished?  when we get Z_STREAM_END, but
				// we may still have a block to hand out
				if( r == Z_STREAM_END ) break ;
				else if( r != Z_OK ) throw zs_.msg ;
			}

			*data = (const void*)next_undelivered_ ;
			*size = zs_.next_out - next_undelivered_ ;
			next_undelivered_ = zs_.next_out ;
			total_ += *size ;
			return *size ;
		}

		void back_up( int count ) 
		{
		    next_undelivered_ -= count ;
			total_ -= count ;
		}

		bool skip( int count )
		{
			while( count ) {
				const void *p ; int l ;
				if( !next( &p, &l ) ) return false ;
				count -= l ;
			}
			if( count < 0 ) back_up( -count ) ;
			return true ;
		}

	public:
		InflateStream( google::protobuf::io::ZeroCopyInputStream *is ) : is_( is ), total_( 0 )
		{
			zs_.avail_in = 0 ;

			const void* p ; int l ;
			if( !is_->Next( &p, &l ) ) throw "no data" ;
			try {
				if( l < 2 || ((const uint8_t*)p)[0] != 31 || ((const uint8_t*)p)[1] != 139 )
					throw "not a gzip file" ;

				zs_.zalloc = 0 ;
				zs_.zfree = 0 ;
				zs_.opaque = 0 ;
				if( inflateInit2( &zs_, 15+16 ) != Z_OK ) throw zs_.msg ;

				zs_.next_in = (Bytef*)p ;
				zs_.avail_in = l ;
				zs_.next_out = obuf_ ;
				zs_.avail_out = sizeof( obuf_ ) ;
				next_undelivered_ = obuf_ ;
			}
			catch( ... ) { is_->BackUp( l  ) ; throw ; }
		}

		virtual ~InflateStream() { inflateEnd( &zs_ ) ; if( zs_.avail_in ) is_->BackUp( zs_.avail_in ) ; }
		virtual bool Next( const void **data, int *size ) { return next( data, size ) ; }
		virtual void BackUp( int count ) { back_up( count ) ; }
		virtual bool Skip( int count ) { return skip( count ) ; }
		virtual int64_t ByteCount() const { return total_ ; } 
} ;

class DeflateStream : public google::protobuf::io::ZeroCopyOutputStream
{
	private:
		google::protobuf::io::ZeroCopyOutputStream *os_ ;

		z_stream zs_ ;
		Bytef ibuf_[15500] ;
		int64_t total_ ;

		bool next( void **data, int *size )
		{
			// anything in input buffer?  deflate more
			while( zs_.avail_in )
			{
				// output buffer full?  try and flush it
				if( !zs_.avail_out && !os_->Next( (void**)&zs_.next_out, (int*)&zs_.avail_out ) )
					return false ;

				// deflate some; after that we either have room for input or
				// need to flush more output (or we have an error)
				int r = deflate( &zs_, Z_NO_FLUSH ) ;
				if( r != Z_OK ) throw zs_.msg ;
			}

			zs_.next_in = ibuf_ ;
			zs_.avail_in = sizeof( ibuf_ ) ;

			*data = (void*)ibuf_ ;
			*size = zs_.avail_in ;
			total_ += zs_.avail_in ;
			return true ;
		}

		void back_up( int count ) { zs_.avail_in -= count ; }

	public:
		DeflateStream( google::protobuf::io::ZeroCopyOutputStream *os, int level = Z_DEFAULT_COMPRESSION )
			: os_(os), total_(0) 
		{
			zs_.zalloc = 0 ;
			zs_.zfree = 0 ;
			zs_.opaque = 0 ;
			if( deflateInit2( &zs_, level, Z_DEFLATED, 15+16, 8, Z_DEFAULT_STRATEGY ) != Z_OK ) throw zs_.msg ;

			zs_.next_in = 0 ;
			zs_.avail_in = 0 ;
			zs_.avail_out = 0 ;
		}

		virtual ~DeflateStream()
		{
			// compress whatever is left in buffer or internal state
			for(;;)
			{
				if( !zs_.avail_out && !os_->Next( (void**)&zs_.next_out, (int*)&zs_.avail_out ) ) break ;
				int r = deflate( &zs_, Z_FINISH ) ;
				if( r == Z_STREAM_END ) break ;
				else if( r != Z_OK ) throw zs_.msg ;
			} 
			// give back leftover output buffer, free resources
			if( zs_.avail_out ) os_->BackUp( zs_.avail_out ) ;
			deflateEnd( &zs_ ) ;
		}

		virtual bool Next( void **data, int *size ) { return next( data, size ) ; }
		virtual void BackUp( int count ) { back_up( count ) ; }
		virtual int64_t ByteCount() const { return total_ ; }
} ;

#endif

struct BzipError : public Exception {
	int e_ ;
	BzipError( int e ) : e_(e) { print_to( std::cerr ) ; }
	void print_to( std::ostream& s ) const { s << "bzip error " << e_ ; }
} ;

#if HAVE_LIBBZ2
#include <bzlib.h>

//! \brief decompression filter that uses the libbz2
class BunzipStream : public google::protobuf::io::ZeroCopyInputStream
{
	private:
		google::protobuf::io::ZeroCopyInputStream *is_ ;

		bz_stream zs_ ;
		char *next_undelivered_ ;
		char obuf_[15500] ;

		bool next( const void **data, int *size )
		{
			// nothing in output buffer?  inflate more
			while( zs_.next_out == next_undelivered_ )
			{
				zs_.avail_out = sizeof( obuf_ ) ;
				zs_.next_out = obuf_ ;
				next_undelivered_ = obuf_ ;

				// input buffer empty?  try and get more
				if( !zs_.avail_in ) is_->Next( (const void**)&zs_.next_in, (int*)&zs_.avail_in ) ;

				// inflate some; after that we either have output or
				// need to get more input (or have an error)
				int r = BZ2_bzDecompress( &zs_ ) ;
				// when are we finished?  when we get BZ_STREAM_END, but
				// we may still have a block to hand out
				if( r == BZ_STREAM_END ) break ;
				else if( r < BZ_OK ) throw BzipError( r ) ;
			}

			*data = (const void*)next_undelivered_ ;
			*size = zs_.next_out - next_undelivered_ ;
			next_undelivered_ = zs_.next_out ;
			return *size ;
		}

		void back_up( int count ) { next_undelivered_ -= count ; }

		bool skip( int count )
		{
			while( count ) {
				const void *p ; int l ;
				if( !next( &p, &l ) ) return false ;
				count -= l ;
			}
			if( count < 0 ) back_up( -count ) ;
			return true ;
		}

	public:
		BunzipStream( google::protobuf::io::ZeroCopyInputStream *is ) : is_( is )
		{
			zs_.avail_in = 0 ;

			const void* p ; int l ;
			if( !is_->Next( &p, &l ) ) throw "no data" ;
			try {
				const char* q = (const char*)p ;
				if( l < 3 || q[0] != 'B' || q[1] != 'Z' || q[2] != 'h' )
					throw "not a gzip file" ;

				zs_.bzalloc = 0 ;
				zs_.bzfree = 0 ;
				zs_.opaque = 0 ;
				int r = BZ2_bzDecompressInit( &zs_, 0, 0 ) ;
				if( r < BZ_OK ) throw BzipError( r ) ;

				zs_.next_in = (char*)p ;
				zs_.avail_in = l ;
				zs_.next_out = obuf_ ;
				zs_.avail_out = sizeof( obuf_ ) ;
				next_undelivered_ = obuf_ ;
			}
			catch( ... ) { is_->BackUp( l  ) ; throw ; }
		}

		virtual ~BunzipStream() { BZ2_bzDecompressEnd( &zs_ ) ; if( zs_.avail_in ) is_->BackUp( zs_.avail_in ) ; }
		virtual bool Next( const void **data, int *size ) { return next( data, size ) ; }
		virtual void BackUp( int count ) { back_up( count ) ; }
		virtual bool Skip( int count ) { return skip( count ) ; }
		virtual int64_t ByteCount() const { return (int64_t)zs_.total_out_hi32 << 32 | zs_.total_out_lo32 ; }
} ;

class BzipStream : public google::protobuf::io::ZeroCopyOutputStream
{
	private:
		google::protobuf::io::ZeroCopyOutputStream *os_ ;

		bz_stream zs_ ;
		char ibuf_[15500] ;
		// int64_t total_ ;

		bool next( void **data, int *size )
		{
			// anything in input buffer?  deflate more
			while( zs_.avail_in )
			{
				// output buffer full?  try and flush it
				if( !zs_.avail_out && !os_->Next( (void**)&zs_.next_out, (int*)&zs_.avail_out ) )
					return false ;

				// deflate some; after that we either have room for input or
				// need to flush more output (or we have an error)
				int r = BZ2_bzCompress( &zs_, BZ_RUN ) ;
				if( r < BZ_OK ) throw BzipError( r ) ;
			}

			zs_.next_in = ibuf_ ;
			zs_.avail_in = sizeof( ibuf_ ) ;

			*data = (void*)ibuf_ ;
			*size = zs_.avail_in ;
			// total_ += zs_.avail_in ;
			return true ;
		}

		void back_up( int count ) { zs_.avail_in -= count ; }

	public:
		BzipStream( google::protobuf::io::ZeroCopyOutputStream *os )
			: os_(os) //, total_(0) 
		{
			zs_.bzalloc = 0 ;
			zs_.bzfree = 0 ;
			zs_.opaque = 0 ;
			int r = BZ2_bzCompressInit( &zs_, 9, 0, 0 ) ;
			if( r < BZ_OK ) throw BzipError( r ) ;

			zs_.next_in = 0 ;
			zs_.avail_in = 0 ;
			zs_.avail_out = 0 ;
		}

		virtual ~BzipStream()
		{
			// compress whatever is left in buffer or internal state
			for(;;)
			{
				if( !zs_.avail_out && !os_->Next( (void**)&zs_.next_out, (int*)&zs_.avail_out ) ) break ;
				int r = BZ2_bzCompress( &zs_, BZ_FINISH ) ;
				if( r == BZ_STREAM_END ) break ;
				else if( r < BZ_OK ) throw BzipError( r ) ;
			} 
			// give back leftover output buffer, free resources
			if( zs_.avail_out ) os_->BackUp( zs_.avail_out ) ;
			BZ2_bzCompressEnd( &zs_ ) ;
		}

		virtual bool Next( void **data, int *size ) { return next( data, size ) ; }
		virtual void BackUp( int count ) { back_up( count ) ; }
		virtual int64_t ByteCount() const { return ((int64_t)zs_.total_in_hi32 << 32) | zs_.total_in_lo32 ; }
} ;

#endif

inline google::protobuf::io::ZeroCopyInputStream *decompress( google::protobuf::io::ZeroCopyInputStream *s )
{
#if HAVE_LIBBZ2
	try { return new BunzipStream( s ) ; } catch( ... ) {}
#endif
#if HAVE_LIBZ
	try { return new InflateStream( s ) ; } catch( ... ) {}
#endif
	return new IdInputStream( s ) ;
}

inline google::protobuf::io::ZeroCopyOutputStream *compress_small( google::protobuf::io::ZeroCopyOutputStream *s )
{
#if HAVE_LIBBZ2
	try { return new BzipStream( s ) ; } catch( ... ) {}
#endif
#if HAVE_LIBZ
	try { return new DeflateStream( s, Z_BEST_COMPRESSION ) ; } catch( ... ) {}
#endif
	return new IdOutputStream( s ) ;
}

inline google::protobuf::io::ZeroCopyOutputStream *compress_fast( google::protobuf::io::ZeroCopyOutputStream *s )
{
#if HAVE_LIBZ
	try { return new DeflateStream( s, Z_BEST_SPEED ) ; } catch( ... ) {}
#endif
	return new IdOutputStream( s ) ;
}

#endif
