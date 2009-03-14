#ifndef INCLUDED_COMPRESS_STREAM_H
#define INCLUDED_COMPRESS_STREAM_H

#include <google/protobuf/io/zero_copy_stream.h>
#include <memory>

#include <zlib.h>
#include <bzlib.h>

//! \brief decompression filter that doesn't actually decompress
//! A placeholder for cases where a filter is needed that does nothing.
//! Forwards all calls to another stream
class IdStream : public google::protobuf::io::ZeroCopyInputStream
{
	private:
		google::protobuf::io::ZeroCopyInputStream *is_ ;

	public:
		IdStream( google::protobuf::io::ZeroCopyInputStream *is ) : is_( is ) {}
		~IdStream() {}

		bool Next( const void **data, int *size ) { return is_->Next( data, size ) ; }
		void BackUp( int count ) { is_->BackUp( count ) ; }
		bool Skip( int count ) { return is_->Skip( count ) ; }
		int64_t ByteCount() const { return is_->ByteCount() ; }
} ;

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

//! \brief decompression filter that uses the libbz2
class BunzipStream : public google::protobuf::io::ZeroCopyInputStream
{
	public:
		//  XXX: needs to be handled
		struct BzipError {
			int e_ ;
			BzipError( int e ) : e_(e) {}
		} ;

	private:
		google::protobuf::io::ZeroCopyInputStream *is_ ;

		bz_stream zs_ ;
		char *next_undelivered_ ;
		// int64_t total_ ;
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
				else if( r != BZ_OK ) throw BzipError( r ) ;
			}

			*data = (const void*)next_undelivered_ ;
			*size = zs_.next_out - next_undelivered_ ;
			next_undelivered_ = zs_.next_out ;
			// total_ += *size ;
			return *size ;
		}

		void back_up( int count ) 
		{
		    next_undelivered_ -= count ;
			// total_ -= count ;
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
		BunzipStream( google::protobuf::io::ZeroCopyInputStream *is ) : is_( is ) // , total_( 0 )
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
				if( r != BZ_OK ) throw BzipError( r ) ;

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

inline google::protobuf::io::ZeroCopyInputStream *decompress( google::protobuf::io::ZeroCopyInputStream *s )
{
	try { return new BunzipStream( s ) ; } catch( ... ) {}
	try { return new InflateStream( s ) ; } catch( ... ) {}
	return new IdStream( s ) ;
}

#endif
