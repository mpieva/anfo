#include "outputfile.h"
#include "compress_stream.h"
#include "util.h"

#include <iostream>
#include <fcntl.h>

using namespace google::protobuf::io ;
using namespace output ;
using namespace std ;

AnfoFile::AnfoFile( const std::string& name )
	: iis_( throw_errno_if_minus1( open( name.c_str(), O_RDONLY ), "opening ", name.c_str() ) )
	, zis_( decompress( &iis_ ) ), error_( false ), name_( name )
{
	iis_.SetCloseOnDelete( true ) ;
	check_valid_file() ;
}

AnfoFile::AnfoFile( const std::string& name, int fd )
	: iis_( fd ), zis_( decompress( &iis_ ) ), error_( false ), name_( name )
{
	iis_.SetCloseOnDelete( true ) ;
	check_valid_file() ;
}

void AnfoFile::check_valid_file()
{
	std::string tag ;
	CodedInputStream cis( zis_.get() ) ;

	if( !cis.ReadString( &tag, 4 ) || tag != "ANFO" ) {
		clog << "\033[K" << name_ << ": not an ANFO file" << endl ;
		error_ = true ;
	}
	foot_.set_exit_code(0) ;
}

Header AnfoFile::read_header()
{
	uint32_t tag ;
	Header hdr ;
	CodedInputStream cis( zis_.get() ) ;
	if( cis.ReadVarint32( &tag ) && tag == 10 )
	{
		if( cis.ReadVarint32( &tag ) )
		{
			int lim = cis.PushLimit( tag ) ;
			if( hdr.ParseFromCodedStream( &cis ) ) {
				cis.PopLimit( lim ) ;
				return hdr ;
			}
		}
	}

	clog << "\033[K" << name_ << ": deserialization error" << endl ;
	error_ = true ;
	return hdr ;
}

Result AnfoFile::read_result()
{
	uint32_t tag = 0 ;
	CodedInputStream cis( zis_.get() ) ;
	if( !error_ ) {
		if( cis.ExpectAtEnd() ) {
			clog << "\033[K" << name_ << ": end of stream" << endl ;
			return Result() ;
		}
		if( tag = cis.ReadTag() ) {
			uint32_t size ;
			std::string buf ;
			if( cis.ReadVarint32( &size ) && cis.ReadString( &buf, size ) ) {
				if( tag == 18 ) {
					Result res ;
					if( res.ParseFromString( buf ) && res.has_seqid() ) {
						return res ;
					}
				}
				if( tag == 26 ) {
					if( foot_.ParseFromString( buf ) ) return Result() ;
				}
				clog << "\033[K" << name_ << ": deserialization error" << endl ;
				error_ = true ;
			}
		}
	}
	return Result() ;
}


void merge_sensibly( Header& lhs, const Header& rhs )
{
	//! \todo need to reduce junk headers
	lhs.MergeFrom( rhs ) ;
}

//! \brief merges two results, keeping the best hit
void merge_sensibly( Result& lhs, const Result& rhs )
{
	// How to merge what...
	// - seqid, description, sequence, trimpoints: all equal, no merging needed

	// - reason: do the sensible thing...
	if( lhs.reason() != aligned && rhs.reason() != no_policy && rhs.reason() != no_seeds )
	{
		if( rhs.reason() == aligned ) lhs.clear_reason() ;
		else if( lhs.reason() == no_seeds || lhs.reason() == no_policy ) lhs.set_reason( rhs.reason() ) ;
		else if( lhs.reason() == too_many_seeds && rhs.reason() == bad_alignment ) lhs.set_reason( rhs.reason() ) ;
	}

	// - num_xxx: just add them
	lhs.set_num_raw_seeds( lhs.num_raw_seeds() + rhs.num_raw_seeds() ) ;
	lhs.set_num_grown_seeds( lhs.num_grown_seeds() + rhs.num_grown_seeds() ) ;
	lhs.set_num_clumps( lhs.num_clumps() + rhs.num_clumps() ) ;

	// - best_hit, diff_to_next_species, diff_to_next_order: TODO
	
	// - best_to_genome: take better hit, recalculate diff_to_next{,_chromosome{,_class}}
	if( rhs.has_best_to_genome() )
	{
		if( lhs.has_best_to_genome() )
		{
			// two hits, this is work...
			if( lhs.best_to_genome().score() <= rhs.best_to_genome().score() )
			{
				// left is better
				if( !lhs.has_diff_to_next() ||
						lhs.best_to_genome().score() + lhs.diff_to_next() > rhs.best_to_genome().score() )
					lhs.set_diff_to_next( rhs.best_to_genome().score() - lhs.best_to_genome().score() ) ;

				//! \todo diff to chromosome, chromosome class? dunno...
			}
			else
			{
				// right is better
				if( !rhs.has_diff_to_next() ||
						rhs.best_to_genome().score() + rhs.diff_to_next() > lhs.best_to_genome().score() )
					lhs.set_diff_to_next( lhs.best_to_genome().score() - rhs.best_to_genome().score() ) ;

				*lhs.mutable_best_to_genome() = rhs.best_to_genome() ;
				//! \todo diff to chromosome, chromosome class? dunno...
			}
		}
		else
		{
			// no hit at left side --> just assign
			*lhs.mutable_best_to_genome() = rhs.best_to_genome() ;
			if( rhs.has_diff_to_next() )
				lhs.set_diff_to_next( rhs.diff_to_next() ) ;

			if( rhs.has_diff_to_next_chromosome() )
				lhs.set_diff_to_next_chromosome( rhs.diff_to_next_chromosome() ) ;

			if( rhs.has_diff_to_next_chromosome_class() )
				lhs.set_diff_to_next_chromosome_class( rhs.diff_to_next_chromosome_class() ) ;
		}
	}
}

//! \brief merges two footers
//! Anything unexpected is simply merged, the resulting exit code is the
//! logical or of the two inputs.
void merge_sensibly( output::Footer& lhs, const output::Footer& rhs )
{
	int exit_code = lhs.exit_code() | rhs.exit_code() ;
	lhs.MergeFrom( rhs ) ;
	lhs.set_exit_code( exit_code ) ;
}


