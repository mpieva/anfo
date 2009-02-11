//! \page anfo_combine Combine output files for different indices.
//! This executable combines results for the same sequences, but from
//! different indices into one.  Results should be roughly in the same
//! order in all files to keep memory consumption low, every file should
//! contain a result for every sequence.
//!
//! \note This will also read output files with untagged fields
//! (deprecated, no longer generated) and add the tags.  Warnings are
//! printed on old files, deserialization errors, and footers with error
//! codes.  Exit code is the logical OR of all error codes found in
//! footers, if a deserialization error occured, this is again OR'ed
//! with 1.

#include "outputfile.h"

#include <fstream>
#include <iostream>
#include <map>

using namespace google::protobuf::io ;
using namespace output ;
using namespace std ;

void merge_results( Result& lhs, const Result& rhs )
{
	// How to merge what...
	// - seqid, description, sequence, trimpoints: all equal, no merging needed

	// - reason: do something sensible...
	if( lhs.reason() != aligned && rhs.reason() != no_policy && rhs.reason() != no_seeds )
	{
		if( rhs.reason() == aligned ) lhs.clear_reason() ;
		else if( lhs.reason() == no_seeds || lhs.reason() == no_policy ) lhs.set_reason( rhs.reason() ) ;
		else if( lhs.reason() == too_many_seeds && rhs.reason() == bad_alignment ) lhs.set_reason( rhs.reason() ) ;
	}

	// - num_xxx: add them
	lhs.set_num_raw_seeds( lhs.num_raw_seeds() + rhs.num_raw_seeds() ) ;
	lhs.set_num_grown_seeds( lhs.num_grown_seeds() + rhs.num_grown_seeds() ) ;
	lhs.set_num_clumps( lhs.num_clumps() + rhs.num_clumps() ) ;

	// - best_hit, diff_to_next_species, diff_to_next_order: TODO
	// - best_to_genome: take better hit, recalculate diff_to_next{,_chromosme{,_class}}
	if( rhs.has_best_to_genome() )
	{
		if( lhs.has_best_to_genome() )
		{
			// two hits, this is work...
			if( lhs.best_to_genome().score() <= rhs.best_to_genome().score() )
			{
				// left is better
				if( !lhs.has_diff_to_next() || lhs.best_to_genome().score() + lhs.diff_to_next() > rhs.best_to_genome().score() )
					lhs.set_diff_to_next( rhs.best_to_genome().score() - lhs.best_to_genome().score() ) ;
				//! \todo to chromosome, chromosome class? dunno...
			}
			else
			{
				// right is better
				if( !rhs.has_diff_to_next() || rhs.best_to_genome().score() + rhs.diff_to_next() > lhs.best_to_genome().score() )
					lhs.set_diff_to_next( lhs.best_to_genome().score() - rhs.best_to_genome().score() ) ;
				*lhs.mutable_best_to_genome() = rhs.best_to_genome() ;
				//! \todo to chromosome, chromosome class? dunno...
			}
		}
		else
		{
			// no hit at left side --> just assign
			*lhs.mutable_best_to_genome() = rhs.best_to_genome() ;
			if( rhs.has_diff_to_next() ) lhs.set_diff_to_next( rhs.diff_to_next() ) ;
			if( rhs.has_diff_to_next_chromosome() ) lhs.set_diff_to_next_chromosome( rhs.diff_to_next_chromosome() ) ;
			if( rhs.has_diff_to_next_chromosome_class() ) lhs.set_diff_to_next_chromosome_class( rhs.diff_to_next_chromosome_class() ) ;
		}
	}
}

int main_( int argc, const char** argv )
{
    vector<AnfoFile*> inputs ;
    for( const char** arg = argv+1 ; arg != argv+argc ; ++arg )
    {
	clog << "open " << *arg << endl ;
	inputs.push_back( new AnfoFile( *arg ) ) ;
    }
    clog << inputs.size() << " input files" << endl ;

    OstreamOutputStream oos( &cout ) ;
    CodedOutputStream cos( &oos ) ;
    cos.WriteRaw( "ANFO", 4 ) ;
    Header hdr ;
    for( size_t i = 0 ; i != inputs.size() ; ++i )
	hdr.MergeFrom( inputs[i]->read_header() ) ;
    write_delimited_message( cos, 1, hdr ) ;

    typedef map< string, pair< size_t, Result > > Buffer ;
	Buffer buffer ;
    bool cont ;
    int nread = 0, nwritten = 0 ;
	do {
		cont = false ;
		for( size_t i = 0 ; i != inputs.size() ; ++i )
		{
			Result r = inputs[i]->read_result() ;
			if( r.has_seqid() ) {
				cont = true ;
				nread++ ;
				pair< size_t, Result > &p = buffer[ r.seqid() ] ;
				p.first++ ;
				if( p.second.has_seqid() ) merge_results( p.second, r ) ;
				else p.second = r ;

				if( p.first == inputs.size() ) 
				{
					nwritten++ ;
					write_delimited_message( cos, 2, p.second ) ;
					buffer.erase( r.seqid() ) ;
				}
			}
		}
		if( ((nread + nwritten) & 0xFFFF) == 0 ) {
			clog << "\033[KRead " << nread << ", wrote " << nwritten << ", buffering "
				<< buffer.size() << '\r' << flush ;
		}
	} while( cont ) ;
		
	for( Buffer::const_iterator b = buffer.begin(), e = buffer.end() ; b!=e ; ++b )
		write_delimited_message( cos, 2, b->second.second ) ;

    Footer f ;
    f.set_exit_code( 0 ) ;
    for( size_t i = 0 ; i != inputs.size() ; ++i )
	{
		Footer f_ = inputs[i]->read_footer() ;
		if( f_.has_exit_code() ) f.set_exit_code( f.exit_code() | f_.exit_code() ) ;
		delete inputs[i] ;
	}
    write_delimited_message( cos, 3, f ) ;
    return f.exit_code() ;
}


