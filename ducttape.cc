#include "ducttape.h"

namespace streams {

void DuctTaper::put_header( const Header& h ) 
{
	state_ = need_input ;
	if( !h.has_is_sorted_by_coordinate() || 
			h.is_sorted_by_coordinate() != g_ )
		throw "need sorted input for duct taping" ;

	// XXX: probably need to store part of the config?
}

// Strategy:
// If a gap is more likely than a base, we call a deletion and do not
// output anything else.  Else we call the most likely base and a
// quality score.  (We can also call an ambiguity code, choosing it so
// the highest amount of information is retained.)  Then synthesize a
// result to return and change state accordingly.
void DuctTaper::flush_contig()
{
	// XXX  res_ =
	cigar_.clear() ;
	for( int i = 0 ; i != 5 ; ++i ) observed_[i].clear() ;
	state_ = have_output ;
}

void DuctTaper::put_footer( const Footer& f ) 
{
	if( !cigar_.empty() ) flush_contig() ;
	else state_ = end_of_stream ;
}

// Strategy:
// Walk along both cigar strings, observed_ and the sequence of r.  If
// both cigars match, update observed_.  If both cigars insert, update
// observed.  If only r inserts, insert in observed_ and cigar_.  If r
// deletes, update observed_.  That means, we can properly call
// deletions by majority vote, and all inserts are forced to overlap and
// will lead to a majority vote, too (probably wreaking havoc).
void DuctTaper::put_result( const Result& r )
{
	if( has_hit_to( r, g_ ) ) {
		if( !cigar_.empty() && (cur_genome_ != hit_to( r, g_ ).genome_name()
					|| cur_sequence_ != hit_to( r, g_ ).sequence() ) )
			flush_contig() ;

		// XXX need to scan
	}
} 

} // namespace

