#include "ducttape.h"

#include <iostream>

namespace streams {

void DuctTaper::put_header( const Header& h ) 
{
	state_ = need_input ;
	if( !h.has_is_sorted_by_coordinate() || 
			h.is_sorted_by_coordinate() != (g_?g_:"") )
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
	if( is_ins_.empty() ) return ;
	// XXX  res_ =

	// if the last column was added, but nothing was observed, leave it
	// out
	if( observed_.back() == Acc() ) 
	{
		observed_.pop_back() ;
		is_ins_.pop_back() ;
		--contig_end_ ;
	}

	std::clog << "(flush) contig " << cur_genome_ << ':'
		<< cur_sequence_ << ':' << contig_start_ << '-' << contig_end_ 
		<< " (" << (contig_end_ - contig_start_) << ") @ " << nreads_ << std::endl ;

	is_ins_.clear() ;
	observed_.clear() ;
	nreads_ = 0 ;
	// state_ = have_output ;
}

Result DuctTaper::fetch_result()
{
	std::clog << __PRETTY_FUNCTION__ << foot_.IsInitialized() << std::endl ;
	state_ = foot_.IsInitialized() ? end_of_stream : need_input ;
	return res_ ;
}

void DuctTaper::put_footer( const Footer& f ) 
{
	foot_ = f ;
	state_ = end_of_stream ;
	flush_contig() ;
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
	if( !has_hit_to( r, g_ ) ) return ;
	const Hit& h = hit_to( r, g_ ) ;

	// XXX is this good enough?
	if( h.has_diff_to_next() && h.diff_to_next() < 60 ) return ;

	++nreads_ ;
	if( !is_ins_.empty() && (cur_genome_ != h.genome_name()
				|| cur_sequence_ != h.sequence()
				|| h.start_pos() >= contig_end_ ) )
	{
		flush_contig() ;
		cur_genome_ = h.genome_name() ;
		cur_sequence_ = h.sequence() ;
		contig_start_ = contig_end_ = h.start_pos() ;
	}

	// XXX Darn, need to treat RC'd alignments differently.
	bool rev = h.aln_length() < 0 ;
	std::vector< int > cigar( h.cigar().begin(), h.cigar().end() ) ;
	std::string qual = r.read().quality() ;
	std::string seq = r.read().sequence() ;
	if( rev ) {
		std::reverse( cigar.begin(), cigar.end() ) ;
		std::reverse( qual.begin(), qual.end() ) ;
		std::reverse( seq.begin(), seq.end() ) ;
		for( size_t i = 0 ; i != seq.size() ; ++i )
			seq[i] =
				seq[i] == 'A' ? 2 : seq[i] == 'C' ? 3 : 
				seq[i] == 'T' ? 0 : seq[i] == 'G' ? 1 : -1 ;
	}
	else for( size_t i = 0 ; i != seq.size() ; ++i )
	{
		seq[i] = 
			seq[i] == 'A' ? 0 : seq[i] == 'C' ? 1 : 
			seq[i] == 'T' ? 2 : seq[i] == 'G' ? 3 : -1 ;
	}

	std::string::const_iterator p_seq = seq.begin() ;
	size_t column = 0 ;
	size_t cigar_maj = 0 ;
	size_t cigar_min = 0 ;
	while( p_seq != seq.end() && cigar_maj != cigar.size() )
	{
		if( cigar_len( cigar[ cigar_maj ] ) == cigar_min )
		{
			++cigar_maj ;
			cigar_min = 0 ;
		}

		if( column == is_ins_.size() ) {
			is_ins_.push_back( false ) ;
			observed_.push_back( Acc() ) ;
			++contig_end_ ;
		}

		switch( cigar_op( cigar[ cigar_maj ] ) )
		{
			case Hit::Match:
			case Hit::Mismatch:
				if( !is_ins_[column] ) {
					if( *p_seq != -1 ) ++observed_[column].seen[(int)*p_seq] ;

					++cigar_min ;
					++p_seq ;
				}
				break ;

			case Hit::Insert:
				if( !is_ins_[column] ) {
					is_ins_.insert( is_ins_.begin()+column, true ) ;
					observed_.insert( observed_.begin()+column, Acc() ) ;
				}
				if( *p_seq != -1 ) ++observed_[column].seen[ (int)*p_seq ] ;

				++cigar_min ;
				++p_seq ;
				break ;

			case Hit::Delete:
				++observed_[column].seen[4] ; // == gap

				++cigar_min ;
				break ;

			default:
				// all unused
				throw "unexpected CIGAR operation" ;
				break ;
		}
		++column ;
	}

	// std::clog << "contig " << cur_genome_ << ':'
		// << cur_sequence_ << ':' << contig_start_ << '-' << contig_end_ << std::endl ;
} 

} // namespace

