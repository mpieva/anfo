#include "ducttape.h"

#include <algorithm>
#include <numeric>
#include <sstream>

namespace streams {

void DuctTaper::put_header( const Header& h ) 
{
	if( !h.has_is_sorted_by_coordinate() || 
			h.is_sorted_by_coordinate() != (g_?g_:"") )
		throw "need sorted input for duct taping" ;

	hdr_ = h ;
	state_ = need_input ;
}

// Strategy:
// If a gap is more likely than a base, we call a deletion and do not
// output anything else.  Else we call the most likely base and a
// quality score.  (We can also call an ambiguity code, choosing it so
// the highest amount of information is retained.)  Then synthesize a
// result to return and change state accordingly.
//
// XXX: need a rule to decide whether sth. was inserted and how much
//
// We need to synthesize a 'Result' here including a 'Read' and a 'Hit'.

void DuctTaper::flush_contig()
{
	if( is_ins_.empty() ) return ;


	// if the last column was added, but nothing was observed, leave it
	// out
	if( observed_.back().pristine() )
	{
		observed_.pop_back() ;
		is_ins_.pop_back() ;
		--contig_end_ ;
	}

	res_.Clear() ;

	Read &rd = *res_.mutable_read() ;
	std::stringstream ss1, ss2 ;
	ss1 << cur_genome_ << ':' << cur_sequence_ << ':' << contig_start_ ;
	rd.set_seqid( ss1.str() ) ;
	
	ss2 << nreads_ << " reads joined" ;
	rd.set_description( ss2.str() ) ;

	// XXX remove columns where sth. was inserted, but most reads show a gap
	
	int eff_length = 0 ;
	for( size_t i = 0 ; i != is_ins_.size() ; ++i )
		if( std::accumulate( observed_[i].seen, observed_[i].seen+4, 0 ) > observed_[i].seen[4] )
			++eff_length ;

    for( int i = 0 ; i != 4 /*10*/ ; ++i ) rd.add_likelihoods() ;

	std::vector<unsigned> cigar ;
	for( size_t i = 0, ie = 0 ; i != is_ins_.size() ; ++i )
	{
		int cov = std::accumulate( observed_[i].seen, observed_[i].seen+4, 0 ) ;
		if( cov > observed_[i].seen[4] ) {
			if( is_ins_[i] ) push_i( cigar, 1 ) ; else push_m( cigar, 1 ) ;

			Logdom lk_tot = observed_[i].lk[0] + observed_[i].lk[1] + observed_[i].lk[2] + observed_[i].lk[3] ;
			int maxlk = 0 ;
			for( int j = 0 ; j != 4 ; ++j )
			{
				rd.add_seen_bases( observed_[i].seen[j] ) ;
				if( observed_[i].lk[j] > observed_[i].lk[maxlk] ) maxlk = j ;
			}

			rd.add_depth( cov + observed_[i].seen[4] ) ;
			rd.mutable_sequence()->push_back( "ACTG"[maxlk] ) ;

			Logdom q = maxlk == 0 ? observed_[i].lk[1] : observed_[i].lk[0] ;
			for( int j = maxlk == 0 ? 2 : 1 ; j != 4 ; ++j )
				if( j != maxlk ) q += observed_[i].lk[j] ;

			rd.mutable_quality()->push_back( (q / lk_tot).to_phred_byte() ) ;

			for( int j = 0 ; j != 4 ; ++j )
				rd.mutable_likelihoods(j)->push_back( 
						(observed_[i].lk[j] / observed_[i].lk[maxlk]).to_phred_byte() ) ;

			++ie ;
		}
		else push_d( cigar, 1 ) ;
	}


	Hit &hit = *res_.add_hit() ;
	hit.set_genome_name( cur_genome_ ) ;
	hit.set_sequence( cur_sequence_ ) ;
	hit.set_start_pos( contig_start_ ) ;
	hit.set_aln_length( contig_end_ - contig_start_ ) ;
	hit.set_score( 0 ) ; // XXX cannot calculate the friggen' score 
	std::copy( cigar.begin(), cigar.end(), RepeatedFieldBackInserter( hit.mutable_cigar() ) ) ;
	// set diff_to_next as RMS of all the diff_to_nexts that went in

	is_ins_.clear() ;
	observed_.clear() ;
	nreads_ = 0 ;
	state_ = have_output ;
}

Result DuctTaper::fetch_result()
{
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
//
// XXX likelihoods do not regard BJ model
void DuctTaper::put_result( const Result& r )
{
	if( !has_hit_to( r, g_ ) ) return ;
	const Hit& h = hit_to( r, g_ ) ;

	++nreads_ ;
	if( cur_genome_ != h.genome_name()
				|| cur_sequence_ != h.sequence()
				|| h.start_pos() > contig_end_ )
	{
		flush_contig() ;
		cur_genome_ = h.genome_name() ;
		cur_sequence_ = h.sequence() ;
		contig_start_ = contig_end_ = h.start_pos() ;
	}

	// Damn, need to treat RC'd alignments differently.
	bool rev = h.aln_length() < 0 ;
	std::vector< int > cigar( h.cigar().begin(), h.cigar().end() ) ;
	std::string qual = r.read().has_quality() ? r.read().quality() : std::string( r.read().sequence().size(), 30 ) ;
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
	std::string::const_iterator q_seq = qual.begin() ;
	size_t column = 0 ;
	size_t cigar_maj = 0 ;
	size_t cigar_min = 0 ;

	for( int offs = h.start_pos() - contig_start_ ; offs ; ++column )
		if( !is_ins_[column] ) --offs ;

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

		// XXX: likelihood calculation doesn't take BJ model into
		// account
		switch( cigar_op( cigar[ cigar_maj ] ) )
		{
			case Hit::Match:
			case Hit::Mismatch:
				if( !is_ins_[column] ) {
					if( *p_seq != -1 ) {
						++observed_[column].seen[(int)*p_seq] ;
						for( int k = 0 ; k != 4 ; ++k )
						{
							observed_[column].lk[k] *= 
								k == *p_seq ? (Logdom::from_float(1) - Logdom::from_phred( *q_seq ))
								            : Logdom::from_phred( *q_seq ) / Logdom::from_float(3) ;
						}
					}

					++cigar_min ;
					++p_seq ;
					++q_seq ;
				}
				break ;

			case Hit::Insert:
				if( !is_ins_[column] ) {
					is_ins_.insert( is_ins_.begin()+column, true ) ;
					observed_.insert( observed_.begin()+column, Acc() ) ;
				}
				if( *p_seq != -1 ) {
					++observed_[column].seen[ (int)*p_seq ] ;
					for( int k = 0 ; k != 4 ; ++k )
					{
						observed_[column].lk[k] *= 
							k == *p_seq ? (Logdom::from_float(1) - Logdom::from_phred( *q_seq ))
							: Logdom::from_phred( *q_seq ) / Logdom::from_float(3) ;
					}
				}

				++cigar_min ;
				++p_seq ;
				++q_seq ;
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
} 

} // namespace

