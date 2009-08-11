#include "ducttape.h"

#include <algorithm>
#include <cmath>
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

void DuctTaper::flush_contig()
{
	if( observed_.empty() ) return ;


	// if the last column was added, but nothing was observed, leave it
	// out
	if( observed_.back().pristine() )
	{
		observed_.pop_back() ;
		--contig_end_ ;
	}

	res_.Clear() ;
	res_.set_num_reads( nreads_ ) ;

	Read &rd = *res_.mutable_read() ;
	std::stringstream ss1 ;
	ss1 << cur_genome_ << ':' << cur_sequence_ << ':' << contig_start_ ;
	rd.set_seqid( ss1.str() ) ;
	
	// remove cols where sth. was inserted, but most reads show a gap
	for( Accs::iterator i = observed_.begin() ; i != observed_.end() ; )
	{
		if( std::accumulate( i->seen, i->seen+4, 0 ) <= i->seen[4] && i->is_ins )
			i = observed_.erase(i) ;
		else ++i ;
	}

    for( int i = 0 ; i != 10 ; ++i ) rd.add_likelihoods() ;

	std::vector<unsigned> cigar ;
	for( Accs::iterator i = observed_.begin() ; i != observed_.end() ; ++i )
	{
		int cov = std::accumulate( i->seen, i->seen+4, 0 ) ;
		if( cov > i->seen[4] ) {
			if( i->is_ins ) push_i( cigar, 1 ) ; else push_m( cigar, 1 ) ;

			Logdom lk_tot = std::accumulate( i->lk+1, i->lk+10, i->lk[0] ) ;
			int maxlk = 0 ;
			for( int j = 0 ; j != 4 ; ++j )
			{
				rd.add_seen_bases( i->seen[j] ) ;
				if( i->lk[j] > i->lk[maxlk] ) maxlk = j ;
			}

			rd.add_depth( cov + i->seen[4] ) ;
			rd.mutable_sequence()->push_back( "ACGTMRWSYK"[maxlk] ) ;

			Logdom q = maxlk == 0 ? i->lk[1] : i->lk[0] ;
			for( int j = maxlk == 0 ? 2 : 1 ; j != 4 ; ++j )
				if( j != maxlk ) q += i->lk[j] ;

			rd.mutable_quality()->push_back( (q / lk_tot).to_phred_byte() ) ;

			for( int j = 0 ; j != 10 ; ++j )
				rd.mutable_likelihoods(j)->push_back( 
						(i->lk[j] / i->lk[maxlk]).to_phred_byte() ) ;

		}
		else push_d( cigar, 1 ) ;
	}


	Hit &hit = *res_.add_hit() ;
	hit.set_genome_name( cur_genome_ ) ;
	hit.set_sequence( cur_sequence_ ) ;
	hit.set_start_pos( contig_start_ ) ;
	hit.set_aln_length( contig_end_ - contig_start_ ) ;
	hit.set_diff_to_next( 0.5 + std::sqrt( mapq_accum_ / nreads_ ) ) ;
	std::copy( cigar.begin(), cigar.end(), RepeatedFieldBackInserter( hit.mutable_cigar() ) ) ;

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

	int mapq = h.has_diff_to_next() ? h.diff_to_next() : 254 ;
	mapq_accum_ += mapq*mapq ;

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
				seq[i] == 'G' ? 0 : seq[i] == 'T' ? 1 : -1 ;
	}
	else for( size_t i = 0 ; i != seq.size() ; ++i )
	{
		seq[i] = 
			seq[i] == 'A' ? 0 : seq[i] == 'C' ? 1 : 
			seq[i] == 'G' ? 2 : seq[i] == 'T' ? 3 : -1 ;
	}

	std::string::const_iterator p_seq = seq.begin() ;
	std::string::const_iterator q_seq = qual.begin() ;
	Accs::iterator column = observed_.begin() ;
	size_t cigar_maj = 0 ;
	size_t cigar_min = 0 ;

	for( int offs = h.start_pos() - contig_start_ ; offs ; ++column )
		if( !column->is_ins ) --offs ;

	while( p_seq != seq.end() && cigar_maj != cigar.size() )
	{
		if( cigar_len( cigar[ cigar_maj ] ) == cigar_min )
		{
			++cigar_maj ;
			cigar_min = 0 ;
		}

		if( column == observed_.end() ) {
			column = observed_.insert( column, Acc() ) ;
			++contig_end_ ;
		}

		// XXX: likelihood calculation doesn't take BJ model into
		// account
		switch( cigar_op( cigar[ cigar_maj ] ) )
		{
			case Hit::Match:
			case Hit::Mismatch:
			case Hit::Insert:
				if( cigar_op( cigar[ cigar_maj ] ) == Hit::Insert )
				{
					if( !column->is_ins ) column = observed_.insert( column, Acc(true, column->crossed ) ) ;
				} 
				else
				{
					if( column->is_ins ) {
						++column->gapped ;
						break ;
					}
				}

				if( *p_seq != -1 ) {
					int base = *p_seq ;
					++column->seen[ base ] ;

					static const bool halfmat[4][10] = {
						{ 1, 1, 1, 0, 0, 0 },
						{ 1, 0, 0, 1, 1, 0 },
						{ 0, 1, 0, 1, 0, 1 },
						{ 0, 0, 1, 0, 1, 1 } } ;

					Logdom l_mat = Logdom::from_float(1) - Logdom::from_phred( *q_seq ) ;
					Logdom l_mismat = Logdom::from_phred( *q_seq ) / Logdom::from_float(3) ;
					Logdom l_half = (l_mat + l_mismat) / Logdom::from_float(2) ;

					for( int k = 0 ; k != 10 ; ++k )
						column->lk[k] *= k == base ? l_mat : halfmat[base][k-4] ? l_half : l_mismat ;
				}

				++cigar_min ;
				++p_seq ;
				++q_seq ;
				break ;

			case Hit::Delete:
				++column->gapped ;
				++cigar_min ;
				break ;

			default: 	// all others are unused
				throw "unexpected CIGAR operation" ;
				break ;
		}
		if( p_seq != seq.begin() ) ++column->crossed ;
		++column ;
	}
} 

} // namespace

