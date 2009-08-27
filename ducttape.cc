#include "ducttape.h"

#include "align.h"
#include "index.h"

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
	simple_adna::configure( h.config().aligner(), 0 ) ;
	state_ = need_input ;
}

//! \brief generates an actual contig from accumulated data
//! \internal
//! The strategy: If a gap is more likely than a base, we call a
//! deletion and do not output anything else.  Else we compute
//! likelihoods for all possible dialleles.  From that we call the most
//! likely base or an ambiguity code if a heterozygote is more likely.
//! Finally store the new contig as an \c output::Result and update
//! internal state.

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

			Logdom lk_tot = std::accumulate( i->lk, i->lk+10, Logdom::null() ) ;
			int maxlk = 0 ;
			for( int j = 0 ; j != 4 ; ++j )
			{
				rd.add_seen_bases( i->seen[j] ) ;
				if( i->lk[j] > i->lk[maxlk] ) maxlk = j ;
			}

			rd.add_depth( cov + i->seen[4] ) ;
			rd.mutable_sequence()->push_back( "ACGTMRWSYK"[maxlk] ) ;

			Logdom q = Logdom::null() ;
			for( int j = 0 ; j != 4 ; ++j ) if( j != maxlk ) q += i->lk[j] ;

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
	hit.set_diff_to_next( int( 0.5 + std::sqrt( mapq_accum_ / nreads_ ) ) ) ;
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


// A bit annoying:  we need to treat RC'd RC'd alignments differently.
class AlnIter
{
	private:
		bool fwd_ ;
		bool has_qual_ ;
		google::protobuf::RepeatedField<unsigned>::const_iterator cigar_ ;
		size_t cigar_minor_ ;
		std::string::const_iterator seq_, qual_ ;

		void norm_cigar()
		{
			if( fwd_ ) while( cigar_len( *cigar_ ) == cigar_minor_ )
			{
				++cigar_ ;
				cigar_minor_ = 0 ;
			}
			else while( cigar_len( cigar_[-1] ) == cigar_minor_ )
			{
				--cigar_ ;
				cigar_minor_ = 0 ;
			}
		}


	public:
		AlnIter( const output::Read& r, const output::Hit& h ) :
			fwd_( h.aln_length() >= 0 ), has_qual_( r.has_quality() ),
			cigar_( fwd_ ? h.cigar().begin() : h.cigar().end() ), cigar_minor_( 0 ),
			seq_( fwd_ ? r.sequence().begin() : r.sequence().end() )
		{
			if( has_qual_ ) qual_ = fwd_ ? r.quality().begin() : r.quality().end() ;
			norm_cigar() ;
		}

		AlnIter( const output::Read& r, const output::Hit& h, int ) :
			fwd_( h.aln_length() >= 0 ), has_qual_( r.has_quality() ),
			cigar_( fwd_ ? h.cigar().end() : h.cigar().begin() ), cigar_minor_( 0 ),
			seq_( fwd_ ? r.sequence().end() : r.sequence().begin() )
		{
			if( has_qual_ ) qual_ = fwd_ ? r.quality().end() : r.quality().begin() ;
		}

		bool operator != ( const AlnIter& rhs )
		{
			return seq_ != rhs.seq_ && qual_ != rhs.qual_ && 
				( cigar_ != rhs.cigar_ || cigar_minor_ != rhs.cigar_minor_ ) ;
		}

		AlnIter& operator ++ ()
		{
			if( cigar_op() != Hit::Delete ) 
			{
				if( fwd_ ) { ++seq_ ; ++qual_ ; }
				else       { --seq_ ; --qual_ ; }
			}
			++cigar_minor_ ;
			norm_cigar() ;
			return *this ;
		}

		Hit::Operation cigar_op() const { return streams::cigar_op( fwd_ ? cigar_[0] : cigar_[-1] ) ; }
	
		int base() const { 
			return fwd_ ? ( seq_[ 0] == 'A' ? 0 : seq_[ 0] == 'C' ? 1 : 
			                seq_[ 0] == 'G' ? 3 : seq_[ 0] == 'T' ? 2 : -1 )
						: ( seq_[-1] == 'A' ? 2 : seq_[-1] == 'C' ? 3 : 
				            seq_[-1] == 'G' ? 1 : seq_[-1] == 'T' ? 0 : -1 ) ;
		}

		Logdom qual() const { return Logdom::from_phred( fwd_ ? qual_[0] : qual_[-1] ) ; }
} ;

//! \brief updates likelihood information
//! The strategy:  Walk along both cigar strings, observed_ and the
//! sequence of r.  If both cigars do a match, update observed_.  If
//! both cigars insert, update observed_.  If only r inserts, insert in
//! observed_ and cigar_.  If r deletes, update observed_.  That means,
//! we can properly call deletions by majority vote, and all inserts are
//! forced to overlap and will lead to a majority vote, too (probably
//! wreaking havoc in the rare case of apparently polymorphic inserts).
//!
//! In observed_ we store likelihoods for observations assuming a given
//! base.  The final base call inverts this by dividing by the total
//! likelihood.  Incorporation of the BJ model is quite easy: we just
//! multiply the likelihood of a C becoming a T onto whatever the
//! quality scores yield where appropriate.  (That also means we need to
//! multiply the probability of a C staying a C onto the appropriate
//! likelihoods.)  The likelihood depends on how likely we are to be in
//! a single stranded region.
//!
//! The likelihood of the
//! alignment with a given overhang length is simply the alignment
//! score, by dividing by the total likelihood, we get a likelihood for
//! that overhang length.
//!

//! Likelihood for having an ss overhang of at least n is just the sum
//! of match scores in first n positions; therefore probability is
//! L_ss(n) / (L_ss(n)+L_ds(n)) * P_ss(n), the latter being the prior of
//! overhang_ext_penalty*n + overhang_enter_penalty.   Incremental
//! calculation is therfore easy in the forward direction.  C->T
//! probabilty is then the weighted sum of the two deamination rates, C
//! staying C is the converse.
//
//! Backward direction:  same calculation, but we need to initialize
//! to the sum of all the match scores, which requires a preliminary
//! pass over everything.

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

	Accs::iterator column = observed_.begin() ;
	for( int offs = h.start_pos() - contig_start_ ; offs ; ++column )
		if( !column->is_ins ) --offs ;

	Logdom lk_ss_5 = simple_adna::overhang_enter_penalty, lk_ds_5 ;
	Logdom lk_ss_3 = simple_adna::overhang_enter_penalty, lk_ds_3 ;

	DnaP ref = Metagenome::find_sequence( h.genome_name(), h.sequence() ).find_pos( h.sequence(), h.start_pos() ) ;
	for( AlnIter aln_i( r.read(), h ), aln_e( r.read(), h, 1 ) ; aln_i != aln_e ; ++aln_i )
	{
		if( aln_i.cigar_op() == Hit::Match || aln_i.cigar_op() == Hit::Mismatch )
		{
			lk_ss_3 *= simple_adna::ss_mat[ 1 << *ref ][ 1 << aln_i.base() ] ;
			lk_ds_3 *= simple_adna::ds_mat[ 1 << *ref ][ 1 << aln_i.base() ] ;
		}
		if( aln_i.cigar_op() != Hit::Insert && aln_i.cigar_op() != Hit::SoftClip ) ++ref ;
		if( aln_i.cigar_op() != Hit::Delete ) lk_ds_3 *= simple_adna::overhang_ext_penalty ;
	}

	ref = Metagenome::find_sequence( h.genome_name(), h.sequence() ).find_pos( h.sequence(), h.start_pos() ) ;
	for( AlnIter aln_b( r.read(), h ), aln_i( aln_b ), aln_e( r.read(), h, 1 ) ;
			aln_i != aln_e ; ++column )
	{
		if( column == observed_.end() ) {
			column = observed_.insert( column, Acc() ) ;
			++contig_end_ ;
		}

		// XXX: likelihood calculation doesn't take BJ model into
		// account
		switch( aln_i.cigar_op() )
		{
			case Hit::SoftClip:
			case Hit::Insert:
				if( !column->is_ins ) column = observed_.insert( column, Acc(true, column->crossed ) ) ;
				goto no_match ;

			case Hit::Match:
			case Hit::Mismatch:
				if( column->is_ins ) {
					++column->gapped ;
					break ;
				}
				++ref ;

				lk_ss_5 *= simple_adna::ss_mat[ complement(1 << *ref) ][ complement(1 << aln_i.base()) ] ;
				lk_ds_5 *= simple_adna::ds_mat[ complement(1 << *ref) ][ complement(1 << aln_i.base()) ] ;
				lk_ss_3 /= simple_adna::ss_mat[ 1 << *ref ][ 1 << aln_i.base() ] ;
				lk_ds_3 /= simple_adna::ds_mat[ 1 << *ref ][ 1 << aln_i.base() ] ;

no_match:
				lk_ds_3 /= simple_adna::overhang_ext_penalty ;
				lk_ds_5 /= simple_adna::overhang_ext_penalty ;

				if( aln_i.base() != -1 ) {
					++column->seen[ aln_i.base() ] ;

					Logdom prob_ss_5 = lk_ss_5 / (lk_ss_5 + lk_ds_5) ;
					Logdom prob_ss_3 = lk_ss_3 / (lk_ss_3 + lk_ds_3) ;
					Logdom l_mat = 1 - aln_i.qual(), l_mismat = aln_i.qual() / 3 ;

					// lk_base[x] is probability of seeing aln_i.base()
					// given that an x was in the real sequence.  It is
					// the products of the probability of x being
					// modified to y times the probability of y being
					// seen as aln_i.base() summed over all possible y.
					// The former comes directly from blended subst
					// matrices, the latter is l_mat or l_mismat.
					// Probabilities for dialleles are simply the
					// average of probabilities for the constituents.
					
					Logdom lk_base[4] ;
					for( int k = 0 ; k != 4 ; ++k )
					{
						lk_base[k] = Logdom::null() ;
						for( int l = 0 ; l != 4 ; ++l )
						{
							lk_base[k] += 
								( prob_ss_3 > prob_ss_5 
								  ? lerp( prob_ss_3, simple_adna::ss_mat[1<<k][1<<l], simple_adna::ds_mat[1<<k][1<<l] )
								  : lerp( prob_ss_3, simple_adna::ss_mat[complement(1<<k)][complement(1<<l)],
									  simple_adna::ds_mat[complement(1<<k)][complement(1<<l)] ) ) *
								( l == aln_i.base() ? l_mat : l_mismat ) ;
						}
					}

					// Translation of indices and calculation of het
					// likelihoods.  A bit irregular, so spelled out in
					// full.
					column->lk[0] *= lk_base[0] ; // A
					column->lk[1] *= lk_base[1] ; // C
					column->lk[2] *= lk_base[3] ; // G(!)
					column->lk[3] *= lk_base[2] ; // T(!)

					column->lk[4] *= lerp( 0.5, lk_base[0], lk_base[1] ) ; // AC
					column->lk[5] *= lerp( 0.5, lk_base[0], lk_base[3] ) ; // AG
					column->lk[6] *= lerp( 0.5, lk_base[0], lk_base[2] ) ; // AT
					column->lk[7] *= lerp( 0.5, lk_base[1], lk_base[3] ) ; // CG
					column->lk[8] *= lerp( 0.5, lk_base[1], lk_base[2] ) ; // CT
					column->lk[9] *= lerp( 0.5, lk_base[3], lk_base[2] ) ; // GT
				}

				++aln_i ;
				break ;

			case Hit::Delete:
				++column->gapped ;
				++aln_i ;
				++ref ;
				break ;

			default: 	// all others are unused
				throw "unexpected CIGAR operation" ;
				break ;
		}
		if( aln_i != aln_b ) ++column->crossed ;
	}
} 

} // namespace

