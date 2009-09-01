#include "config.h"
#include "ducttape.h"

#include "align.h"
#include "index.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
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
//! likely base.  Calling ambiguity codes for heterozygous sites is not
//! attempted, as that would be pointless without applying a prior for
//! the rate of het sites.  Finally store the new contig as an \c
//! output::Result and update internal state.

void DuctTaper::flush_contig()
{
	static int external_to_internal_base[4] = { 0, 1, 3, 2 } ;
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
	report_( Console::info, ss1.str() ) ;
	
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
				rd.add_seen_bases( i->seen[ external_to_internal_base[j] ] ) ;
				if( i->lk[j] > i->lk[maxlk] ) maxlk = j ;
			}

			rd.add_depth( cov + i->seen[4] ) ;
			rd.mutable_sequence()->push_back( "ACGT"[maxlk] ) ;

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

		Logdom qual() const { return Logdom::from_phred( (uint8_t)( fwd_ ? qual_[0] : qual_[-1] ) ) ; }
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
//! The likelihood of the alignment with a given overhang length is
//! simply the alignment score, by dividing by the total likelihood, we
//! get a likelihood for that overhang length.
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

void DuctTaper::put_result( const Result& r )
{
	if( !has_hit_to( r, g_ ) ) return ;
	const Hit& h = hit_to( r, g_ ) ;

	Logdom rate_ss = Logdom::from_float( hdr_.config().aligner().rate_of_ss_deamination() ), 
		   rate_ds = Logdom::from_float( hdr_.config().aligner().rate_of_ds_deamination() ) ; 

	int mapq = h.has_diff_to_next() ? h.diff_to_next() : 254 ;
	mapq_accum_ += mapq*mapq ;

	if( cur_genome_ != h.genome_name()
				|| cur_sequence_ != h.sequence()
				|| h.start_pos() > contig_end_ )
	{
		flush_contig() ;
		cur_genome_ = h.genome_name() ;
		cur_sequence_ = h.sequence() ;
		contig_start_ = contig_end_ = h.start_pos() ;
	}
	++nreads_ ;

	if( h.start_pos() < contig_start_ ) 
		throw "ducttaping: input was not sorted" ;

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
				lk_ds_5 *= simple_adna::overhang_ext_penalty ;
				lk_ds_3 /= simple_adna::overhang_ext_penalty ;

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
					
#if 0
					// Original code, unrolled and simplified below to
					// avoid numerical difficulties.
					Logdom lk_base[4]
					for( int k = 0 ; k != 4 ; ++k )
					{
						lk_base[k] = Logdom::null() ;
						for( int l = 0 ; l != 4 ; ++l )
						{
							lk_base[k] +=
								( l == aln_i.base() ? l_mat : l_mismat ) *
								( prob_ss_5 > prob_ss_3 
								  ? lerp( prob_ss_5,
									  Logdom::from_float( ss_mat[k][l] ),
									  Logdom::from_float( ds_mat[k][l] ) ) 
								  : lerp( prob_ss_3,
									  Logdom::from_float( ss_mat[k^2][l^2] ),
									  Logdom::from_float( ds_mat[k^2][l^2] ) ) ) ;
						}
					}
#endif

					Logdom lk_base[4] = {
						( 0 == aln_i.base() ? l_mat : l_mismat ),			// k == 0
						( 1 == aln_i.base() ? l_mat : l_mismat ) *			// k == 1
							( 1 - lerp( prob_ss_5, rate_ss, rate_ds ) ) +	// l == 1

							( 2 == aln_i.base() ? l_mat : l_mismat ) *		// l == 2
							lerp( prob_ss_5, rate_ss, rate_ds ),

						( 2 == aln_i.base() ? l_mat : l_mismat ),			// k == 2
						
						( 0 == aln_i.base() ? l_mat : l_mismat ) *			// k == 3
							lerp( prob_ss_3, rate_ss, rate_ds ) +			// l == 0

							( 3 == aln_i.base() ? l_mat : l_mismat ) *		// l == 3
							( 1 - lerp( prob_ss_3, rate_ss, rate_ds ) )
					} ;
						
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

namespace {
	template< typename T > T *with_close( T* t )
	{ t->SetCloseOnDelete( true ) ; return t ; }
} ;

GlzWriter::GlzWriter( int fd ) : gos_( with_close( new google::protobuf::io::FileOutputStream( fd ) ) ) {}
GlzWriter::GlzWriter( const char* fn ) : gos_( with_close(
			new google::protobuf::io::FileOutputStream(
				throw_errno_if_minus1( creat( fn, 0666 ), "creating", fn ) ) ) ) {}

void GlzWriter::put_result( const Result& rr )
{
	static uint8_t dna_to_glf_base[] = { 0, 1, 2, 3, 8, 9, 10, 11, 4, 5, 6, 7, 12, 13, 14, 15 } ;

	const Read& r = rr.read() ;
	if( r.likelihoods_size() == 10 && has_hit_to( rr, 0 ) ) {
		chan_( Console::info, r.seqid() ) ;
		google::protobuf::io::CodedOutputStream c( &gos_ ) ;
		// Per chromosome
		// int   chrNameLen ;         /* includes terminating 0 */
		// char* chrName ;            /* chrNamelen chars including 0 */
		// int   chrLen ;
		c.WriteLittleEndian32( r.seqid().size()+1 ) ;
		c.WriteString( r.seqid() ) ;
		c.WriteTag( 0 ) ;
		c.WriteLittleEndian32( r.sequence().size() ) ;

		const Hit& h = hit_to( rr, 0 ) ;
		DnaP ref = Metagenome::find_sequence( h.genome_name(), h.sequence() ).find_pos( h.sequence(), h.start_pos() ) ;
		char buf[12] ;
		int i = 0 ;
		uint8_t min_lk ;

		for( AlnIter beg( r, h ), end( r, h, 1 ) ; beg != end ; ++beg )
		{
			switch( beg.cigar_op() )
			{
				case Hit::Match:
				case Hit::Mismatch:
					// Per base
					//   unsigned char ref:4, dummy:4 ; /* ref A=1,C=2,G=4,T=8,N=15 etc.  */
					//   unsigned char max_mapQ ;       /* maximum mapping quality */
					//   unsigned char lk[10] ;         /* log likelihood ratio, max 255 */
					//   unsigned min_lk:8,             /* minimum lk capped at 255
					//   depth:24 ;            			/* and the number of mapped reads */

					// I'm not sure what's the business about min_lk,
					// I'll read it literally and take the smallest
					// likelihood value.
					min_lk = r.likelihoods(0)[i] ;
					for( int j = 1 ; j != 10 ; ++j )
						if( min_lk > r.likelihoods(j)[i] ) 
							min_lk = r.likelihoods(j)[i] ; 

					buf[0] = dna_to_glf_base[ *ref ] ;
					buf[1] = h.has_diff_to_next() ? h.diff_to_next() : 254 ;
					for( int j = 0 ; j != 10 ; ++j ) buf[2+j] = (uint8_t)(r.likelihoods(j)[i]) - min_lk ;
					c.WriteRaw( buf, 12 ) ;
					c.WriteLittleEndian32( (unsigned(r.depth(i)) << 8) | min_lk ) ;
					++ref ;
					++i ;
					break ;

				case Hit::Delete:
					++ref ;
					break ;

				case Hit::Insert:
				case Hit::SoftClip:
					++i ;
					break ;

				default:
					break ;
			}
		}
	}
}

// Need to walk along any number of alignments (hg18, pt2, ..., nt),
// producing either gaps or bases at each position.  If *all* positions
// are gapped, the position is skipped.  Else produce symbols (gaps or
// base, either from genome or from majority sequence), summary
// information and likelihoods.
//
// XXX: for the time being, walk only one alignment
void ThreeAlnWriter::put_result( const Result& res )
{
	const Read& r = res.read() ;
	const Hit& h = hit_to( res, 0 ) ;
	DnaP ref = Metagenome::find_sequence( h.genome_name(), h.sequence() ).find_pos( h.sequence(), h.start_pos() ) ;

	std::stringstream ss ;
	ss << name_ << ": " << h.genome_name() << '/' << h.sequence() << '@' << h.start_pos() ; 
	chan_( Console::info, ss.str() ) ;

	out_ << '>' << r.seqid() << ' ' << r.description() << '\n' ;
	// XXX more header lines like this:
	// out_ << ';' << h.genome_name() << ' ' << h.sequence() << ' '
		// << ( h.aln_length() < 0 ? '-' : '+' ) << h.start_pos() << std::endl ;

	AlnIter aln( r, h ), aln_e( r, h, 1 ) ;
	int offs = 0 ;
	for( ; aln != aln_e ; ++aln ) 
	{
		switch( aln.cigar_op() )
		{
			case Hit::Match:
			case Hit::Mismatch:
			case Hit::Delete:
				out_ << from_ambicode( *ref ) ;
				break ;
			case Hit::Insert:
			case Hit::SoftClip:
				out_ << '-' ;
				break ;
			default:
				break ;
		}
		int ngaps = -std::accumulate(
				r.seen_bases().begin() + 4*offs,
				r.seen_bases().begin() + 4*offs +4,
				-r.depth(offs) ) ;

		switch( aln.cigar_op() )
		{
			case Hit::Match:
			case Hit::Mismatch:
			case Hit::Insert:
			case Hit::SoftClip:
				out_ << r.sequence()[offs] << std::setw(4) << (int)(uint8_t)r.quality()[offs]
					<< std::setw(5) << r.depth(offs) << std::setw(5) << ngaps ;

				out_ << "   " ;
				for( int i = 0 ; i != 4 ; ++i )
					out_ << std::setw(4) << r.seen_bases( i + 4*offs ) ;

				out_ << "   " ;
				for( int i = 0 ; i != 4 ; ++i )
					out_ << std::setw(4) << (int)(uint8_t)r.likelihoods(i)[offs] ;
				break ;

			case Hit::Delete:
				out_ << '-' ;
				break ;
			default:
				break ;
		}
		out_ << '\n' ;
		switch( aln.cigar_op() )
		{
			case Hit::Match:
			case Hit::Mismatch:
				++ref ;

			case Hit::Insert:
			case Hit::SoftClip:
				++offs ;
				break ;

			default:
				++ref ;
				break ;
		}
	}
}

} // namespace

