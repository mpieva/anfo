#ifndef INCLUDED_DUCTTAPE_H
#define INCLUDED_DUCTTAPE_H

#include "stream.h"
#include <string>
#include <vector>

namespace streams {

//! \brief calls a genome wide consensus
//! Duct-taping something is almost the same as assembling it, it's just
//! not as durable, and that's what we do to a genome here.  Input is a
//! stream ordered by genome coordinate.  For every position covered at
//! least once, we calculate the likelihood for every possible
//! combination of two alleles.  In principle, those bases can be
//! written out as soon as reads stop covering them, however, since we
//! will write contigs, we can only write them when we no longer have
//! coverage (or reach the end of a contig in the reference?), as the
//! GLZ format wants to know the contig length in advance.
//!
//! \todo What to do at edges?  Sometimes reads will overlap the ends of
//!       contigs, and in principle, we might want to call a consensus
//!       there, too.  Probably nothing at all...
//! \todo What about indels?  How to deal with the shifting coordinates
//!       if our consensus actually contains indels?

class DuctTaper : public Stream
{
	private:
		const char* g_ ;	// will be 'assembled'
		int maxq_ ;			// clamp QScores to this

		// tracking of current reference sequence; we need to start a
		// new contig if one of these changes
		std::string cur_genome_ ;
		std::string cur_sequence_ ;

		// start of current contig on the reference; we're always on the
		// forward strand
		uint32_t contig_start_ ;
		uint32_t contig_end_ ;

		// Collects number of observed bases or gaps per position in
		// contig (to allow calling of indels, and because someone might
		// ask for it) and accumulated likelihoods for each base.
		// Accumulated quantity is \f$ P(\omega | x) * P(y->x) * P(y)
		// \f$, for we can easily calculate a probability from that.
		// Order is A,C,T,G,- 
		struct Acc {
			int seen[5] ;
			Logdom lk[4] ; 

			Acc() { seen[0] = seen[1] = seen[2] = seen[3] = seen[4] = 0 ; }
			bool pristine() const { for( int i = 0 ; i != 5 ; ++i ) if( seen[i] ) return false ; return true ; }
		} ;
		std::vector< Acc > observed_ ;

		// Marks columns that are insertions.  If empty, no contig has
		// been started yet.
		std::vector< bool > is_ins_ ;

		// for silly statistics
		size_t nreads_ ; 

		void flush_contig() ;

	public:
		DuctTaper( const char* g, int maxq ) 
			: g_(g), maxq_(maxq), contig_start_(0), contig_end_(0), nreads_(0) {}
		virtual ~DuctTaper() {}

		virtual void put_header( const Header& ) ;
		virtual void put_result( const Result& ) ;
		virtual void put_footer( const Footer& ) ;
		virtual Result fetch_result() ;
} ;

} // namespace

#endif
