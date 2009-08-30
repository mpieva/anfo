#ifndef INCLUDED_DUCTTAPE_H
#define INCLUDED_DUCTTAPE_H

#include "stream.h"
#include "util.h"
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
		// Order is A,C,G,T,-.  Also tracks whether a column is inserted
		// relative to the reference.
		struct Acc {
			int seen[4] ;		// # of times A,C,G,T was seen
			int gapped ;		// # of times a gap was seen
			int crossed ;		// # of times a read overlapped the gap before this col.
			bool is_ins ;		// was this iserted rel. to reference?
			Logdom lk[10] ; 	// likelihoods for allele pairs

			Acc( bool ins = false, int gap = 0 ) : gapped(gap), crossed(0), is_ins(ins)
			{ seen[0] = seen[1] = seen[2] = seen[3] = 0 ; }

			bool pristine() const
			{ for( int i = 0 ; i != 5 ; ++i ) if( seen[i] ) return false ; return gapped == 0 ; }
		} ;
		typedef std::vector< Acc > Accs ;
		Accs observed_ ;

		// for silly statistics
		size_t nreads_ ; 

		// for MAPQ -- GLF wants MAPQ to be the RMS of the MAPQs that
		// went into a contig.
		int mapq_accum_ ;

		void flush_contig() ;

	public:
		DuctTaper( const char* g, int maxq ) 
			: g_(g), maxq_(maxq), contig_start_(0), contig_end_(0), nreads_(0), mapq_accum_(0) {}
		virtual ~DuctTaper() {}

		virtual void put_header( const Header& ) ;
		virtual void put_result( const Result& ) ;
		virtual void put_footer( const Footer& ) ;
		virtual Result fetch_result() ;
} ;

class GlzWriter : public Stream
{
	private:
		DeflateStream gos_ ;
		Chan chan_ ;

	public:
		GlzWriter( int fd ) ;
		GlzWriter( const char* fn ) ;
		virtual ~GlzWriter() {}
		virtual void put_result( const Result& ) ;
} ;

//! \brief writes consensuus in text format
//! Format as agreed upon internally:  
//! > "hg18" chromosome ±start
//! ; "pt2" chromosome ±start
//! [repeat as necessary?]
//! <hg18-base> <pt2-base> <nt-majority-base (A,C,G,T)>
//!   <Q-score> <depth> <#gaps> <#A> <#C> <#G> <#T> <4 likelihood-ratios>
//! [repeat as necessary]
//!
//! \todo reconstruction of the involved whole genome alignments is not
//!       yet possible
class ThreeAlnWriter : public Stream
{
	private:
		std::ofstream out_ ;
		std::string name_ ;
		Chan chan_ ;

	public:
		ThreeAlnWriter( const char* ) ;
		virtual ~ThreeAlnWriter() {} 
		virtual void put_result( const Result& ) ;
} ;

} // namespace

#endif
