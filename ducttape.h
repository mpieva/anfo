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
		google::protobuf::io::FileOutputStream fos_ ;
		DeflateStream out_ ;
		const char* g_ ;

		// tracking of current reference sequence; we need to start a
		// new contig if one of these changes
		std::string cur_genome_ ;
		std::string cur_sequence_ ;

		// start of current contig on the reference; we're always on the
		// forward strand
		uint32_t contig_start_ ;

		// stop-gap: number of observed bases per position in contig
		// will soon be replaced by somthing that can be used to
		// calculate posterior probabilities
		std::vector< int > observed_[5] ;

		// current CIGAR line
		// This shouldn't ever contain a deletion, but inserts will
		// happen.  If cigar_ is empty, no contig has been started yet.
		std::vector< int > cigar_ ;

		void flush_contig() ;

	public:
		DuctTaper( const char* fn, const char* g ) 
			: fos_( throw_errno_if_minus1( creat( fn, 0666 ), "writing to", fn ) )
			, out_( &fos_ ), g_(g)
			{ fos_.SetCloseOnDelete(true) ; }

		DuctTaper( int fd, const char* g )
			: fos_( fd ), out_( &fos_ ), g_(g) {}

		virtual ~DuctTaper() {}

		virtual void put_header( const Header& ) ;
		virtual void put_result( const Result& ) ;
		virtual void put_footer( const Footer& ) ;
		virtual Result fetch_result() { state_ = need_input ; return res_ ; }
} ;

} // namespace

#endif
