//    Copyright 2009 Udo Stenzel
//    This file is part of ANFO
//
//    ANFO is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    Anfo is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with Anfo.  If not, see <http://www.gnu.org/licenses/>.

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include "index.h"
#include "util.h"

#include <cmath>
#include <new>

#include <fcntl.h>
#include <glob.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>

using namespace config ;
using namespace std ; 

static inline void cleanup_genome_name( string& gn )
{
	size_t gnl = gn.length() ;
	for( size_t i = 0 ; i != gnl ; ++i ) gn[i] = tolower( gn[i] ) ;
	if( gnl > 4 && gn.substr( gnl-4 ) == ".dna" ) gn = gn.substr( 0, gnl-4 ) ;
}

CompactGenome::CompactGenome( const std::string &name )
	: base_(), file_size_(0), length_(0), fd_(-1), contig_map_(), posn_map_(), g_(), refcount_(0)
{
	try
	{
		fd_ = throw_errno_if_minus1( open( name.c_str(), O_RDONLY ), "opening", name.c_str() ) ;

		struct stat the_stat ;
		throw_errno_if_minus1( fstat( fd_, &the_stat ), "statting", name.c_str() ) ;
		file_size_ = the_stat.st_size ;
		void *p = Metagenome::mmap( 0, file_size_, &fd_, 0 ) ;
		throw_errno_if_minus1( p, "mmapping", name.c_str() ) ;
		base_.assign( (uint8_t*)p ) ;

		if( ((uint32_t const*)p)[0] != signature ) 
			throw name + string(" does not have 'DNA1' signature") ;

		uint32_t meta_off = ((uint32_t const*)p)[1] ;
		uint32_t meta_len = ((uint32_t const*)p)[2] ;

		if( !g_.ParseFromArray( (const char*)p + meta_off, meta_len ) )
			throw "error parsing meta information" ;

		cleanup_genome_name( *g_.mutable_name() ) ;
		length_  = meta_off ;

		for( int i = 0 ; i != g_.sequence_size() ; ++i )
		{
			const Sequence &s = g_.sequence(i) ;
			posn_map_[ s.name() ].first = i ;
			PosnMap1& pm = posn_map_[ s.name() ].second ;
			for( int j = 0 ; j != s.contig_size() ; ++j )
			{
				const Contig &c = s.contig(j) ;
				contig_map_[ c.offset() ] = std::make_pair( i, j ) ;
				pm[ c.range_start() ] = j ;
			}
		}
		contig_map_[ g_.total_size() ] = make_pair( -1, -1 ) ;
	}
	catch(...) 
	{
		if( base_ ) munmap( (void*)base_.unsafe_ptr(), file_size_ ) ;
		if( fd_ != -1 ) close( fd_ ) ;
		throw ;
	}
}

CompactGenome::~CompactGenome()
{
	if( base_ ) munmap( (void*)base_.unsafe_ptr(), file_size_ ) ;
	if( fd_ != -1 ) close( fd_ ) ;
}


void CompactGenome::report( uint32_t o, uint32_t l, const char* msg ) {
	if( msg )
		clog << '\r' << msg << ": at offset " << o << " (of " << 2*l
			<< ", " << round((50.0*o)/l) << "%)\e[K" << flush ;
}

FixedIndex::FixedIndex( const std::string& name )
	: base(0), secondary(0), first_level_len(0), length(0), fd_(-1), refcount_(0)
{
	void *p = 0 ;
	try 
	{
		fd_ = path_open( name, "idx", "ANFO_PATH" ) ;

		struct stat the_stat ;
		throw_errno_if_minus1( fstat( fd_, &the_stat ), "statting", name.c_str() ) ;
		length = the_stat.st_size ;
		p = Metagenome::mmap( 0, length, &fd_, 0 ) ;
		throw_errno_if_minus1( p, "mmapping", name.c_str() ) ;
		p_ = p ;

		if( *(const uint32_t*)p != signature )
			throw name + " does not have 'IDX3' signature" ;

		uint32_t meta_len = ((const uint32_t*)p)[1] ;
		if( !meta_.ParseFromArray( (const char*)p + 8, meta_len ) )
			throw "error parsing meta information" ;

		cleanup_genome_name( *meta_.mutable_genome_name() ) ;
		first_level_len = (1 << (2 * meta_.wordsize())) + 1 ;
		base = (const uint32_t*)( (const char*)p + 8 + ((3+meta_len) & ~3) ) ;
		secondary = base + first_level_len ; 
	}
	catch(...)
	{
		if( p_ ) munmap( (void*)p_, length ) ;
		if( fd_ != -1 ) close( fd_ ) ;
		throw ;
	}
}

//! \brief directly looks up an oligo
//! A single oligo is looked up and results are appended to a vector.
//! Seeds that occur too often can be ignored, but for statistical
//! purposes, they are always counted in the result.
//!
//! \note Setting both stride and cutoff will cause unexpected
//!       interactions, but since it's only intended for debugging, I
//!       won't fix it.
//!
//! \param o oligo to be looked up.
//! \param v receiver for the resulting seeds
//! \param cutoff discard oligos more frequent than this
//! \param offs offset value to be placed in the seeds
//! \param num_useless filled in with the number of seeds ignored due to
//!                    being too frequent
//! \return number of seeds found, including repetitive ones

unsigned FixedIndex::lookup1( Oligo o, PreSeeds& v, const FixedIndex::LookupParams& p, int32_t offs, int *num_useless ) const 
{
	Oligo o_min =     o   << ( 2*( metadata().wordsize() - p.wordsize ) ) ;
	Oligo o_max = ( (o+1) << ( 2*( metadata().wordsize() - p.wordsize ) ) ) ;
	if( o >= first_level_len )
		throw "oligo number out of index' range" ;

	uint32_t n = base[o_max] - base[o_min] ;
	if( n < p.cutoff ) v.post( secondary + base[o_min], secondary + base[o_max], offs, p.wordsize, p.stride ) ;
	else if( num_useless ) *num_useless += n ;
	return n ;
} 

//! \brief looks up an oligo with up to one mismatch
//! The oligo itself will first be looked for, then all its
//! one-substituion variants are generated (by successively xor'ing
//! every position with 01, 10 and 11) and looked up, too.
//!
//! \param o oligo to be looked up.
//! \param v receiver for the resulting seeds
//! \param cutoff discard oligos more frequent than this
//! \param offs offset value to be placed in the seeds
//! \param num_useless filled in with the number of seeds ignored due to
//!                    being too frequent
//! \return number of seeds found, including repetitive ones

unsigned FixedIndex::lookup1m(
		Oligo o,
		PreSeeds& v,
		const FixedIndex::LookupParams& p,
		int32_t offs,
		int *num_useless ) const 
{
	return lookup1m( o, v, p, offs, num_useless, p.wordsize ) ;
}

unsigned FixedIndex::lookup1m(
		Oligo o,
		PreSeeds& v,
		const FixedIndex::LookupParams& p,
		int32_t offs,
		int *num_useless,
		size_t imax ) const 
{
	unsigned total = lookup1( o, v, p, offs, num_useless ) ;
	Oligo m1 = 1, m2 = 2, m3 = 3 ;
	for( size_t i = 0 ; i != imax ; ++i )
	{
		total += lookup1( o ^ m1, v, p, offs, num_useless ) ;
		total += lookup1( o ^ m2, v, p, offs, num_useless ) ;
		total += lookup1( o ^ m3, v, p, offs, num_useless ) ;
		m1 <<= 2 ;
		m2 <<= 2 ;
		m3 <<= 2 ;
	}
	return total ;
} 

unsigned FixedIndex::lookup2m(
		Oligo o,
		PreSeeds& v,
		const FixedIndex::LookupParams& p,
		int32_t offs,
		int *num_useless ) const 
{
	unsigned total = lookup1m( o, v, p, offs, num_useless ) ;
	Oligo m1 = 1, m2 = 2, m3 = 3 ;
	for( size_t i = 0 ; i != p.wordsize ; ++i )
	{
		total += lookup1m( o ^ m1, v, p, offs, num_useless, i ) ;
		total += lookup1m( o ^ m2, v, p, offs, num_useless, i ) ;
		total += lookup1m( o ^ m3, v, p, offs, num_useless, i ) ;
		m1 <<= 2 ;
		m2 <<= 2 ;
		m3 <<= 2 ;
	}
	return total ;
} 

//! \brief looks up a whole sequence
//! The sequence is split into words as appropriate for the index, then
//! each one of them is looked up.  This method can be implemented for
//! any kind of index, whether based on fixed words or not.
//!
//! \param dna sequence to search
//! \param v receiver for resulting seeds
//! \param near_perfect set to one to search for seeds with one mismatch
//! \param num_useless is filled in with the number of seeds thrown away
//!                    due to being too frequent
//! \param cutoff disregard oligos more frequent than this
//! \return number of seeds found

unsigned FixedIndex::lookupS( const std::string& dna, PreSeeds& v,
		const FixedIndex::LookupParams &p, int *num_useless ) const
{
	Oligo o_f = 0, o_r = 0 ;
	Oligo mask = ~( ~0 << (p.wordsize * 2) ) ;

	int s = 2 * p.wordsize - 2 ;
	unsigned filled = 0 ;
	unsigned total = 0 ;

	for( int32_t offset = 0 ; offset != (int32_t)dna.length() ; ++offset )
	{
		o_f >>= 2 ;
		o_r = (o_r << 2) & mask ;
		++filled ;

		switch( dna[offset] )
		{
			case 'a': case 'A': o_f |= 0 << s ; o_r |= 2 ; break ;
			case 'c': case 'C': o_f |= 1 << s ; o_r |= 3 ; break ;
			case 't': case 'T': o_f |= 2 << s ; o_r |= 0 ; break ;
			case 'g': case 'G': o_f |= 3 << s ; o_r |= 1 ; break ;
			default: filled = 0 ; break ;
		}
		if( filled >= p.wordsize ) 
		{
			int32_t off_f = offset - p.wordsize + 2 ;
			int32_t off_r = - offset - 1 ;

			if( p.allow_mismatches >= 2 )
				total += lookup2m( o_f, v, p, off_f, num_useless ) 
					   + lookup2m( o_r, v, p, off_r, num_useless ) ;
			else if( p.allow_mismatches )
				total += lookup1m( o_f, v, p, off_f, num_useless ) 
					   + lookup1m( o_r, v, p, off_r, num_useless ) ;
			else
				total += lookup1( o_f, v, p, off_f, num_useless ) 
					   + lookup1( o_r, v, p, off_r, num_useless ) ;
		}
	}
	return total ;
}


std::map< std::string, FixedIndex* > MetaIndex::map_ ;

const config::Sequence *CompactGenome::translate_back( DnaP pos, uint32_t& offset ) const 
{
	ContigMap::const_iterator high = contig_map_.upper_bound( pos.abs() - base_ ) ;
	if( pos.abs() < base_ ) return 0 ; // before start
	if( high == contig_map_.end() ) return 0 ; // after end
	--high ;

	const config::Sequence &sequ = g_.sequence( high->second.first ) ;
	const config::Contig &contig = sequ.contig( high->second.second ) ;

	offset = pos.abs() - base_ - contig.offset() + contig.range_start() ;
	return &sequ ;
}

DnaP CompactGenome::find_pos( const std::string& seq, uint32_t pos ) const
{
	PosnMap::const_iterator a = posn_map_.find( seq ) ;
	if( a != posn_map_.end() && !a->second.second.empty() ) {
		PosnMap1::const_iterator b = a->second.second.upper_bound( pos ) ;
		--b ;
		
		const config::Contig& c = g_.sequence( a->second.first ).contig( b->second ) ;
		if( c.range_start() <= pos && c.range_end() >= pos )
			return base_ + c.offset() + (pos - c.range_start()) ;
	}
	return DnaP(0) ;
}

Metagenome Metagenome::the_metagenome( getenv("ANFO_PATH") ) ;

void Metagenome::make_room()
{
	// make room by forgetting about some genome
	int nephemeral = 0 ;
	for( Genomes::iterator gi = the_metagenome.genomes.begin(), ge = the_metagenome.genomes.end() ; gi != ge ; ++gi )
		if( gi->second->refcount_ == 0 ) ++nephemeral ;

	throw_errno_if_null( nephemeral, "mmap()ing" ) ;

	Genomes::iterator gi = the_metagenome.genomes.begin() ;
	for( int i = random() % nephemeral ; i || gi->second->refcount_ ; ++gi ) if( gi->second->refcount_ == 0 ) --i ;

	for( SeqMap::iterator i = the_metagenome.seq_map.begin(), ie = the_metagenome.seq_map.end() ; i != ie ; ++i )
	{
		for( SeqMap1::iterator j = i->second.begin(), je = i->second.end() ; j != je ; )
		{
			SeqMap1::iterator k = j ; ++j ;
			if( k->second == gi->second ) i->second.erase( k ) ;
		}
	}

	console.output( Console::info, "Metagenome: forgot about " + gi->first ) ;
	delete gi->second ;
	the_metagenome.genomes.erase( gi ) ;
}

Metagenome::Metagenome( const char* p )
{
	if( p ) 
	{
		const char* q = p ;
		while( *p )
		{
			while( *q && *q != ':' ) ++q ;
			path.push_back( std::string( p, q ) ) ;

			if( *q ) ++q ;
			p = q ;
		}
	}
	std::set_new_handler( &Metagenome::make_room ) ;
}

//! \brief searches genome path for matching files
//! This functions looks for exact matches (supporting index files that
//! name the genome file) and for files starting with the genome name
//! and ending in ".dna" (to support lookup of named genomes).

glob_t Metagenome::glob_path( const std::string& genome ) 
{
	glob_t the_glob = { 0 } ;
	int glob_append = 0 ;
	for( std::list< std::string >::iterator j = the_metagenome.path.begin() ; j != the_metagenome.path.end() ; ++j )
	{
		std::string p = *j + "/" + genome ;
		if( !genome.empty() ) {
			glob( p.c_str(), GLOB_NOSORT | glob_append | GLOB_TILDE, 0, &the_glob ) ;
			glob_append = GLOB_APPEND ;
		}
		p += "*.dna" ;
		glob( p.c_str(), GLOB_NOSORT | glob_append | GLOB_TILDE, 0, &the_glob ) ;
	}
	return the_glob ;
}

GenomeHolder Metagenome::find_sequence( const std::string& genome, const std::string& seq )
{
	SeqMap1 &m = the_metagenome.seq_map[ genome ] ;

	// check if we got the sequence somewhere already
	SeqMap1::iterator i = m.find( seq ) ;
	if( i == m.end() )
	{
		glob_t the_glob = glob_path( genome ) ;

		// sequence missing.  load genomes until we find the sequence
		for( size_t j = 0 ; j != the_glob.gl_pathc && i == m.end() ; ++j )
		{
			try
			{
				const char *n = the_glob.gl_pathv[j] ; 
				Genomes::iterator gi = the_metagenome.genomes.find( n ) ;
				if( gi == the_metagenome.genomes.end() )
				{
					console.output( Console::info, "Metagenome: trying file " + std::string(n) ) ;

					CompactGenome* g = new CompactGenome( n ) ;
					the_metagenome.genomes[ n ] = g ;

					console.output( Console::info, "Metagenome: loaded (part of) genome " + g->name() ) ;

					SeqMap1 &m_ = the_metagenome.seq_map[ g->name() ] ;
					for( int k = 0 ; k != g->g_.sequence_size() ; ++k )
						m_[ g->g_.sequence(k).name() ] = g ;

					if( g->name() == genome ) i = m.find( seq ) ;
				}
			}
			catch( const std::string& e ) { perr( e ) ; }
			catch( const char *e ) { perr( e ) ; }
			catch( char *e ) { perr( e ) ; }
			catch( const Exception& e ) { perr( e ) ; }
			catch( const std::exception& e ) { perr( e.what() ) ; }
			catch( ... ) { perr( "Oh noes!" ) ; }
		}
		globfree( &the_glob ) ;
	}
	return i == m.end() ? GenomeHolder() : i->second ;
}

GenomeHolder Metagenome::find_genome( const std::string& genome )
{
	// would be nicer to check if we already have a suitable genome, but
	// this is called so rarely, I won't bother

	glob_t the_glob = glob_path( genome ) ;
	if( !the_glob.gl_pathc ) throw "no genome file found for " + genome ;

	const char *n = the_glob.gl_pathv[0] ; 
	Genomes::iterator gi = the_metagenome.genomes.find( n ) ;
	if( gi == the_metagenome.genomes.end() )
	{
		CompactGenome *g = new CompactGenome( n ) ;
		the_metagenome.genomes[ n ] = g ;
		globfree( &the_glob ) ;
		return g ;
	}
	globfree( &the_glob ) ;
	return gi->second ;
}

bool Metagenome::translate_to_genome_coords( DnaP pos, uint32_t &xpos, const config::Sequence** s_out, const config::Genome** g_out )
{
	for( Genomes::const_iterator g = the_metagenome.genomes.begin(), ge = the_metagenome.genomes.end() ; g != ge ; ++g )
	{
		if( const Sequence *sequ = g->second->translate_back( pos, xpos ) )
		{
			if( s_out ) *s_out = sequ ;
			if( g_out ) *g_out = &g->second->g_ ;
			return true ;
		}
	}
	return false ;
}

static int get_zero_device()
{
	static int fdz = throw_errno_if_minus1(
			open( "/dev/zero", O_RDWR ), "opening", "/dev/zero" ) ;
	return fdz ;
}

int Metagenome::nommap = false ;

void *Metagenome::mmap( void *start, size_t length, int *fd, off_t offset )
{
	for(;;)
	{
		if( nommap ) {
			void *p = ::mmap( 0, length, PROT_READ | PROT_WRITE,
							  MAP_PRIVATE, get_zero_device(), 0 ) ;
			if( p != (void*)(-1) ) {
				myread( *fd, p, length ) ;
				close( *fd ) ;
				*fd = -1 ;
				return p ;
			}
		} else {
			void *p = ::mmap( start, length, PROT_READ, 
					MAP_PRIVATE | MAP_NORESERVE, *fd, offset ) ;
			if( p != (void*)(-1) ) return p ;
		}
		make_room() ;
	}
}

