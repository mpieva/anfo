#include "ducttape.h"

namespace streams {

void DuctTaper::put_header( const Header& h ) 
{
	state_ = need_input ;
	// XXX: probably need to store part of the config
}

void DuctTaper::put_footer( const Footer& f ) 
{
	state_ = end_of_stream ;
	// XXX: need to flush contig
}

void DuctTaper::put_result( const Result& r )
{
	// XXX need to scan
} 

} // namespace

