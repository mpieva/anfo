//    Copyright 2009 Udo Stenzel
//    This file is part of ANFO
//
//    ANFO is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    Foobar is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

#include "config.h"

//! \page dump_params Parameter Extraction
//! dump-params expects exactly one file argument, which must be either
//! an ANFO output file or a textual configuration file.  The file is
//! opened, a configuration extracted, and if it contains the parameters
//! for an ancient DNA aligner, the aligner is configured and its
//! internal parameters (substitution matrices and associated penalties)
//! are dumped to stdout.

#include "align.h"
#include "conffile.h"
#include "stream.h"
#include "util.h"

#include "output.pb.h"
#include <iostream>
#include <memory>

int main_( int argc, const char** argv )
{
	if( argc != 2 ) throw "expected exactly one argument" ;

	config::Config conf ;
	try {
		std::auto_ptr< streams::Stream > af( streams::make_input_stream( argv[1] ) ) ;
		conf = af->fetch_header().config() ;
	} 
	catch( const streams::AnfoReader::ParseError& )
	{
		conf = parse_text_config( argv[1] ) ;
	}
	
	if( !conf.has_aligner() ) throw "couldn't find aligner configuration" ;
	simple_adna::configure( conf.aligner(), &std::cout ) ;
	return 0 ;
}

