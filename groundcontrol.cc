
#include "compress_stream.h"
#include "conffile.h"
#include "outputfile.h"
#include "sequence.h"

#include "config.pb.h"
#include "output.pb.h"

#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>

#include <enet/enet.h>
#include <iostream>
#include <memory>

//! \todo Teardown of workers is not handled gracefully.  When we
//!       decide to go down, we should explicitly disconnect every
//!       worker in turn.
//! \todo Unannounced disconnection of workers is not handled correctly.
//!       We need to recover messages in flight and distribute then to
//!       other workers.
void qs_to_read( const QSequence& qs, output::Read& rd ) 
{
	rd.set_seqid( qs.get_name() ) ;
	rd.set_description( qs.get_descr() ) ;
	rd.clear_sequence() ;
	rd.clear_quality() ;
	
	for( size_t i = 0 ; i != qs.length() ; ++i )
	{
		rd.mutable_sequence()->push_back( from_ambicode( qs[i].ambicode ) ) ;
		rd.mutable_quality( )->push_back(                qs[i].qscore     ) ;
	}
}

int main_( int argc, const char**argv )
{
	// for simplicity: arguments are port number and config file
	// input comes from stdin (FASTA or FASTQ), output goes to stdout
	// (ANFO, compressed)
	
	if( argc != 3 ) {
		std::clog << "Usage: " << argv[0] << " <lport> <config-file>\n"
			"  where lport is the local port to use, config-file is the configuration\n"
			"  file, input comes from stdin and output goes to stdout.  Debug\n"
			"  messages go to stderr." << std::endl ;
		return 1 ;
	}

	// Outline: get ENet up and running, read config, set up streams.
	// The rest is event driven:
	//
	// - if someone connects, send him the config and preload with
	//   sequences,
	// - XXX: if someone disconnects, enqueue his workload,
	// - if we receive a result, write it out and send a replacement,
	// - if we're out of sequences, send quit instead,
	// - on timeout, simply bail out.
	
	config::Config conf = get_default_config( argv[2] ) ;
	std::clog << "got configuration" << std::endl ;

	google::protobuf::io::FileOutputStream fos( 1 ) ;
	std::auto_ptr< google::protobuf::io::ZeroCopyOutputStream > zos( compress_fast( &fos ) ) ;

	output::Header ohdr ;
	google::protobuf::io::CodedOutputStream cos( zos.get() ) ;
	cos.WriteRaw( "ANFO", 4 ) ; // signature
	*ohdr.mutable_config() = conf ;
	ohdr.set_version( PACKAGE_VERSION ) ;
	write_delimited_message( cos, 1, ohdr ) ;
	std::clog << "wrote header" << std::endl ;

	google::protobuf::io::FileInputStream raw_inp( 0 ) ;
	std::auto_ptr<google::protobuf::io::ZeroCopyInputStream> inp( decompress( &raw_inp ) ) ;

	ENetAddress comsat_address = { ENET_HOST_ANY, atoi( argv[1] ) } ;
	throw_if_negative( enet_initialize(), "initializing ENet" ) ;
	ENetHost* comsat = throw_errno_if_null( enet_host_create( &comsat_address, 1024, 0, 0 ), "creating comsat" ) ;
	std::clog << "comsat deployed" << std::endl ;

	output::Footer ofoot ;
	ofoot.set_exit_code( 0 ) ;

	int num_peers = 0, reads_in_flight = 0, reads_in = 0, reads_out = 0 ;
	do {
		ENetEvent event ;
		throw_if_negative( enet_host_service( comsat, &event, 30000 ), "servicing host" ) ;
		if( event.type == ENET_EVENT_TYPE_CONNECT ) 
		{
			std::clog << "\r\e[Kcomsat: worker connected, sending config and reads" << std::endl ;
			++num_peers ;

			ENetPacket *p = enet_packet_create( 0, 1+conf.ByteSize(), ENET_PACKET_FLAG_RELIABLE ) ;
			p->data[0] = packet_config ;
			conf.SerializeToArray( p->data+1, p->dataLength-1 ) ;
			enet_peer_send( event.peer, 0, p ) ;

			for( int i = 0 ; i != 8 ; ++i )
			{
				QSequence qs ;
				if( read_fastq( inp.get(), qs ) )
				{
					output::Read rd ;
					qs_to_read( qs, rd ) ;

					ENetPacket *p = enet_packet_create( 0, 1+rd.ByteSize(), ENET_PACKET_FLAG_RELIABLE ) ;
					p->data[0] = packet_read ;
					rd.SerializeToArray( p->data+1, p->dataLength-1 ) ;
					enet_peer_send( event.peer, 0, p ) ;
					++reads_in_flight ;
					++reads_out ;
					std::clog << "\r\e[Ksent message " << qs.get_name() << "; " 
						<< num_peers << " workers, " << reads_in_flight << " reads in flight." 
						<< std::endl ;
				}
			}
		}
		else if( event.type == ENET_EVENT_TYPE_DISCONNECT ) 
		{
			std::clog << "\r\e[Kcomsat: worker disconnected (this may not be good)" << std::endl ;
			--num_peers ;
		}
		else if( event.type == ENET_EVENT_TYPE_RECEIVE )
		{
			if( event.packet->dataLength && event.packet->data[0] == packet_result )
			{
				output::Result res ;
				res.ParseFromArray( event.packet->data+1, event.packet->dataLength-1 ) ;
				write_delimited_message( cos, 2, res ) ;
				--reads_in_flight ;
				++reads_in ;
				std::clog << "\r\e[Kreceived message; " 
					<< num_peers << " workers, " << reads_in_flight << " reads in flight, " 
					<< reads_in << " in, " << reads_out << " out." << std::endl ;

				QSequence qs ;
				if( read_fastq( inp.get(), qs ) )
				{
					output::Read rd ;
					qs_to_read( qs, rd ) ;
					ENetPacket *p = enet_packet_create( 0, 1+rd.ByteSize(), ENET_PACKET_FLAG_RELIABLE ) ;
					p->data[0] = packet_read ;
					rd.SerializeToArray( p->data+1, p->dataLength-1 ) ;
					enet_peer_send( event.peer, 0, p ) ;
					++reads_in_flight ;
					++reads_out ;

					std::clog << "\r\e[Ksent message " << qs.get_name() << "; " 
						<< num_peers << " workers, " << reads_in_flight << " reads in flight, " 
						<< reads_in << " in, " << reads_out << " out." << std::endl ;
				}
			}
			enet_packet_destroy( event.packet ) ;
		}
		else if( event.type == ENET_EVENT_TYPE_NONE ) 
		{
			ofoot.set_exit_code( 1 ) ;
			break ;
		}
	} while( num_peers && reads_in_flight ) ;

	write_delimited_message( cos, 3, ofoot ) ;

	enet_host_destroy( comsat ) ;
	enet_deinitialize() ;
	return ofoot.exit_code() ;
}



