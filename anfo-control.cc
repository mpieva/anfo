
#include "compress_stream.h"
#include "conffile.h"
#include "sequence.h"
#include "stream.h"

#include "config.pb.h"
#include "output.pb.h"

#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>

#include <enet/enet.h>
#include <iostream>
#include <memory>

static const int max_timeouts = 8 ;

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
	// - if someone disconnects, enqueue his workload,
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
	streams::write_delimited_message( cos, 1, ohdr ) ;
	std::clog << "wrote header" << std::endl ;

	google::protobuf::io::FileInputStream raw_inp( 0 ) ;
	std::auto_ptr<google::protobuf::io::ZeroCopyInputStream> inp( decompress( &raw_inp ) ) ;

	ENetAddress comsat_address = { ENET_HOST_ANY, atoi( argv[1] ) } ;
	throw_if_negative( enet_initialize(), "initializing ENet" ) ;
	ENetHost* comsat = throw_errno_if_null( enet_host_create( &comsat_address, 1024, 0, 0 ), "creating comsat" ) ;
	std::clog << "comsat deployed" << std::endl ;

	output::Footer ofoot ;
	ofoot.set_exit_code( 0 ) ;

	int reads_in_flight = 0, reads_in = 0, reads_out = 0, timeouts = 0 ;
	std::deque< QSequence > incoming_queue ;

	typedef std::map< std::string, QSequence > SeqDir ;
	typedef std::map< ENetPeer*, SeqDir > PeerDir ;
	PeerDir peers ;

	do {
		ENetEvent event ;
		throw_if_negative( enet_host_service( comsat, &event, 30000 ), "servicing host" ) ;
		if( event.type == ENET_EVENT_TYPE_CONNECT ) 
		{
			std::clog << "\r\e[Kcomsat: worker connected, sending config and reads" << std::endl ;
			std::map< std::string, QSequence > &seqs = peers[ event.peer ] ;

			ENetPacket *p = enet_packet_create( 0, 1+conf.ByteSize(), ENET_PACKET_FLAG_RELIABLE ) ;
			p->data[0] = packet_config ;
			conf.SerializeToArray( p->data+1, p->dataLength-1 ) ;
			enet_peer_send( event.peer, 0, p ) ;

			for( int i = 0 ; i != 8 ; ++i )
			{
				QSequence qs ;
				if( !incoming_queue.empty() || read_fastq( inp.get(), qs ) )
				{
					if( !incoming_queue.empty() )
					{
						qs = incoming_queue.front() ;
						incoming_queue.pop_front() ;
					}

					output::Read rd ;
					qs_to_read( qs, rd ) ;

					ENetPacket *p = enet_packet_create( 0, 1+rd.ByteSize(), ENET_PACKET_FLAG_RELIABLE ) ;
					p->data[0] = packet_read ;
					rd.SerializeToArray( p->data+1, p->dataLength-1 ) ;
					enet_peer_send( event.peer, 0, p ) ;
					++reads_in_flight ;
					++reads_out ;
					seqs[ qs.get_name() ] = qs ;

					std::clog << "\r\e[Ksent message " << qs.get_name() << "; " 
						<< peers.size() << " workers, " << reads_in_flight << " reads in flight, " 
						<< seqs.size() << " on this peer." << std::endl ;
				}
			}
		}
		else if( event.type == ENET_EVENT_TYPE_DISCONNECT ) 
		{
			SeqDir sd = peers[ event.peer ] ; 
			std::clog << "\r\e[Kcomsat: worker disconnected (" ;
			if( sd.empty() ) std::clog << "not a problem" ;
			else std::clog << sd.size() << " sequences re-enqueued" ;
			std::clog << ')' << std::endl ;

			for( SeqDir::const_iterator l = sd.begin(), r = sd.end() ; l != r ; ++l )
				incoming_queue.push_back( l->second ) ;
			reads_in_flight -= sd.size() ;
			peers.erase( event.peer ) ;
		}
		else if( event.type == ENET_EVENT_TYPE_RECEIVE )
		{
			if( event.packet->dataLength && event.packet->data[0] == packet_result )
			{
				std::map< std::string, QSequence > &seqs = peers[ event.peer ] ;
				output::Result res ;
				res.ParseFromArray( event.packet->data+1, event.packet->dataLength-1 ) ;
				streams::write_delimited_message( cos, 2, res ) ;
				--reads_in_flight ;
				++reads_in ;
				seqs.erase( res.seqid() ) ;

				std::clog << "\r\e[Kreceived message; " 
					<< peers.size() << " workers, " << reads_in_flight << " reads in flight, " 
					<< reads_in << " in, " << reads_out << " out." << std::endl ;

				QSequence qs ;
				if( !incoming_queue.empty() || read_fastq( inp.get(), qs ) )
				{
					if( !incoming_queue.empty() )
					{
						qs = incoming_queue.front() ;
						incoming_queue.pop_front() ;
					}

					output::Read rd ;
					qs_to_read( qs, rd ) ;
					ENetPacket *p = enet_packet_create( 0, 1+rd.ByteSize(), ENET_PACKET_FLAG_RELIABLE ) ;
					p->data[0] = packet_read ;
					rd.SerializeToArray( p->data+1, p->dataLength-1 ) ;
					enet_peer_send( event.peer, 0, p ) ;
					++reads_in_flight ;
					++reads_out ;
					seqs[ qs.get_name() ] = qs ;

					std::clog << "\r\e[Ksent message " << qs.get_name() << "; " 
						<< peers.size() << " workers, " << reads_in_flight << " reads in flight, " 
						<< reads_in << " in, " << reads_out << " out, "
						<< seqs.size() << " on this peer." << std::endl ;
				}
			}
			enet_packet_destroy( event.packet ) ;
		}

		if( event.type == ENET_EVENT_TYPE_NONE ) 
		{
			std::clog << "Timeout, " << peers.size() << " workers still connected." << std::endl ;
			++timeouts ;
			if( timeouts == max_timeouts )
				throw peers.empty() ? "no workers available, giving up" 
					                : "network communication failed, aborting" ;
		}
		else timeouts = 0 ;
	} while( reads_in_flight || !incoming_queue.empty() ) ;

	streams::write_delimited_message( cos, 3, ofoot ) ;

	for( PeerDir::const_iterator l = peers.begin(), r = peers.end() ; l != r ; ++l )
		enet_peer_disconnect_later( l->first, 0 ) ;

	enet_host_flush( comsat ) ;
	enet_host_destroy( comsat ) ;
	enet_deinitialize() ;
	return ofoot.exit_code() ;
}



