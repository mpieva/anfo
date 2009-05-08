#include "config.h"

#include "anfo_common.h"
#include "util.h"

#include "config.pb.h"
#include "output.pb.h"

#include <enet/enet.h>

#include <deque>
#include <iostream>

#include <signal.h>

#if HAVE_UNISTD_H
#include <unistd.h>
#endif

//! \page anfo_enet ANFO with networking support
//!
//! Here's the Outline:  Since the Grid Engine will get in the way, we
//! need to run in separate processes that have to communicate somehow,
//! and that's naturally the network.
//!
//! So, any sge process will assume it has access to ~16GB RAM and four
//! nodes.  The master process daemon-forks off four workers, then waits
//! for their completion.  While we're at it, workload is taken from the
//! network, from a central service.  In summary, the master:
//!
//! - initializes ENet, 
//! - acquires a configuration,
//! - forks once,
//! - waits a bit,
//! - waits for all peers to disconnect,
//! - exits.
//!
//! The first child:
//!
//! - shuts down ENet,
//! - initializes shared data structures,
//! - detaches from parent,
//! - forks the workers,
//! - exits.
//!
//! The workers:
//!
//! - initialize ENet,
//! - connect to grandparent,
//! - connect to central server,
//! - exit on command or when either connection is dropped.

int run_worker( enet_uint16 lport, const ENetAddress *remote_address )
{
	throw_if_negative( enet_initialize(), "initializing ENet" ) ;
	ENetHost* host = throw_errno_if_null(
			enet_host_create( 0, 2, 0, 0 ), "creating ENet host" ) ;

	ENetAddress master_addr = { ENET_HOST_ANY, lport } ;
	enet_address_set_host( &master_addr, "localhost" ) ;
	ENetPeer* grampa = throw_errno_if_null(
			enet_host_connect( host, &master_addr, 1 ), "connecting to master" ) ;

	ENetPeer* gcontrol = throw_errno_if_null(
			enet_host_connect( host, remote_address, 1 ), "connecting to ground control" ) ;

	// from here on, everything has to be event driven...
	// - As soon as a configuration structure comes in, we init ANFO.
	// - If a sequence comes in, we align it...
	// - then send out the result.
	// - If someone disconnects, we exit (we just lost our purpose in
	//   life).
	// - If nothing happens for... say... 60 seconds, we also exit (in
	//   this case something is _very_ wrong)

	std::auto_ptr< Mapper > mapper ;
	for(;;)
	{
		ENetEvent event ;
		throw_if_negative( enet_host_service( host, &event, 60000 ), "servicing host" ) ;
		if( event.type == ENET_EVENT_TYPE_NONE ) break ;
		if( event.type == ENET_EVENT_TYPE_DISCONNECT && event.peer == gcontrol ) break ;
		if( event.type == ENET_EVENT_TYPE_DISCONNECT && event.peer == grampa ) break ;
		if( event.type == ENET_EVENT_TYPE_RECEIVE && event.peer == gcontrol )
		{
			// okay, got a packet.  we'll switch on the first byte to
			// see what it is
			if( event.packet->dataLength && event.packet->data[0] == packet_quit )
			{
				std::clog << "worker exiting on quit message" << std::endl ;
				enet_packet_destroy( event.packet ) ;
				break ;
			}
			if( event.packet->dataLength && event.packet->data[0] == packet_config )
			{
				// got configuration 
				config::Config conf ;
				conf.ParseFromArray( event.packet->data+1, event.packet->dataLength-1 ) ;
				mapper.reset( new Mapper( conf ) ) ;
			}
			if( event.packet->dataLength && event.packet->data[0] == packet_read && mapper.get() )
			{
				// got sequence
				output::Read read ;
				read.ParseFromArray( event.packet->data+1, event.packet->dataLength-1 ) ;

				QSequence ps( read.sequence().c_str(), (const uint8_t*)read.quality().c_str(),
						read.seqid(), read.description() ) ;
				output::Result res ;
				std::deque< alignment_type > ol ;
				int pmax = mapper->index_sequence( ps, res, ol ) ;
				if( pmax != INT_MAX ) mapper->process_sequence( ps, pmax, ol, res ) ;

				ENetPacket *p = enet_packet_create( 0, 1+res.ByteSize(), ENET_PACKET_FLAG_RELIABLE ) ;
				p->data[0] = packet_result ;
				res.SerializeToArray( p->data+1, p->dataLength-1 ) ;
				enet_peer_send( gcontrol, 0, p ) ;
			}
			enet_packet_destroy( event.packet ) ;
		}
	}
	std::clog << "worker exiting" << std::endl ;

	// uncontrolled teardown; ground control has to deal with
	// unannounced crashes anyway, so this is already more nice than
	// necessary
	enet_peer_disconnect( grampa, 0 ) ;
	enet_peer_disconnect( gcontrol, 0 ) ;
	enet_host_flush( host ) ;

	enet_host_destroy( host ) ;
	enet_deinitialize() ;
	return 0 ;
}

//! \brief runs the middle child in a daemon fork
//! Here we just detach from the parent, fork workers and exit.
//! \param lport local port of master process
//! \param nworkers number of workers to fork
//! \param remote_address network adress of ground control process
//! \return 0 on success
int run_middle_child( enet_uint16 lport, unsigned nworkers, const ENetAddress *remote_address )
{
	// detach from parent (is chdir() and umask() needed, and if so, why?)
	throw_errno_if_minus1( setsid(), "creating session" ) ;
	// chdir("/");
	// umask(0);

	for( unsigned i = 0 ; i != nworkers ; ++i )
		if( !throw_errno_if_minus1( fork(), "forking worker" ) )
			return run_worker( lport, remote_address ) ;
	return 0 ;
}

//! \brief runs local master task
//! The master doesn't have to do anything, it just needs to be there in
//! case the grid engine kills it.  It will exit once all workers are
//! gone, so we simply count connections here.  We also if no one
//! connected after the first minute.
//! \param lport local port to bind to
//! \param nworkers number of worker connections to expect
//! \return 0 on success
int run_master( enet_uint16 lport, unsigned nworkers )
{
	signal( SIGCHLD, SIG_IGN ) ;

	throw_if_negative( enet_initialize(), "initializing ENet" ) ;
	ENetAddress local_address = { ENET_HOST_ANY, lport } ;
	ENetHost* host = throw_errno_if_null( enet_host_create( &local_address, nworkers, 0, 0 ), "creating ENet host" ) ;

	int num_peers = 0 ;
	do {
		ENetEvent event ;
		throw_if_negative( enet_host_service( host, &event, 30000 ), "servicing host" ) ;
		switch( event.type )
		{
			case ENET_EVENT_TYPE_CONNECT:    
				std::clog << "local master: worker connected" << std::endl ;
				++num_peers ;
				break ;
			case ENET_EVENT_TYPE_DISCONNECT:
				std::clog << "local master: worker disconnected" << std::endl ;
				--num_peers ;
				break ;
			case ENET_EVENT_TYPE_RECEIVE:
				enet_packet_destroy( event.packet ) ;
				break ;
			case ENET_EVENT_TYPE_NONE:
				break ;
		}
	} while( num_peers ) ;

	enet_host_destroy( host ) ;
	enet_deinitialize() ;
	return 0 ;
}

int main_( int argc, const char**argv )
{
	if( argc != 5 ) {
		std::clog << "Usage: " << argv[0] << " <lport> <nworkers> <rhost> <rport>\n"
			"  where lport is the local port to use, nworkers is the number of\n"
			"  workers to fork, rhost is the centrol controlling host and rhost is\n"
			"  the port to connect to on rhost." << std::endl ;
		return 1 ;
	}

	ENetAddress remote_address = { ENET_HOST_ANY, atoi( argv[4] ) } ;
	enet_address_set_host( &remote_address, argv[3] ) ;
	unsigned nworkers = atoi( argv[2] ) ;
	enet_uint16 lport = atoi( argv[1] ) ;

	return throw_errno_if_minus1( fork(), "forking first child" )
		? run_master( lport, nworkers ) 
		: run_middle_child( lport, nworkers, &remote_address ) ;
}
