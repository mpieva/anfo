#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "anfo_common.h"
#include "queue.h"
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
//! network, from a central service.  The general course of action is
//! thus for the local supervisor:
//!
//! - forks the intermediate child,
//! - initializes ENet, 
//! - waits a bit (enet_host_service with long timeout),
//! - services ENet until no peers are left,
//! - exits.
//!
//! The intermediate child:
//!
//! - detaches from parent,
//! - forks the workers,
//! - exits.
//!
//! The workers:
//!
//! - initialize ENet,
//! - connect to local supervisor,
//! - connect to central server,
//! - acquire a configuration from the central server and initialize
//!   shared data structures,
//! - receive workloads from the central server and send results back,
//! - exit when either connection is dropped.
//!
//! This allows for some pretty nifty features:  We can reconfigure at
//! runtime just by supplying a new configuration.  We can also use a
//! side channel of ENet to send messages to ground control, that way we
//! can avoid littering SGE output everywhere...  Teardown of the
//! network is easy here:  We exit once a network connection disappears.
//! That way, SGE keeps a bit of control (we can still be killed).
//!
//! Since long alignments kill the ENet subsystem, we need to serve the
//! ENet host regularly.  Since messing with _every_ calculation to
//! include enet support is cumbersome, we'll use one thread for enet
//! comms and one for calculation.  Communication is done through two
//! queues of pointers to raw packets; the actual worker goes into the
//! background, by sending a null pointer, we can tell it to shut down.
//! It will respond by sending a null pointer itself.  We'll limit the
//! que size to 16, which is twice a smuch as we ever expect to get
//! sent.  If the incoming queu overflows, we will defend ourselves by
//! simply dropping the packet.

class Worker 
{
	private:
		int id_ ;
		Queue< ENetPacket*, 16 > incoming_ ;
		Queue< ENetPacket*, 16 > outgoing_ ;

		ENetHost *host_ ;
		ENetPeer *grampa_, *gcontrol_ ;
		pthread_t worker_thread_ ;
		std::auto_ptr< Mapper > mapper ;

	public:
		//! \brief actual worker thread
		void run_worker()
		{
			while( ENetPacket *packet = incoming_.dequeue() )
			{
				if( packet->dataLength && packet->data[0] == packet_config )
				{
					// got configuration 
					config::Config conf ;
					conf.ParseFromArray( packet->data+1, packet->dataLength-1 ) ;
					mapper.reset( new Mapper( conf ) ) ;
				}
				if( packet->dataLength && packet->data[0] == packet_read && mapper.get() )
				{
					// got sequence
					output::Read read ;
					read.ParseFromArray( packet->data+1, packet->dataLength-1 ) ;

					QSequence ps( read.sequence().c_str(), (const uint8_t*)read.quality().c_str(),
							read.seqid(), read.description() ) ;
					output::Result res ;
					std::deque< alignment_type > ol ;
					int pmax = mapper->index_sequence( ps, res, ol ) ;
					if( pmax != INT_MAX ) mapper->process_sequence( ps, pmax, ol, res ) ;

					ENetPacket *p = enet_packet_create( 0, 1+res.ByteSize(), ENET_PACKET_FLAG_RELIABLE ) ;
					p->data[0] = packet_result ;
					res.SerializeToArray( p->data+1, p->dataLength-1 ) ;
					outgoing_.enqueue( p ) ;
				}
				enet_packet_destroy( packet ) ;
			}
			outgoing_.enqueue(0) ;
		}

		static void *run_worker_wrap( void *p ) { ((Worker*)p)->run_worker() ; return 0 ; }

		//! \brief initialization including network setup
		//! \param id numerical id of this worker (used to identify it
		//!           in log messages)
		//! \param lport local port of supervisor
		//! \param remote_address address of ground control
		Worker( int id, enet_uint16 lport, const ENetAddress *remote_address )
			: id_(id), incoming_()
		{
			throw_if_negative( enet_initialize(), "initializing ENet" ) ;
			host_ = throw_errno_if_null(
					enet_host_create( 0, 2, 0, 0 ), "creating ENet host" ) ;

			ENetAddress master_addr = { ENET_HOST_ANY, lport } ;
			enet_address_set_host( &master_addr, "localhost" ) ;
			grampa_ = throw_errno_if_null(
					enet_host_connect( host_, &master_addr, 1 ), "connecting to master" ) ;

			gcontrol_ = throw_errno_if_null(
					enet_host_connect( host_, remote_address, 1 ), "connecting to ground control" ) ;

			throw_if_not_null( pthread_create( &worker_thread_, 0, run_worker_wrap, this ), "creating worker thread" ) ;
		}
		
		//! \brief tries to shutdown
		//! Whatever happened, we will exit anyway.  At this point, we try
		//! to disconnect everything, whether that works or not is
		//! irrelevant.
		~Worker() {
			enet_peer_disconnect_later( grampa_, 0 ) ;
			enet_peer_disconnect_later( gcontrol_, 0 ) ;
			enet_host_flush( host_ ) ;
			enet_host_destroy( host_ ) ;
			enet_deinitialize() ;
		}

		//! \brief serves ENet events in a loop
		//! ENet dictates that everything has to be event driven:
		//! - If a packet comes in, we enqueue it.
		//! - If someone disconnects, we return false (thus causing an early exit).
		//! - If nothing happens for... say... 60 seconds, we also exit (in
		//!   this case something is _very_ wrong)
		//! \param timeout timeout in milliseconds if we should wait for
		//!                an event
		//! \return return 0 if processing shall continue, else an
		//!                exit code xored with INT_MIN

		int run()
		{
			int timeout = 1000 ;
			for(;;)
			{
				ENetEvent event ;
				throw_if_negative( enet_host_service( host_, &event, timeout ), "servicing host" ) ;

				// is a timeout an error? may be dangerous...
				// if( event.type == ENET_EVENT_TYPE_NONE ) {} // return 0 ; // timeout ? 1 ^ INT_MIN : 0 ;
				if( event.type == ENET_EVENT_TYPE_DISCONNECT && event.peer == gcontrol_ ) break ; // return INT_MIN ;
				if( event.type == ENET_EVENT_TYPE_DISCONNECT && event.peer == grampa_ ) break ; //  return 1 ^ INT_MIN ;
				if( event.type == ENET_EVENT_TYPE_RECEIVE && event.peer == gcontrol_ )
				{
					// just enqueue it, but if it doesn't fit, destroy
					// it in self-defense
					if( !incoming_.try_enqueue( event.packet ) ) enet_packet_destroy( event.packet ) ;
				}
				else if( event.type == ENET_EVENT_TYPE_RECEIVE ) 
					enet_packet_destroy( event.packet ) ;

				ENetPacket *p ;
				if( outgoing_.try_dequeue( p ) )
				{
					enet_peer_send( gcontrol_, 0, p ) ;
					timeout = 0 ;
				}
				else timeout = 1000 ;
			}

			exit_with = 1 ; // in case alignments are still running
			incoming_.enqueue( 0 ) ; // signal worker to exit
			pthread_join( worker_thread_, 0 ) ;
			return 0 ;
		}
} ;

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
			case ENET_EVENT_TYPE_CONNECT:    ++num_peers ; break ;
			case ENET_EVENT_TYPE_DISCONNECT: --num_peers ; break ;
			case ENET_EVENT_TYPE_RECEIVE:    enet_packet_destroy( event.packet ) ; break ;
			case ENET_EVENT_TYPE_NONE:       break ;
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

	if( throw_errno_if_minus1( fork(), "forking first child" ) ) return run_master( lport, nworkers ) ;
	else {
		// middle child: detach from parent (is chdir() and umask() needed, and if so, why?)
		throw_errno_if_minus1( setsid(), "creating session" ) ;
		throw_errno_if_minus1( chdir("/"), "changing working directory" ) ;
		umask( 0 ) ;
		
		throw_errno_if_null( freopen( "/dev/null", "r", stdin), "redirecting stdin" ) ;
		throw_errno_if_null( freopen( "/dev/null", "w", stdout), "redirecting stdout" ) ;
		throw_errno_if_null( freopen( "/dev/null", "w", stderr), "redirecting stderr" ) ;

		for( unsigned i = 0 ; i != nworkers ; ++i )
			if( !throw_errno_if_minus1( fork(), "forking worker" ) )
				return Worker( i, lport, &remote_address ).run() ;
		return 0 ;
	}
}
