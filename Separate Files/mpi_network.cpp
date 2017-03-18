#include "mpi_network.h"

MPI_Network::MPI_Network(): Network() {

}

MPI_Network::MPI_Network(MPI_Network const & source) {

	initialized = false;
	copy(source);
	initialized = true;

}

MPI_Network::MPI_Network(Network const & source): Network(source) {

}

MPI_Network::MPI_Network(sacNetwork const & source): Network(source) {
	
}