#ifndef MPI_NETWORK_H
#define MPI_NETWORK_H
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <math.h>
#include <random>
#include <time.h>
#include <vector>
#include <unistd.h>
#include <algorithm>
#include <ctime>
#include <chrono>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <stddef.h>
#include "network.h"
#include "sac_network.h"
using namespace std;

#define NSYNC 100 // syncs plate_forces after NSYNC iterations

/**
@file mpi_network.h
\brief Wrapper around the Network and sacNetwork classes to allow for MPI implementation
*/

class MPI_Network : virtual public Network, public sacNetwork{

public:

	MPI_Network();
	MPI_Network(string&);
	MPI_Network(MPI_Network const & source);
	MPI_Network(Network const & source);
	MPI_Network(sacNetwork const & source);
	~MPI_Network();
	void copy(MPI_Network const & source);
	void clear();
	void get_forces(bool) override;
	// void optimize(float eta = 0.1, float alpha = 0.9, int max_iter = 800);
	void init_MPI(int world_rank, int world_size);
	//bool notmoving(int nodeid);
	void rewrite_moving_nodes();

	//int* not_moving_nodes;
	//int n_not_moving;
	int * chunk_nodes;
	int * chunk_edges;
	int chunk_edges_len;
	int chunk_nodes_len;
	int n_chunk_edges;





};

#endif
