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
// #include "vel.h"
// #include "gnuplot_i.hpp"
#include "network.h"
#include "sac_network.h"
using namespace std;

class MPI_Network : public Network {

public:

	MPI_Network();
	MPI_Network(MPI_Network const & source);
	MPI_Network(Network const & source);
	MPI_Network(sacNetwork const & source);
	~MPI_Network();
	void copy(MPI_Network const & source);
	void clear();
	//void malloc_network();
	void load_network(string& fname);
	void get_forces(bool);
	void optimize(float eta = 0.1, float alpha = 0.9, int max_iter = 800);
	void init_MPI(int world_rank, int world_size);
	bool notmoving(int nodeid);


	int* not_moving_nodes;
	int n_not_moving;
	int * chunk_nodes;
	int * chunk_edges;
	int chunk_edges_len;
	int chunk_nodes_len;





};

#endif
