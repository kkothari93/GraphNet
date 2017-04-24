#ifndef NETWORK_H
#define NETWORK_H
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
#include "vel.h"

/**
@file network.h
\brief Conceptualizes and builds a polymer network as a graph
G = (V,E) (i.e. a set of nodes (crosslinkers, V) and a set of
edges (polymers, E)).
*/
#include "helper_funcs.h"
#include "crack.h"
using namespace std;

const float PBC_vector[DIM] = {MAXBOUND_X*1.1, 0};
const float vel[DIM] = {vel_x, vel_y};

class Network {
/*!<
\brief Implements a network of discrete chains.

This is the base class that takes in a GMSH .msh file to build 
a network. All edges in the mesh generated from GMSH are taken 
as edges of the network. The nodes in the mesh are taken as 
cross-linker positions and the edges as polymers. 

All properties relevant to the computational experiments are initialized.
Experiments are then member functions that can be called.
*/
public:
	Network();
	Network(Network const & source);
	Network(string& fname);
	virtual ~Network();
	Network const & operator=(Network const & other);
	virtual void build_network();
	void side_nodes();
	void remove_duplicates(int&);
	void apply_crack(Cracklist &);
	virtual void load_network(string&);
	virtual void malloc_network(string&);
	void make_edge_connections(float dely_allowed = 10.0);
	virtual void get_forces(bool);
	virtual void move_top_plate();
	virtual void get_plate_forces(float*, int);
	virtual void optimize(float eta = 0.1, float alpha = 0.9, int max_iter = 800);
	float get_weight();
	float set_weight(float);
	bool get_stats();
	int get_current_edges();
	virtual void plotNetwork(int, bool);
	virtual void clear();
	void copy(Network const & source);
	void dump(int,bool first_time = false);

	bool cracked; 		///<Flag to check if the network has cracks
	int n_nodes; 		///<Stores number of crosslinker nodes in the network
	int n_elems; 		///<Stores number of polymer connections in the network

	int n_rside; 		///<Stores number of crosslinker nodes on the right edge of the network sample
	int n_lside; 		///<Stores number of crosslinker nodes on the left edge of the network sample
	int n_bside; 		///<Stores number of crosslinker nodes on the bottom edge of the network sample
	int n_tside; 		///<Stores number of crosslinker nodes on the top edge of the network sample

	float * R; 			///<Stores the position of all nodes of the graph in a flat n_nodes*DIM size array
	int * edges; 		///<Stores the i,j pairs that form edges, size is n_elems*2

	float * forces; 	///<Stores the forces on all nodes of the network in a flat n_nodes*DIM size array
	float * damage; 	///<Stores the damage in each edge of the network
	float * L; 			///<Stores the contour length of the edge
	bool* PBC; 			///<Tag to check if an edge is a periodic boundary condition edge

	int* lsideNodes; 	///<Stores index of crosslinker nodes on the left edge of the network sample
	int* rsideNodes; 	///<Stores index of crosslinker nodes on the right edge of the network sample
	int* tsideNodes; 	///<Stores index of crosslinker nodes on the top edge of the network sample
	int* bsideNodes; 	///<Stores index of crosslinker nodes on the bottom edge of the network sample
	
	bool initialized;	///<Internal variable to check if Network object is initialized
	
	//add moving nodes to speed up force
	int* moving_nodes;	///<Stores index of all nodes where no boundary condition is applied (i.e. they participate in optimization)
	int n_moving; 		///<Stores number of moving nodes
};

#endif