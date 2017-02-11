#ifndef NETWORK_H
#define NETWORK_H
#include <iostream>
#include <cmath>
#include "crack.h"
#include <cstdlib>
#include <math.h>
#include <random>
#include <time.h>
#include <vector>
#include <algorithm>
#include <ctime>
#include <chrono>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include "crack.h"
#include <stddef.h>
#include "vel.h"
using namespace std;

#define __params__	
#define DIM 2								// Number of dimensions
#define TIME_STEP 1e-4						// Time step
#define SIM_TIME 5.0						// Simulation time
#define TOL 1e-6							// Tolerance
#define STEPS 10// int(SIM_TIME/TIME_STEP)		// Number of time steps
#define L_MEAN 120.0f						// Average for contour length
#define L_STD 6.0f							// Std. deviation for contour lengths
#define Z_MAX 10							// Max coordination number
#define MAXBOUND 500.0f
#define CRACKED false

// Define constants
#define __constants__

#define kB 1.38064852e-5					// Boltzmann constant
#define b 0.1								// Persistence length
#define T 300 								// Temperature
#define ae 0.1 								// Strength of bond - includes activation energy
#define delxe 0.15 							// parameter for breaking crosslink connection
#define BLOCK_SIZE 1024

using namespace std;

const float PBC_vector[DIM] = {MAXBOUND*1.1, 0};
const float vel[DIM] = {vel_x, vel_y};

class Network {

public:
	Network();
	Network(Network const & source);
	Network(string& fname);
	virtual ~Network();
	Network const & operator=(Network const & other);
	void build_network();
	void apply_crack(Crack const & crack);
	void load_network(string&);
	void malloc_network(string&);
	void make_edge_connections(float dely_allowed = 10.0);
	void get_forces(bool);
	void move_top_plate();
	void get_plate_forces(float*, int);
	void optimize(float, float, int);
	float get_weight();
	void set_weight(float);
	bool get_stats();
	int get_current_edges();
	void split_for_MPI(float * R_split, int * edges_split, float * forces, int number_of_procs, int curr_proc_rank);

private:

	void clear();
	void copy(Network const & source);
	
	bool cracked;
	//int DIM;
	int n_nodes;
	int n_elems;
	int n_rside, n_lside, n_bside, n_tside;
	//int num_edges;
	float * R; //n_nodes*DIM
	int * edges; //n_elems*2
	float * forces; //n_nodes*DIM
	float * damage;// n_nodes
	bool ** edge_matrix; 
	float * L; //n_elems
	bool* PBC; //n_elems
	int* lsideNodes; //max_nodes_on_a_side*2
	int* rsideNodes;
	int* tsideNodes;
	int* bsideNodes;
	bool initialized;
	//add moving nodes to speed up force
	int* moving_nodes;
	int n_moving;

	
};

//#include "Network.cpp"
#endif