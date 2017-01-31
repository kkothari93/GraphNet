#ifndef NETWORK_H
#define NETWORK_H
#include <iostream>
#include <cmath>
#include "crack.h"
#include <cstdlib>
#include <random>
#include <time.h>
#include <vector>
#include <algorithm>
#include <ctime>
#include <chrono>
#include <fstream>
#include <sstream>
#include <string>

#define __params__	
#define TIME_STEP 1e-4						// Time step
#define SIM_TIME 10.0							// Std. deviation for contour lengths
#define Z_MAX 10							// Max coordination number
#define TOL 1e-6
#define MAXBOUND 500.0
#define L_MEAN 120.0						// Average for contour length
#define L_STD 4.0	

// Define constants
#define __constants__

#define kB 1.38064852e-5					// Boltzmann constant
#define b 0.1								// Persistence length
#define T 300 								// Temperature
#define ae 0.1 								// Strength of bond - includes activation energy
#define delxe 0.15 							// parameter for breaking crosslink connection
#define BLOCK_SIZE 1024

//using namespace std;

class Network {

public:

	Network();
	Network(Network const & source);
	virtual ~Network();
	Network const & operator=(Network const & other);
	void build_network();
	void apply_crack(Crack const & crack);
	void load_network();
	void make_edge_connections(double dely_allowed = 10.0);
	void get_forces(const float* PBC_vector, bool update_damage = false);


	

private:

	void clear();
	void copy(Network const & source);
	
	bool cracked;
	int DIM;
	int n_nodes;
	int n_elems;
	int n_rside, n_lside, n_bside, n_tside;
	//int num_edges;
	double * R;
	int * edges;
	double * forces;
	double * damage;
	bool ** edge_matrix;
	double * L;
	bool* PBC;
	int* lsideNodes;
	int* rsideNodes;
	int* tsideNodes;
	int* bsideNodes;

	
};

#include "Network.cpp"
#endif