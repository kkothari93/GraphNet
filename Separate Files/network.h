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
#include "gnuplot_i.hpp"
#include "helper_funcs.cpp"
using namespace std;


#define __params__	
#define DIM 2								// Number of dimensions
#define TIME_STEP 1e-4						// Time step
#define SIM_TIME 9.0						// Simulation time
#define TOL 1e-6							// Tolerance
#define STEPS int(SIM_TIME/TIME_STEP)		// Number of time steps
#define L_MEAN 250.0f						// Average for contour length
#define L_STD 150.0f							// Std. deviation for contour lengths
#define MAXBOUND 500.0f
#define SACBONDS false
#define IMPLEMENT_PBC true
#define FNAME_STRING "high_disorder_high_L_"
#define CRACKED false
//#define GENERATOR std::uniform_real_distribution<float>(L_MEAN - L_STD, L_MEAN + L_STD
#if CRACKED
#define PROB_REMOVAL 0.8
#else
#define PROB_REMOVAL 0.0
#endif


// Define constants
#define __constants__

#define kB 1.38064852e-5					// Boltzmann constant
#define b 0.1								// Persistence length
#define T 300 								// Temperature
#define ae 0.1 								// Strength of bond - includes activation energy
#define delxe 0.15 							// parameter for breaking crosslink connection
#define af 0.3
#define delxf 0.25

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
	virtual void build_network();
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

//protected:

	virtual void clear();
	void copy(Network const & source);
	Gnuplot gnu;
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

#endif