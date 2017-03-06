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
//#include "gnuplot_i.hpp"
using namespace std;

#define __params__	
#define DIM 2								// Number of dimensions
#define TIME_STEP 1e-4						// Time step
#define SIM_TIME 5.0						// Simulation time
#define TOL 1e-6							// Tolerance
#define STEPS 10 //int(SIM_TIME/TIME_STEP)		// Number of time steps
#define L_MEAN 120.0f						// Average for contour length
#define L_STD 4.0f							// Std. deviation for contour lengths
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

const float PBC_vector[DIM] = {MAXBOUND*1.2, 0};
const float vel[DIM] = {vel_x, vel_y};

template <typename t>
void write_to_file(string& fname, t* arr, int rows, int cols){

	ofstream logger;
	std::time_t result = std::time(nullptr);


	cout<<"1D pulling of a 2D gel"<<endl;
	cout<<"File created at "<<std::asctime(std::localtime(&result));
	cout<<endl;

	cout<<"Sim dimension : "<<DIM<<endl;
	cout<<"Simulation time : "<<SIM_TIME<<endl;
	cout<<"Velocity : "<<vel_x<<"\t"<<vel_y<<endl;
	cout<<endl;

	cout<<"Disorder characteristics : "<<endl;
	cout<<" -- L_MEAN : "<<L_MEAN<<endl;
	cout<<" -- L_STD : "<<L_STD<<endl;
	cout<<endl;
	
	cout<<"Cracked? : "<<CRACKED<<endl;

	logger.open(fname, ios::trunc|ios_base::out);
	for(int i =0; i < rows; i++){
		for(int j = 0; j< cols; j++){
			logger<<arr[i*cols + j]<<"\t";
		}
		logger<<"\n";
	}
}

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
	//bool get_forces(bool, int lo = 0, int hi = 0);
	bool get_forces(bool, int x_lo = 0, int x_hi = 0, int y_lo = 0, int y_hi = 0);
	void move_top_plate();
	void get_plate_forces(float*, int);
	//void optimize(bool& BROKEN, int lo, int hi, float eta = 0.1, float alpha = 0.9, int max_iter = 800);
	void optimize(bool& BROKEN, int x_lo, int x_hi, int y_lo, int y_hi, float eta = 0.1, float alpha = 0.9, int max_iter = 800);
	void split_for_MPI(float * R_split, int * edges_split, float * forces, int number_of_procs, int curr_proc_rank);
	int get_current_edges();
	bool get_stats();
	float get_weight();
	void set_weight(float weight);
	void plotNetwork(int iter_step, bool first_time);

//private:

	void clear();
	void copy(Network const & source);
	//Gnuplot gnu;
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
	int * chunk_nodes;
	int * chunk_edges;
	int chunk_edges_len;
	int chunk_nodes_len;

	
};

//#include "Network.cpp"
#endif