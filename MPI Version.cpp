#include <mpi.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <random>
#include <fstream>
#include <time.h>
#include <algorithm>
#include <ctime>
#include <chrono>


#include "gnuplot_i.hpp"
#include "input.h"


#define DIM 2								// Number of dimensions
#define TIME_STEP 1e-4						// Time step
#define SIM_TIME 10.0						// Simulation time
#define TOL 1e-6							// Tolerance
#define STEPS int(SIM_TIME/TIME_STEP)		// Number of time steps
#define L_MEAN 120.0f						// Average for contour length
#define L_STD 4.0f							// Std. deviation for contour lengths
#define Z_MAX 10							// Max coordination number
#define MAXBOUND 500.0f
#define CRACKED true

// Define constants
#define kB 1.38064852e-5					// Boltzmann constant
#define b 0.1								// Persistence length
#define T 300 								// Temperature
#define ae 0.1 								// Strength of bond - includes activation energy
#define delxe 0.15 							// parameter for breaking crosslink connection

static double vel[DIM] = {0.0, 100.0};


int main() {

	int n_nodes=1600, n_elems = 11000;
	double * R =  new double[n_nodes*DIM];
	int * edges = new int[Z_MAX*n_nodes*2];
	double * pull_forces = new double[STEPS*DIM];

	vector<double> p_x;
	vector<double> p_y;

	// file to read and write position data
	const string fname = "coordinates.txt";

	// Read in mesh
	cout<<"Reading the mesh...\n";	
	take_input(R, edges, n_nodes, n_elems);
	cout<<"Mesh read successfully!\n";
	cout<<"Number of nodes are: "<<n_nodes<<endl;
	cout<<"Number of elements are: "<<n_elems<<endl;

	int max_nodes_on_a_side = int(sqrt(n_nodes))*2;
	int n_rside = 0, n_lside = 0, n_bside = 0, n_tside = 0;
	// get left and right nodes
	int * lsideNodes = new int[max_nodes_on_a_side];
	int * rsideNodes = new int[max_nodes_on_a_side];
	// get top and bottom nodes
	int * tsideNodes = new int[max_nodes_on_a_side];
	int * bsideNodes = new int[max_nodes_on_a_side];

	side_nodes(R, lsideNodes, rsideNodes, tsideNodes, bsideNodes,\
		n_lside, n_rside, n_tside, n_bside, n_nodes);	

	// Initialize edge properties
	double * damage = new doubel[2*n_elems];
	double * L = new doubel[2*n_elems];
	bool * PBC = new bool[2*n_elems];
	__init__(L, damage, PBC, n_elems);

	// Make PBC connections
	const double PBC_vector[DIM] = {MAXBOUND*1.2, 0};
	// TODO: get lside and rside nodes
	make_edge_connections(R, edges, n_elems, \
		lsideNodes, rsideNodes, n_lside, n_rside, \
		PBC, L, damage, 15.0);
	cout<<"Number after new connections made: "<<n_elems<<endl;

	if(CRACKED){
		double c[] = {MAXBOUND/2.0, MAXBOUND/2.0};
		double a[] = {MAXBOUND/20.0, MAXBOUND/25.0};
		crack(c, a, R, n_nodes, edges, n_elems);
	}
	
	// gnuplots
	//Gnuplot gforces("lines lw 2");
	//Gnuplot gnetwork;

	//plot_network(gnetwork, R, edges, PBC, \
	//			n_nodes, n_elems, 0);
	// Track time for 1000 iterations
	clock_t t = clock(); 
	
	int iter = 0; // needed to write forces later

	MPI_Init(NULL, NULL);

  	// Get the number of processes
  	int world_size;
  	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  	// Get the rank of the process
  	int world_rank;
  	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  	// Get the name of the processor
  	char processor_name[MPI_MAX_PROCESSOR_NAME];
  	int name_len;
  	MPI_Get_processor_name(processor_name, &name_len);






















}