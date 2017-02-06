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
#include "Network.h"


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

	
	MPI_Init(NULL, NULL);

  	// Get the number of processes
  	int world_size;
  	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  	// Get the rank of the process
  	int world_rank;
  	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  	// Get the name of the processor
  	// char processor_name[MPI_MAX_PROCESSOR_NAME];
  	// int name_len;
  	// MPI_Get_processor_name(processor_name, &name_len);

  	Network * main_network = NULL;
  	String fname = "";
	main_network = new Network(fname);


  	// int number_of_edges_per_process = Z_MAX*n_nodes*2/world_size;
  	// int * edges_this_process = new int[number_of_edges_per_process];

  	// int k = 0;
  	// for (int i = number_of_edges_per_process * world_rank; i < min(number_of_edges_per_process * (world_rank+1), Z_MAX*n_nodes*2); i++) {

  	// 	edges_this_process[k] = edges[i];
  	// 	k++;

  	// }

	float * local_R;
	int * local_edges;

	main_network->split_for_MPI(local_R, local_edges, NULL, world_size, world_rank);

	int iter = 0; // needed to write forces later

	for(iter = 0; iter<STEPS; iter++){
		if((iter+1)%1000 == 0){ // +1 required to have values in p_x, p_y
			cout<<(iter+1)<<endl; 
			cout<<"That took "<<(clock()-t)/CLOCKS_PER_SEC<<" s\n";
			t = clock();  // reset clock

			// plot network
			//plot_network(gnetwork, R, edges, PBC, \
			 	//n_nodes, n_elems, iter);
			
			// Plot forces
			//plot_forces(gforces, p_x, p_y, iter);


		}
		optimize(local_R, local_edges, damage, L, n_nodes, n_elems,\
			PBC, PBC_vector, tsideNodes, n_tside, bsideNodes, n_bside,\
			pull_forces, iter);
	
		get_components(p_x, p_y, pull_forces, iter);

		if (world_rank == 0) {
			move_top_nodes(R, tsideNodes, n_tside);
		}
		
		//sync step:
		//receive one above your topmost from prev_proc
		//send your topmost to the prev proc
		//send your bottommost to next proc
		//receive one below your bottommost from next proc

	}
















  	MPI_Finalize();

}