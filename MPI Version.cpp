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


//#include "gnuplot_i.hpp"
#include "input.h"
#include "Network.h"
//compile using: mpic++ MPI\ Version.cpp Network.h Network.cpp crack.h Crack.cpp
//execute using mpirun ./a.out
//makefile not working


// #define DIM 2								// Number of dimensions
// #define TIME_STEP 1e-4						// Time step
// #define SIM_TIME 10.0						// Simulation time
// #define TOL 1e-6							// Tolerance
// #define STEPS int(SIM_TIME/TIME_STEP)		// Number of time steps
// #define L_MEAN 120.0f						// Average for contour length
// #define L_STD 4.0f							// Std. deviation for contour lengths
// #define Z_MAX 10							// Max coordination number
// #define MAXBOUND 500.0f
// #define CRACKED true

// // Define constants
// #define kB 1.38064852e-5					// Boltzmann constant
// #define b 0.1								// Persistence length
// #define T 300 								// Temperature
// #define ae 0.1 								// Strength of bond - includes activation energy
// #define delxe 0.15 							// parameter for breaking crosslink connection

//static double vel[DIM] = {0.0, 100.0};


int main(int argc, char* argv[]) {

	
	MPI_Init(NULL, NULL);

  	// Get the number of processes
  	int world_size;
  	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  	// Get the rank of the process
  	int world_rank;
  	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  	

  	Network * main_network = NULL;
  	string fname = "";
  	fname = argv[1];
  	//TODO: char[] vs string --> test small snippet
	main_network = new Network(fname);
	int chunk_size = ceil((main_network->n_nodes * 2)/world_size);
	
	int lo = world_rank * chunk_size;
	int hi = lo + chunk_size - 1;

	//uniform L and PBC across all processors
	MPI_Bcast(main_network->L, main_network->n_elems, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(main_network->PBC, main_network->n_elems, MPI_C_BOOL, 0, MPI_COMM_WORLD); //not sure if it is randomly generated.
	float * recv_buffer = new float[main_network->n_nodes * DIM * world_size]; //buffer to gather the R from all nodes
	bool * broken_buffer = new bool[world_size];
	int * edges_buffer = new int[main_network->n_elems * DIM * world_size];
	int iter = 0; // needed to write forces later
	clock_t t = clock(); 
	for(iter = 0; iter<STEPS; iter++){
		if((iter+1)%1000 == 0){ // +1 required to have values in p_x, p_y
			cout<<(iter+1)<<endl; 
			cout<<"That took "<<(clock()-t)/CLOCKS_PER_SEC<<" s\n";
			t = clock();  // reset clock
		}
		bool BROKEN = false;
		//TODO: add broken flag
		main_network->optimize(BROKEN, lo, hi);
		//TODO: use Network::get_plate_forces() on root proc

		if (world_rank == 0) {
			//move_top_nodes(R, tsideNodes, n_tside);
			main_network->move_top_plate();
			main_network->get_plate_forces(main_network->forces/** TODO: not sure what to put here **/, iter);
		}
		
		//TODO: check the size of array parameter
		MPI_Gather(main_network->R, main_network->n_nodes * DIM, MPI_FLOAT, recv_buffer, main_network->n_nodes * DIM, MPI_FLOAT, 0, MPI_COMM_WORLD);
		if (world_rank == 0) {
			for (int i = 0; i < main_network->n_nodes * DIM; i++) {
				int proc_to_copy_from = (main_network->R[i]/chunk_size);
				main_network->R[i] = recv_buffer[main_network->n_nodes * DIM * proc_to_copy_from + i];
			}
		}
		
		MPI_Bcast(main_network->R, main_network->n_nodes * DIM, MPI_FLOAT, 0, MPI_COMM_WORLD);


		//checking if edges need to be resynced
		MPI_Allgather(&BROKEN, 1, MPI_C_BOOL, broken_buffer, 1, MPI_C_BOOL, MPI_COMM_WORLD);
		
		bool need_edges_resyncing = false;
		for (int i = 0; i < world_size; i++) {
			need_edges_resyncing = need_edges_resyncing || broken_buffer[i];
			if (need_edges_resyncing) {
				break;
			}
		}
		if (need_edges_resyncing) {
			MPI_Gather(main_network->edges, main_network->n_elems * DIM, MPI_INT, edges_buffer, main_network->n_elems * DIM, MPI_INT, 0, MPI_COMM_WORLD);
			if (world_rank == 0) {
				for (int i = 0; i < main_network->n_elems * DIM; i += DIM) {
					int proc_to_copy_from = i / chunk_size;
					main_network->edges[i] = edges_buffer[main_network->n_elems * DIM * proc_to_copy_from + i];
					main_network->edges[i+1] = edges_buffer[main_network->n_elems * DIM * proc_to_copy_from + i + 1];
				}
			}
			MPI_Bcast(main_network->edges, main_network->n_elems * DIM, MPI_INT, 0, MPI_COMM_WORLD);
			
		}

		MPI_Barrier(MPI_COMM_WORLD);
	}

	MPI_Finalize();

}



//sync step:
		//receive one above your topmost from prev_proc
		//send your topmost to the prev proc
		//send your bottommost to next proc
		//receive one below your bottommost from next proc

		// float * topmost_R_recv_buffer = new float[R_buffer];
		// if (world_rank != 0) {
		// 	MPI_Recv(topmost_R_recv_buffer, R_buffer, MPI_FLOAT, world_rank-1, 0, MPI_COMM_WORLD, NULL);
		// }
		
		// float * topmost_R_send_buffer = new float[R_buffer];
		// if (world_rank != 0) {
		// 	for (int i = 0; i < R_buffer; i++) {
		// 		topmost_R_send_buffer[i] = local_R[R_buffer+i];
		// 	}
		// 	MPI_Send(topmost_R_send_buffer, R_buffer, MPI_FLOAT, world_rank-1, 1, MPI_COMM_WORLD);
		// }
		
		// float * bottommost_R_send_buffer = new float[R_buffer]
		// for (int i = 0; i < R_buffer; i++) {
		// 	bottommost_R_send_buffer[i] = local_R[/** to fill **/];
		// }
		// MPI_Send(bottommost_R_send_buffer, R_buffer, MPI_FLOAT, world_rank+1, 2, MPI_COMM_WORLD);
		// float * bottommost_R_recv_buffer = new float[R_buffer];
		// MPI_Recv(bottommost_R_recv_buffer, R_buffer, MPI_FLOAT, world_rank+1, 3, MPI_COMM_WORLD, NULL);

		// for (int i = 0; i < R_buffer; i++) {
		// 	local_R[i] = topmost_R_recv_buffer[i];
		// 	local_R[/** to fill**/ + i] = bottommost_R_recv_buffer[i];
		// }

		// MPI_Barrier(MPI_COMM_WORLD);

	//aggregation:



	// float *  R_recv_buffer;
	// float *  R_send_buffer = new float[R_split_size];
	// if (world_rank == 0) {
	// 	R_recv_buffer = new float[/**to fill**/];
	// }
	// for (int i = 0; i < R_split_size; i++) {
	// 	R_send_buffer[i] = local_R[R_buffer+i];
	// }

	// MPI_Gather(R_send_buffer, R_buffer, MPI_FLOAT, R_recv_buffer, R_buffer, MPI_FLOAT, 0, MPI_COMM_WORLD);

	// MPI_Barrier(MPI_COMM_WORLD);

// plot network
			//plot_network(gnetwork, R, edges, PBC, \
			 	//n_nodes, n_elems, iter);
			
			// Plot forces
			//plot_forces(gforces, p_x, p_y, iter);