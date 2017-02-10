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
  	//TODO: char[] vs string --> test small snippet
	main_network = new Network(fname);


  	// int number_of_edges_per_process = Z_MAX*n_nodes*2/world_size;
  	// int * edges_this_process = new int[number_of_edges_per_process];

  	// int k = 0;
  	// for (int i = number_of_edges_per_process * world_rank; i < min(number_of_edges_per_process * (world_rank+1), Z_MAX*n_nodes*2); i++) {

  	// 	edges_this_process[k] = edges[i];
  	// 	k++;

  	// }
	//TODO: get_n_elems() -- n_nodes
	int chunk_size = ceil((main_network->get_n_elems()*2)/world_size);
	//if (chunk_size % 2 == 1) { chunk_size += 1;}

	int lo = world_rank * chunk_size;
	int hi = lo + chunk_size - 1;

	//uniform L across all processors
	MPI_Bcast(main_network->L, main_network->get_n_elems(), MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(main_network->PBC, main_network->get_n_elems(), MPI_BOOL, 0, MPI_COMM_WORLD); //not sure if it is randomly generated.
	float * recv_buffer;

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
		//TODO: add broken flag
		main_network->optimize(lo, hi);
		//TODO: use Network::get_plate_forces() on root proc
		get_components(p_x, p_y, pull_forces, iter);

		if (world_rank == 0) {
			move_top_nodes(R, tsideNodes, n_tside);
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

		float * recv_buffer = new float[main_network->get_n_elems()*2*world_size]
		if (world_rank == 0) {

		}
		//TODO: check the size of array parameter
		MPI_Gather(main_network->R, main_network->get_n_elems()*2, MPI_FLOAT, recv_buffer, main_network->get_n_elems()*2, MPI_FLOAT, 0, MPI_COMM_WORLD);

		// for (int i = 0; i < world_size; i++) {
		// 	for (int j = i * chunk_size; j < chunk_size * (i+1); j++) {
		// 		R[j] = recv_buffer[main_network->get_n_elems()*2*i + j];
		// 	}
		// }

		for (int i = 0; i < main_network->get_n_elems()*2; i++) {
			int proc_to_copy_from = (R[i]/chunk_size);
			R[i] = recv_buffer[main_network->get_n_elems()*2*proc_to_copy_from + i];
		}

		MPI_Bcast(main_network->R, main_network->get_n_elems()*2, MPI_FLOAT, 0, MPI_COMM_WORLD);
		//delete[] recv_buffer;
		MPI_Barrier(MPI_COMM_WORLD);
	}

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




















  	MPI_Finalize();

}