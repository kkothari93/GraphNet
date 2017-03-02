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
//execute using mpirun ./a.out <filename>
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

//TODO: Define function 
inline bool Rinset(float* R, float x_lo, float x_hi, float y_lo, float y_hi){
	return (R[0] >= x_lo && R[0] <= x_hi && R[1] >= y_lo && R[1] <= y_hi);
}


int main(int argc, char* argv[]) {

	
	
	MPI_Init(NULL, NULL);

  	// Get the number of processes
  	int world_size;
  	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  	//cout << world_size << endl;
  	// Get the rank of the process
  	int world_rank;
  	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  	

  	Network * main_network = NULL;
  	string fname = "";
  	fname = argv[1];
  	//TODO: char[] vs string --> test small snippet
	main_network = new Network(fname);
	int chunk_size = ceil((main_network->n_nodes * 2)/world_size);
	
	float* plate_forces;
    
	main_network->get_weight();
	bool should_stop = main_network->get_stats();
	if(should_stop || world_size % 2 == 1) {
		//Always take even number of processors
		return 0;
	}
    
	if (world_rank == 0) {
		plate_forces = (float*)malloc(sizeof(float)*DIM*STEPS);
		memset(plate_forces, 0.0, STEPS*DIM*sizeof(*plate_forces));
	}
	
	int y_level = world_rank % 2 == 0? world_rank/2 : (world_rank-1)/2;
	int x_lo = world_rank%2 == 0? 0 : (MAXBOUND/2) + 1;
	int x_hi = world_rank%2 == 0? MAXBOUND/2 : MAXBOUND;
	int y_lo = ((MAXBOUND*2.0/world_size) * y_level);
	int y_hi = (world_rank == world_size - 1) || (world_rank == world_size - 2)? MAXBOUND : ((MAXBOUND*2.0/world_size) * (y_level+1)) - 1;


	//put them in a chunk specific array:
	// changed from /2 to *2
	main_network->chunk_nodes_len = main_network->n_nodes/world_size + sqrt(main_network->n_nodes)*2; 
	main_network->chunk_nodes = new int[main_network->chunk_nodes_len];
	int k = 0;
	// changed from i+=DIM to i+=1, use Rinset
	for (int i = 0; i < main_network->n_nodes; i+=1) {
		if (!Rinset(&(main_network->R[i*DIM]),x_lo, x_hi, y_lo, y_hi)) {
			main_network->chunk_nodes[k] = i;
			k++;
		}
	}
	// setting extra elements to -1
	for (; k < main_network->chunk_nodes_len; k++) {
		main_network->chunk_nodes[k] = -1;
	}

	//TODO: chunk node uniqueness check : sum each and check total against sum of 1+2+3+...+n_nodes

	//moved *DIM from chunk_edges_len calc to new int [] declaration, also changed /2 to *2
	main_network->chunk_edges_len = main_network->n_elems/world_size + sqrt(main_network->n_elems)*2;
	main_network->chunk_edges = new int[main_network->chunk_edges_len];
	k = 0;

	//changed i+=DIM to i+=1
	for (int i = 0; i < main_network->n_elems; i+=1) {

		// check this condition for both edges, use Rinset
		if (!Rinset(&(main_network->R[main_network->edges[i*2]*DIM]),x_lo, x_hi, y_lo, y_hi) && !Rinset(&(main_network->R[main_network->edges[i*2+1]*DIM]),x_lo, x_hi, y_lo, y_hi)) {
			continue;
		}
		else {
			// This is not correct logic: this would not add the counter point (b) for edge i (a,b)
			main_network->chunk_edges[k] = i;
			k++;
		}
	}
	for ( ; k < main_network->chunk_edges_len; k++) {
		main_network->chunk_edges[k] = -1;
	}
	

	//uniform L and PBC across all processors
	MPI_Bcast(main_network->L, main_network->n_elems, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(main_network->PBC, main_network->n_elems, MPI_C_BOOL, 0, MPI_COMM_WORLD); //not sure if it is randomly generated.
	//float * forces_buffer = new float[main_network->n_nodes * DIM * world_size];
	float * R_buffer = new float[main_network->n_nodes * DIM * world_size]; //buffer to gather the R from all nodes
	int * chunk_nodes_buffer = new int[main_network->chunk_nodes_len*world_size];
	MPI_Gather(main_network->chunk_nodes, main_network->chunk_nodes_len, MPI_INT, chunk_nodes_buffer, main_network->chunk_nodes_len, MPI_INT, 0, MPI_COMM_WORLD);
	//bool * broken_buffer = new bool[world_size];
	//int * edges_buffer = new int[main_network->n_elems * DIM * world_size];
	int iter = 0; // needed to write forces later
	clock_t t = clock(); 
	for(iter = 0; iter<STEPS; iter++){
		if((iter+1)%100 == 0){ // +1 required to have values in p_x, p_y
			cout<<(iter+1)<<endl; 
			cout<<"That took "<<(clock()-t)/CLOCKS_PER_SEC<<" s\n";
			t = clock();  // reset clock
			if(world_rank==0){
				main_network->get_stats();
			}
		}
		bool BROKEN = false;
		//TODO: add broken flag
		cout << __LINE__ << endl;
		main_network->optimize(BROKEN, x_lo, x_hi, y_lo, y_hi);
		cout << __LINE__ << endl;
		//TODO: use Network::get_plate_forces() on root proc
		//MPI_Gather(main_network->forces, main_network->n_nodes * DIM, MPI_FLOAT, forces_buffer, main_network->n_nodes * DIM, MPI_FLOAT, 0, MPI_COMM_WORLD);


		//MPI_Barrier(MPI_COMM_WORLD);


		//TODO: check the size of array parameter
		MPI_Gather(main_network->R, main_network->n_nodes * DIM, MPI_FLOAT, R_buffer, main_network->n_nodes * DIM, MPI_FLOAT, 0, MPI_COMM_WORLD);
		
		// syncing R
		if (world_rank == 0) {
			int node_to_sync  = 0;
			for (int i = 0; i < world_size; i += 1) {
				for (int j = i*main_network->chunk_nodes_len; j < (i+1)*main_network->chunk_nodes_len; j++) {
					node_to_sync = chunk_nodes_buffer[j];
					if (node_to_sync != -1) {
						main_network->R[DIM * node_to_sync] = R_buffer[main_network->n_nodes * DIM * i + DIM * node_to_sync];
						main_network->R[DIM * node_to_sync + 1] = R_buffer[main_network->n_nodes * DIM * i + DIM * node_to_sync + 1];
					}
				}
				// int i = (main_network->R[i]/chunk_size);
				// main_network->R[i] = R_buffer[main_network->n_nodes * DIM * i + i];
				//main_network->forces[i] = forces_buffer[main_network->n_nodes * DIM * i + i];
			}
			main_network->plotNetwork(iter, false);
			main_network->move_top_plate();

		}
		
		MPI_Bcast(main_network->R, main_network->n_nodes * DIM, MPI_FLOAT, 0, MPI_COMM_WORLD);

		cout << __LINE__ << endl;
		if (world_rank == 0) {
			main_network->get_plate_forces(plate_forces, iter);
			main_network->get_stats();
			
		}
		cout << __LINE__ << endl;
		MPI_Barrier(MPI_COMM_WORLD);
	}
	cout<<__LINE__<<endl;

	if (world_rank == 0) {
		string file_name = "forcesMPI.txt";
		write_to_file<float>(file_name, plate_forces, STEPS, DIM);
		free(plate_forces);
	}

	cout << __LINE__ << endl;
	MPI_Finalize();
}

// template <class T, class P>
// void mpi_gather_sync_bcast(P* input_array, int len, int world_rank, int world_size) {

// 	P* recv_buffer = new P[len*world_size];
// 	MPI_Gather(input_array, len, T, recv_buffer, len, T, 0, MPI_COMM_WORLD);

// }



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
