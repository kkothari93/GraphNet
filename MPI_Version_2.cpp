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


//#include "input.h"
#include "Network.h"
#define PAD MAXBOUND*1.03
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
	return (R[0] >= x_lo && R[0] < x_hi && R[1] >= y_lo && R[1] < y_hi);
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
  	string fname = "template2d.msh";
  	//fname = argv[1];
  	//TODO: char[] vs string --> test small snippet
	main_network = new Network(fname);
	int chunk_size = ceil((main_network->n_nodes * 2)/world_size);
	
	float* plate_forces = NULL;
    
	// main_network->get_weight();
	// bool should_stop = main_network->get_stats();
	// if(should_stop || world_size % 2 == 1) {
	// 	//Always take even number of processors
	// 	return 0;
	// }
    
	if (world_rank == 0) {
		plate_forces = (float*)malloc(sizeof(float)*DIM*STEPS);
		memset(plate_forces, 0.0, STEPS*DIM*sizeof(*plate_forces));
	}
	
	int y_level = world_rank % 2 == 0? world_rank/2 : (world_rank-1)/2;
	float x_lo = world_rank%2 == 0? 0 : (PAD/2);
	float x_hi = world_rank%2 == 0? PAD/2 : PAD;
	float y_lo = ((PAD*2.0/world_size) * y_level);
	float y_hi = (world_rank == world_size - 1) || (world_rank == world_size - 2)? PAD : ((PAD*2.0/world_size) * (y_level+1));


	//put them in a chunk specific array:
	// changed from /2 to *2
	main_network->chunk_nodes_len = main_network->n_nodes/world_size + sqrt(main_network->n_nodes)*2; 
	main_network->chunk_nodes = new int[main_network->chunk_nodes_len];
	int k = 0, n_chunk_nodes;
	// changed from i+=DIM to i+=1, use Rinset
	for (int i = 0; i < main_network->n_nodes; i+=1) {
		if (Rinset(&(main_network->R[i*DIM]),x_lo, x_hi, y_lo, y_hi)) {
			main_network->chunk_nodes[k] = i;
			k++;
		}
	}
	n_chunk_nodes = k;
	// setting extra elements to -1
	for (; k < main_network->chunk_nodes_len; k++) {
		main_network->chunk_nodes[k] = -1;
	}

	



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
	printf("World rank: %d, x_lo, x_hi, y_lo, y_hi %f, %f, %f, %f\n", world_rank, x_lo, x_hi, y_lo, y_hi);
	printf("World rank: %d, number of edges %d, number of nodes in chunk %d\n", world_rank, k, n_chunk_nodes);
	for ( ; k < main_network->chunk_edges_len; k++) {
		main_network->chunk_edges[k] = -1;
	}
	//cout << __LINE__ << endl;
	
	//cout << __LINE__ << endl;
	//uniform L and PBC across all processors
	MPI_Bcast(main_network->L, main_network->n_elems, MPI_FLOAT, 0, MPI_COMM_WORLD);
	cout << world_rank << "  "<<__LINE__ << endl;
	MPI_Bcast(main_network->PBC, main_network->n_elems, MPI_C_BOOL, 0, MPI_COMM_WORLD); //not sure if it is randomly generated.
	// cout << world_rank << endl;
	// cout << world_rank << "  "<<__LINE__ << endl;
	size_t r_size = main_network->n_nodes * DIM * world_size;
	

	float * R_buffer; 
	R_buffer = (float*)malloc(main_network->n_nodes * DIM * world_size*sizeof(float));//main_network->n_nodes * DIM * world_size]; //buffer to gather the R from all nodes
	cout << world_rank << "  "<<__LINE__ << endl;
	int * chunk_nodes_buffer = new int[main_network->chunk_nodes_len*world_size];
	cout << world_rank << "  "<<__LINE__ << endl;
	//MPI_Gather(main_network->chunk_nodes, main_network->chunk_nodes_len, MPI_INT, chunk_nodes_buffer, main_network->chunk_nodes_len, MPI_INT, 0, MPI_COMM_WORLD);
	// cout << world_rank << "  "<<__LINE__ << endl;
	//TODO: chunk node uniqueness check : sum each and check total against sum of 1+2+3+...+n_nodes
	// if (world_rank == 0) {
	// 	int chunk_sum = 0;
	// 	for (int i = 0; i < main_network->chunk_nodes_len * world_size; i++) {
	// 		if (chunk_nodes_buffer[i] != -1) {
	// 			chunk_sum += chunk_nodes_buffer[i];
	// 		}
	// 	}
	// 	if (chunk_sum != (main_network->n_nodes)*(main_network->n_nodes + 1)/2) {
	// 		cout << "Uneven chunk partitioning" << endl;
	// 	      	//return 0;
	// 	}
	// }
	
	// cout << "Starting the loop" << endl;
	// //bool * broken_buffer = new bool[world_size];
	// //int * edges_buffer = new int[main_network->n_elems * DIM * world_size];
	// int iter = 0; // needed to write forces later
	// clock_t t = clock(); 
	// for(iter = 0; iter<STEPS; iter++){
	// 	if((iter+1)%100 == 0){ // +1 required to have values in p_x, p_y
	// 		cout<<(iter+1)<<endl; 
	// 		cout<<"That took "<<(clock()-t)/CLOCKS_PER_SEC<<" s\n";
	// 		t = clock();  // reset clock
	// 		if(world_rank==0){
	// 			main_network->get_stats();
	// 		}
	// 	}
	// 	bool BROKEN = false;
	// 	//TODO: add broken flag
	// 	cout << __LINE__ << endl;
	// 	main_network->optimize(BROKEN, x_lo, x_hi, y_lo, y_hi);
	// 	cout << __LINE__ << endl;
	// 	//TODO: use Network::get_plate_forces() on root proc
	// 	//MPI_Gather(main_network->f:orces, main_network->n_nodes * DIM, MPI_FLOAT, forces_buffer, main_network->n_nodes * DIM, MPI_FLOAT, 0, MPI_COMM_WORLD);


	// 	//MPI_Barrier(MPI_COMM_WORLD);


	// 	//TODO: check the size of array parameter
	// 	MPI_Gather(main_network->R, main_network->n_nodes * DIM, MPI_FLOAT, R_buffer, main_network->n_nodes * DIM, MPI_FLOAT, 0, MPI_COMM_WORLD);
		
	// 	// syncing R
	// 	if (world_rank == 0) {
	// 		int node_to_sync  = 0;
	// 		for (int i = 0; i < world_size; i += 1) {
	// 			for (int j = i*main_network->chunk_nodes_len; j < (i+1)*main_network->chunk_nodes_len; j++) {
	// 				node_to_sync = chunk_nodes_buffer[j];
	// 				if (node_to_sync != -1) {
	// 					main_network->R[DIM * node_to_sync] = R_buffer[main_network->n_nodes * DIM * i + DIM * node_to_sync];
	// 					main_network->R[DIM * node_to_sync + 1] = R_buffer[main_network->n_nodes * DIM * i + DIM * node_to_sync + 1];
	// 				}
	// 			}
	// 			// int i = (main_network->R[i]/chunk_size);
	// 			// main_network->R[i] = R_buffer[main_network->n_nodes * DIM * i + i];
	// 			//main_network->forces[i] = forces_buffer[main_network->n_nodes * DIM * i + i];
	// 		}
	// 		main_network->plotNetwork(iter, false);
	// 		main_network->move_top_plate();

	// 	}
		
	// 	MPI_Bcast(main_network->R, main_network->n_nodes * DIM, MPI_FLOAT, 0, MPI_COMM_WORLD);

	// 	cout << __LINE__ << endl;
	// 	if (world_rank == 0) {
	// 		main_network->get_plate_forces(plate_forces, iter);
	// 		main_network->get_stats();
			
	// 	}
	// 	cout << __LINE__ << endl;
	// 	MPI_Barrier(MPI_COMM_WORLD);
	// }
	// cout<<__LINE__<<endl;

	// if (world_rank == 0) {
	// 	string file_name = "forcesMPI.txt";
	// 	write_to_file<float>(file_name, plate_forces, STEPS, DIM);
	// 	free(plate_forces);
	// }

	// cout << __LINE__ << endl;
	//MPI_Finalize();
	// cout << __LINE__ << endl;
	free(R_buffer);
	return 1;
}

/**
int main(int argc, char* argv[]) {
	MPI_Init(NULL, NULL);
	// Get the number of processes
  	int world_size;
  	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  	cout << world_size << endl;
  	// Get the rank of the process
  	int world_rank;
  	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  	string fname = "template2d.msh";
  	Network main_network(fname);
  	//delete main_network;
	MPI_Finalize();
}
**/