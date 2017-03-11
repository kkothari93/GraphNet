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

#include "Network.h"
#define PAD MAXBOUND*1.03
#define NSYNC 1
//compile using: mpic++ MPI\ Version.cpp Network.h Network.cpp crack.h Crack.cpp
//execute using mpirun ./a.out <filename>
//makefile not working

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
    
	main_network->get_weight();
	bool should_stop = main_network->get_stats();
	if(should_stop || world_size % 2 == 1) {
		//Always take even number of processors
		return 0;
	}
    
	if (world_rank == 0) {
		plate_forces = (float*)malloc(sizeof(float)*DIM*STEPS);
		memset(plate_forces, 0, STEPS*DIM*sizeof(float));
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
	main_network->chunk_edges_len = main_network->n_elems/world_size + sqrt(main_network->n_elems)*5;
	main_network->chunk_edges = new int[main_network->chunk_edges_len];
	k = 0;
	int n_chunk_edges;
	//changed i+=DIM to i+=1
	for (int i = 0; i < main_network->n_elems; i+=1) {

		// check this condition for both edges, use Rinset
		if (!Rinset(&(main_network->R[main_network->edges[i*2]*DIM]),x_lo, x_hi, y_lo, y_hi) && !Rinset(&(main_network->R[main_network->edges[i*2+1]*DIM]),x_lo, x_hi, y_lo, y_hi)) {
			continue;
		}
		else {
			main_network->chunk_edges[k] = i;
			k++;
		}
	}
	n_chunk_edges = k;
	printf("World rank: %d, x_lo, x_hi, y_lo, y_hi %f, %f, %f, %f\n", world_rank, x_lo, x_hi, y_lo, y_hi);
	printf("World rank: %d, number of edges %d, number of nodes in chunk %d\n", world_rank, n_chunk_edges, n_chunk_nodes);
	for ( ; k < main_network->chunk_edges_len; k++) {
		main_network->chunk_edges[k] = -1;
	}

	//uniform L and PBC across all processors
	MPI_Bcast(main_network->L, main_network->n_elems, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(main_network->PBC, main_network->n_elems, MPI_C_BOOL, 0, MPI_COMM_WORLD); 
	
	size_t r_size = main_network->n_nodes * DIM * world_size;
	float * R_buffer; 
	R_buffer = (float*)malloc(r_size*sizeof(float));//buffer to gather the R from all nodes
	float * forces_buffer; 
	forces_buffer = (float*)malloc(r_size*sizeof(float));//buffer to gather the R from all nodes


	int * chunk_nodes_buffer = new int[main_network->chunk_nodes_len*world_size];
	// TODO: Sync forces every nth iteration

	MPI_Gather(main_network->chunk_nodes, main_network->chunk_nodes_len, MPI_INT, chunk_nodes_buffer, main_network->chunk_nodes_len, MPI_INT, 0, MPI_COMM_WORLD);

	// Uniqueness of partition check
	if (world_rank == 0) {
		int chunk_sum = 0;
		for (int i = 0; i < main_network->chunk_nodes_len * world_size; i++) {
			if (chunk_nodes_buffer[i] != -1) {
				chunk_sum += chunk_nodes_buffer[i];
			}
		}
		int nn = main_network->n_nodes;
		int id;
		for(int d = 0 ; d<n_chunk_edges; d++){
			id = main_network->chunk_edges[d];
			if(main_network->edges[2*id] >= nn || main_network->edges[2*id + 1] >= nn){
						cout<<"Node is "<<main_network->edges[2*id]<<" for index "<<id<<endl;
				}	
		}
		if (chunk_sum != (nn*nn - nn)/2 ) {
			cout << chunk_sum << " | "<< (nn*nn - nn)/2<<endl; 
			cout << "Uneven chunk partitioning" << endl;

		}
	}
	
	cout << "World rank proc "<<world_rank << " starting the loop:" << endl;

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
		// main_network->plotNetwork(iter, false);
		main_network->optimize();
		MPI_Barrier(MPI_COMM_WORLD);



		//TODO: check the size of array parameter
		MPI_Gather(main_network->R, main_network->n_nodes * DIM, MPI_FLOAT, R_buffer, main_network->n_nodes * DIM, MPI_FLOAT, 0, MPI_COMM_WORLD);
		if((iter+1)%NSYNC == 0){
			MPI_Gather(main_network->forces, main_network->n_nodes * DIM, MPI_FLOAT, forces_buffer, main_network->n_nodes * DIM, MPI_FLOAT, 0, MPI_COMM_WORLD);
		}
		// syncing R
		if (world_rank == 0) {
			int node_to_sync  = 0;
			for (int i = 0; i < world_size; i += 1) {
				for (int j = i*main_network->chunk_nodes_len; j < (i+1)*main_network->chunk_nodes_len; j++) {
					node_to_sync = chunk_nodes_buffer[j];
					if (node_to_sync == -1) {
						break;
					}
					else{
						main_network->R[DIM * node_to_sync] = R_buffer[main_network->n_nodes * DIM * i + DIM * node_to_sync];
						main_network->R[DIM * node_to_sync + 1] = R_buffer[main_network->n_nodes * DIM * i + DIM * node_to_sync + 1];
						if((iter+1)%NSYNC == 0){
							main_network->forces[DIM * node_to_sync] = forces_buffer[main_network->n_nodes * DIM * i + DIM * node_to_sync];
							main_network->forces[DIM * node_to_sync + 1] = forces_buffer[main_network->n_nodes * DIM * i + DIM * node_to_sync + 1];
						}
					}
				}
				// int i = (main_network->R[i]/chunk_size);
				// main_network->R[i] = R_buffer[main_network->n_nodes * DIM * i + i];
				// main_network->forces[i] = forces_buffer[main_network->n_nodes * DIM * i + i];
			}
			
			main_network->move_top_plate();

		}
		
		MPI_Bcast(main_network->R, main_network->n_nodes * DIM, MPI_FLOAT, 0, MPI_COMM_WORLD);

		if (world_rank == 0) {
			main_network->get_plate_forces(plate_forces, iter);
			//main_network->get_stats();
			
		}
		//MPI_Barrier(MPI_COMM_WORLD);
	}

	if (world_rank == 0) {
		string file_name = "forcesMPI.txt";
		write_to_file<float>(file_name, plate_forces, STEPS, DIM);
		free(plate_forces);
	}

	free(R_buffer);
	cout << "Made it to the end!" << endl;
	MPI_Finalize();
	return 1;
}
