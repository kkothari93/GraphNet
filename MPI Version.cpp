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

<<<<<<< HEAD
	//uniform L and PBC across all processors
	MPI_Bcast(main_network->L, main_network->n_elems, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(main_network->PBC, main_network->n_elems, MPI_C_BOOL, 0, MPI_COMM_WORLD); //not sure if it is randomly generated.
	float * forces_buffer = new float[main_network->n_nodes * DIM * world_size];
	float * recv_buffer = new float[main_network->n_nodes * DIM * world_size]; //buffer to gather the R from all nodes
	bool * broken_buffer = new bool[world_size];
	int * edges_buffer = new int[main_network->n_elems * DIM * world_size];
=======

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

>>>>>>> BW
	int iter = 0; // needed to write forces later

	for(iter = 0; iter<STEPS; iter++){
		if((iter+1)%1000 == 0){ // +1 required to have values in p_x, p_y
			cout<<(iter+1)<<endl; 
			cout<<"That took "<<(clock()-t)/CLOCKS_PER_SEC<<" s\n";
			t = clock();  // reset clock

<<<<<<< HEAD
		main_network->optimize(BROKEN, lo, hi);
		//TODO: use Network::get_plate_forces() on root proc
		MPI_Gather(main_network->forces, main_network->n_nodes * DIM, MPI_FLOAT, forces_buffer, main_network->n_nodes * DIM, MPI_FLOAT, 0, MPI_COMM_WORLD);


		MPI_Barrier(MPI_COMM_WORLD);


		//TODO: check the size of array parameter
		MPI_Gather(main_network->R, main_network->n_nodes * DIM, MPI_FLOAT, recv_buffer, main_network->n_nodes * DIM, MPI_FLOAT, 0, MPI_COMM_WORLD);
		
		if (world_rank == 0) {
			for (int i = 0; i < main_network->n_nodes * DIM; i++) {
				int proc_to_copy_from = (main_network->R[i]/chunk_size);
				main_network->R[i] = recv_buffer[main_network->n_nodes * DIM * proc_to_copy_from + i];
				main_network->forces[i] = forces_buffer[main_network->n_nodes * DIM * proc_to_copy_from + i];
			}
		}
		
		MPI_Bcast(main_network->R, main_network->n_nodes * DIM, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Bcast(main_network->forces, main_network->n_nodes * DIM, MPI_FLOAT, 0, MPI_COMM_WORLD);
		// if (world_rank == 0) {
		// 	for (int i = 0; i < main_network->n_nodes; i++) {
		// 		printf("%d:\t%f\t%f\n",i, main_network->forces[2*i],main_network->forces[2*i+1]);
		// 	}
		// }
		//checking if edges need to be resynced
		MPI_Allgather(&BROKEN, 1, MPI_C_BOOL, broken_buffer, 1, MPI_C_BOOL, MPI_COMM_WORLD);
		//cout << __LINE__ << endl;
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

		if (world_rank == 0) {
			//move_top_nodes(R, tsideNodes, n_tside);
			main_network->move_top_plate();
			main_network->get_plate_forces(plate_forces, iter);
			main_network->get_stats();
			
		}

		MPI_Barrier(MPI_COMM_WORLD);
=======
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

>>>>>>> BW
	}

<<<<<<< HEAD
	if (world_rank == 0) {
		string file_name = "forcesMPI.txt";
		write_to_file<float>(file_name, plate_forces, STEPS, DIM);
		free(plate_forces);
	}
=======
>>>>>>> BW















  	MPI_Finalize();

}