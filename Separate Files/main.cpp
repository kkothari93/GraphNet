#include <iostream>
#include <cmath>
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
#include <mpi.h>

#include "mpi_network.h"
using namespace std;

//#define USE_MPI false


int main() {

	//string path = "/media/konik/Research/2D sacrificial bonds polymers/cpp11_code_with_cuda/template2d.msh";
	string path = "./template2d_z4.msh";
	#if SACBONDS
	#define DECL_NET sacNetwork test_network(path)
	#else
	#define DECL_NET Network test_network(path)
	#endif
	DECL_NET;

	cout << "Do you want to use MPI? Input \'y\' to use MPI, else enter any other character." << endl;
	char USE_MPI;
	cin >> USE_MPI;

	if (USE_MPI != 'y') {
		float weight_goal = 1.03754e6; // weight of similarly sized triangular mesh network
	
		float weight_multiplier;
		float weight = test_network.get_weight();
		if (weight<weight_goal){weight_multiplier = test_network.set_weight(weight_goal);}
		bool should_stop = test_network.get_stats();
		int old_n_edges = test_network.get_current_edges();
		int curr_n_edges = old_n_edges;


		if(should_stop){return 0;}
		float* plate_forces;
		plate_forces = (float*)malloc(sizeof(float)*DIM*STEPS);
		memset(plate_forces, 0.0, STEPS*DIM*sizeof(*plate_forces));

		if(CRACKED){
			Cracklist alist(4, MAXBOUND);
			test_network.apply_crack(alist);
		}


		test_network.plotNetwork(0, true);

		clock_t t = clock();
		cout<<"\n Will run for "<<STEPS<<":\n";
		
		for(int i = 0; i<STEPS; i++ ){
			
			test_network.optimize();
			test_network.move_top_plate();
			test_network.get_plate_forces(plate_forces, i);
			if((i+1)%100 == 0){
				should_stop = test_network.get_stats();
				if(should_stop){break;}
				curr_n_edges = test_network.get_current_edges();
				if(curr_n_edges<old_n_edges){
						test_network.plotNetwork(i, false);
				}
				cout<<"Step "<<(i+1)<<" took "<<float(clock()-t)/CLOCKS_PER_SEC<<" s\n";
				t = clock();  // reset clock
			}

		}
		string sb = SACBONDS ? "true" : "false" ; 
		string fname = FNAME_STRING + std::to_string(L_STD/L_MEAN) + "_" + sb + ".txt";
		write_to_file<float>(fname, plate_forces, STEPS, DIM);

		free(plate_forces);
	}
	else {

		MPI_Init(NULL, NULL);
	
		// Get the number of processes
	  	int world_size;
	  	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	  	// Get the rank of the process
	  	int world_rank;
	  	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	  	MPI_Network * main_network = new MPI_Network(test_network);

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

		main_network->init_MPI(world_rank, world_size);
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
			for(int d = 0 ; d<main_network->n_chunk_edges; d++){
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
			cout << world_rank << " " << __LINE__ << " iter: " << iter << endl;
			main_network->plotNetwork(iter, false);
			main_network->optimize();
			MPI_Barrier(MPI_COMM_WORLD);
			//cout <<  world_rank<< "  " <<__LINE__ << endl;

			//MPI_Barrier(MPI_COMM_WORLD);
			MPI_Gather(main_network->R, main_network->n_nodes * DIM, MPI_FLOAT, R_buffer, main_network->n_nodes * DIM, MPI_FLOAT, 0, MPI_COMM_WORLD);
			if((iter+1)%NSYNC == 0){
				MPI_Gather(main_network->forces, main_network->n_nodes * DIM, MPI_FLOAT, forces_buffer, main_network->n_nodes * DIM, MPI_FLOAT, 0, MPI_COMM_WORLD);
				if (world_rank == 0) {
					cout << "Synced forces" << endl;
				}
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
		//cout<<__LINE__<<endl;
		//sync forces at the end??
		MPI_Gather(main_network->forces, main_network->n_nodes * DIM, MPI_FLOAT, forces_buffer, main_network->n_nodes * DIM, MPI_FLOAT, 0, MPI_COMM_WORLD);
		if (world_rank == 0) {
			int node_to_sync  = 0;
			for (int i = 0; i < world_size; i += 1) {
				for (int j = i*main_network->chunk_nodes_len; j < (i+1)*main_network->chunk_nodes_len; j++) {
					node_to_sync = chunk_nodes_buffer[j];
					if (node_to_sync == -1) {
						break;
					}
					else{
						main_network->forces[DIM * node_to_sync] = forces_buffer[main_network->n_nodes * DIM * i + DIM * node_to_sync];
						main_network->forces[DIM * node_to_sync + 1] = forces_buffer[main_network->n_nodes * DIM * i + DIM * node_to_sync + 1];
					}
				}
			}
			main_network->get_plate_forces(plate_forces, STEPS);
			//main_network->move_top_plate();
		}
		


		if (world_rank == 0) {
			string file_name = "forcesMPI.txt";
			write_to_file<float>(file_name, plate_forces, STEPS, DIM);
			free(plate_forces);
			plate_forces = NULL;
		}
		// Needed to not have double free, corruption error
		free(R_buffer);
		R_buffer = NULL;
		free(forces_buffer);
		forces_buffer = NULL;
		free(chunk_nodes_buffer);
		chunk_nodes_buffer = NULL;
		delete main_network;
		main_network = NULL;
		cout << "Made it to the end! Exiting Now." << endl;
		MPI_Finalize();

	}

	return 0;
}
