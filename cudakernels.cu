#include "cudakernels.h"
#include <stdio.h>
#include "vel.h"
#include <time.h>

#ifndef __params__
#define Z_MAX 10
#define DIM 2
#endif

#ifndef __constants__
#define kB 1.38064852e-5					// Boltzmann constant
#define b 0.1								// Persistence length
#define T 300 								// Temperature
#define ae 0.1 								// Strength of bond - includes activation energy
#define delxe 0.15 							// parameter for breaking crosslink connection
#define TIME_STEP 1e-4
#endif

__device__ static float vel[2] = {vel_x, vel_y};


__device__ float kfe_cuda(float force){
	return expf(force*delxe/kB/T)*TIME_STEP;
}

__device__ float force_wlc_cuda(float x, float L){
	float t = x / L;
	if (t < 0.99){ return kB*T / b * (t + 1.0 / 4.0 / powf((1.0 - t), 2.0) - 1.0 / 4.0); }
	else { return 999999.0; }
}

__device__ bool notmember(int id, const int* bnodes, int n_bnodes, \
	const int* tnodes, int n_tnodes){
	bool out = false;
	for(int i=0;i<max(n_tnodes, n_bnodes);i++){
		if(i<n_tnodes){
			if(id==tnodes[i]){out=true; break;}
		}
		if(i<n_bnodes){
			if(id==bnodes[i]){out=true; break;}
		}
	}
	return out;
}

__global__ void optimize_cuda(float*R, int* edges, float* damage_integral, float* forces, \
	const float* chain_len, const int num_nodes, const int num_edges, \
	const bool* PBC_STATUS, const float* PBC_vector, \
	const int* tnodes, int n_tnodes, const int* moving_nodes, int n_moving, \
	float* plate_force, int n_steps,\
	float eta = 0.1, float alpha = 0.9, int max_iter_opt = 1000){

	// Get indices
	int tx = threadIdx.x; 
	int bx = blockIdx.x;
	int tid = tx + bx*BLOCK_SIZE;
	

	float rms_history[2] = {0.0, 0.0};
	float delR[2] = {0.0, 0.0};
	float grad[2] = {0.0, 0.0};

	int pair, edge_num, n1, n2, n_t;
	float L, x1, x2, y1, y2, dist, force, diss_energy, top_force_x, top_force_y;
	float unitvector[DIM];
	for(int iter = 0; iter<n_steps; iter++){
		for(int step = 0; step < max_iter_opt; step++){
			///////////////////////////////////////////////
			// Force calculations
			// 
			// Here each edge is assigned one thread
			//
			///////////////////////////////////////////////

			// Assign threads to edges
			pair = tid * 2;
			edge_num = tid; 
			if(tid<num_edges){
			L = chain_len[edge_num];

			// read the nodes that the thread has been assigned
			n1 = edges[pair];
			n2 = edges[pair+1];

			// Check if connection exists 
			if(n1!=SPCL_NUM || n2 != SPCL_NUM){
				// read the positions of the crosslinkers
				x1 = R[n1*DIM];
				y1 = R[n1*DIM + 1];
				x2 = R[n2*DIM];
				y2 = R[n2*DIM + 1];

				// Calculate distance, unit vector and force
				// Shared memory is per block. If num_edges*DIM is too large each block
				// can be held responsible for separate pairs and then atomic adds can 
				// be done. Another approach could be to have each thread implement force
				// calc for one node to avoid atomic adds but that requires each thread to 
				// run through all edges and figure out which ones to add. That will be order
				// n whereas atomic adds should be order z extra work

				// check for PBC;
				if(PBC_STATUS[edge_num]==true){
					dist = hypot(x1-x2-PBC_vector[0], y1-y2-PBC_vector[1]);
				}
				else{
					dist = hypot(x1-x2, y1-y2);
				}

				// add unitvector 3 for DIM 3. __in future use for loop here
				unitvector[0] = (x1 - x2)/dist;
				unitvector[1] = (y1 - y2)/dist;


				// calculate force
				force = force_wlc_cuda(dist, L);

				// zero all forces
				if(tid<num_nodes*DIM){
					forces[tid] = 0.0;
				}

				//required before next step
				__syncthreads();

				// add the forces calculated to the nodes (atomic add)
				atomicAdd(&forces[n1], force*unitvector[0]);
				atomicAdd(&forces[n1+1], force*unitvector[1]);

				atomicAdd(&forces[n2], force*unitvector[0]);
				atomicAdd(&forces[n2+1], force*unitvector[1]);

				__syncthreads();
			}

			/////////////////////////////////////////////////
			//
			// Optimization step
			//
			/////////////////////////////////////////////////

			// Assign each thread to nodes*DIM
			if(tid<n_moving){
				n1 = moving_nodes[tid];
				grad[0] = forces[n1*2];
				grad[1] = forces[n1*2 + 1];
				
				rms_history[0] = alpha*rms_history[0] + (1-alpha)*grad[0]*grad[0];
				rms_history[1] = alpha*rms_history[1] + (1-alpha)*grad[1]*grad[1];

				delR[0] = eta*__frsqrt_rn(1.0/(rms_history[0] + 1.0e-6)) * grad[0];
				delR[1] = eta*__frsqrt_rn(1.0/(rms_history[1] + 1.0e-6)) * grad[1];
				
				R[n1*2] += delR[0];
				R[n1*2 + 1] += delR[1];
			}
			__syncthreads();
		}

		// Update damage integral
		diss_energy = kfe_cuda(force)*TIME_STEP;
		damage_integral[edge_num] += diss_energy;

		// Update edges acc. to damage
		if(damage_integral[edge_num] > 1.0){
			edges[pair] = SPCL_NUM;
			edges[pair + 1] = SPCL_NUM;
		}

		// update the force in the array
		if(tid<n_tnodes){
			n_t = tnodes[tid];
			top_force_x = forces[DIM*n_t];
			top_force_y = forces[DIM*n_t + 1];
			atomicAdd(&plate_force[iter*DIM], top_force_x);
			atomicAdd(&plate_force[iter*DIM + 1], top_force_y);
		}

		// move top nodes acc. to velocity
		if(tid < 2*n_tnodes && tid >= n_tnodes){
			n_t = tnodes[tid - n_tnodes];
			R[DIM*n_t] += vel[0]*TIME_STEP;
			R[DIM*n_t + 1] += vel[1]*TIME_STEP;
		}
		__syncthreads();
	}	
	//	if(iter%500 == 0 && tid == 1){
	//		printf("Completed %d iterations...\n",iter);
	//		printf("%d took %0.5f s\n", iter, float(clock()-t)/CLOCKS_PER_SEC);
	//		t = clock();
	//	}
	}
}


void pull_CUDA(hostvars* vars, int n_steps){
	// Pass all host variables in a struct
	
	// Initialize nodes and edges
	printf("Got in! \n");
	float* R_d; int* edges_d;
	float* forces_d;
	int* tsideNodes_d;
	int* moving_nodes_d;
	bool* PBC_d;
	float* L_d; float* damage_d;
	float* pull_forces_d;
	float* PBC_vector_d;
	int n_nodes = vars->n_nodes;
	int n_elems = vars->n_elems;
	int n_tside = vars->n_tnodes;
	int n_moving = vars->n_moving;
	
	// Check if the transfers are ok
	printf("n_nodes: %d \t", n_nodes);
	printf("n_elems: %d \t", n_elems);
	printf("n_tnodes: %d \t", n_tside);
	printf("n_moving: %d \n", n_moving);

	// GPU allocations
	cudaMalloc((void**)&R_d, n_nodes*DIM*sizeof(float));
	cudaMalloc((void**)&forces_d, n_nodes*DIM*sizeof(float));
	cudaMalloc((void**)&edges_d, Z_MAX*n_nodes*2*sizeof(int));
	cudaMalloc((void**)&tsideNodes_d, n_tside*sizeof(int));
	cudaMalloc((void**)&moving_nodes_d, n_moving*sizeof(float));
	cudaMalloc((void**)&PBC_d, 2*n_elems*sizeof(bool));
	cudaMalloc((void**)&PBC_vector_d, DIM*sizeof(float));
	cudaMalloc((void**)&L_d, n_elems*sizeof(float));
	cudaMalloc((void**)&damage_d, n_elems*sizeof(float));
	cudaMalloc((void**)&pull_forces_d, n_steps*DIM*sizeof(float));
	printf("Malloc successful! \n");

	// Copy host to device
	cudaMemcpy(PBC_vector_d, vars->PBC_vector, DIM*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(PBC_d, vars->PBC, n_elems*sizeof(bool), cudaMemcpyHostToDevice);
	cudaMemcpy(R_d, vars->R, n_nodes*DIM*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(forces_d, vars->forces, n_nodes*DIM*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(edges_d, vars->edges, 2*n_elems*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(moving_nodes_d, vars->moving_nodes, n_moving*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(tsideNodes_d, vars->tsideNodes, n_tside*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(L_d, vars->L, 2*n_elems*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(damage_d, vars->damage, 2*n_elems*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(pull_forces_d, vars->pull_forces, n_steps*DIM*sizeof(float), cudaMemcpyHostToDevice);
	printf("Transfer successful! \n");
	
	// Define grid and block size
	// Launch atleast as many threads as edges
	dim3 gridsize((n_elems-1)/BLOCK_SIZE + 1);
	dim3 blocksize(BLOCK_SIZE);

	// Launch timer code
	clock_t t = clock();

	printf("Launching kernel...\n");
	optimize_cuda<<< gridsize, blocksize >>>(
		R_d, edges_d, damage_d, forces_d, \
		L_d, n_nodes, n_elems, PBC_d, \
		PBC_vector_d, tsideNodes_d, n_tside, \
		moving_nodes_d, n_moving, \
		pull_forces_d, n_steps);
	cudaDeviceSynchronize();
	
	printf("%d steps took %0.5f s\n", n_steps, float(clock()-t)/CLOCKS_PER_SEC);

	// Copy device to host
	cudaMemcpy(vars->R, R_d,  n_nodes*DIM*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(vars->forces, forces_d,  n_nodes*DIM*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(vars->damage, damage_d,  2*n_elems*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(vars->pull_forces, pull_forces_d,  n_steps*DIM*sizeof(float), cudaMemcpyDeviceToHost);

	// Free up global memory
	cudaFree(R_d);
	cudaFree(forces_d);
	cudaFree(edges_d);
	cudaFree(moving_nodes_d);
	cudaFree(tsideNodes_d);
	cudaFree(PBC_d);
	cudaFree(PBC_vector_d);
	cudaFree(L_d);
	cudaFree(damage_d);
	cudaFree(pull_forces_d);

	return;
}
