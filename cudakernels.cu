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


#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

__device__ float kfe_cuda(float force){
	return ae*expf(force*delxe/kB/T)*TIME_STEP;
}

__device__ float force_wlc_cuda(float x, float L){
	float t = x / L;
	if (t < 0.99){ return kB*T / b * (t + 1.0 / 4.0 / powf((1.0 - t), 2.0) - 1.0 / 4.0); }
	else { return 999999.0; }
}

__global__ void optimize_cuda(float* R, int* edges, float* damage_integral, float* forces, \
	const float* chain_len, const int num_nodes, const int num_edges, \
	const bool* PBC_STATUS, const float* PBC_vector, \
	const int* tnodes, int n_tnodes, const int* moving_nodes, int n_moving, \
	float* plate_force, int n_steps,\
	float eta = 0.01, float alpha = 0.9, int max_iter_opt = 1000){

	int tid = threadIdx.x + blockIdx.x*BLOCK_SIZE;
	float pbc_v[2] = {PBC_vector[0], PBC_vector[1]}; 
	float rms_history[2];
	float delR[2];
	float grad[2];

	int pair, n1, n2, n_t;
	float L, x1, x2, y1, y2, dist, force;
	float unitvector[DIM];
	bool pbc;

	//if(tid==1){printf}

	if(tid<num_edges){
		damage_integral[tid] = 0.0;
		pbc = PBC_STATUS[tid];

		// Assign threads to edges
		pair = tid * 2;
		L = chain_len[tid];

		// read the nodes that the thread has been assigned
		n1 = edges[pair];
		n2 = edges[pair+1];
	}
	__syncthreads();

	for(int iter = 0; iter<n_steps; iter++){
		
		// reset grads for next iteration
		rms_history[0] = 0.0;
		rms_history[1] = 0.0;
		delR[0] = 0.0;
		delR[1] = 0.0;
		grad[0] = 0.0;
		grad[1] = 0.0;

		for(int step = 0; step < max_iter_opt; step++){
			///////////////////////////////////////////////
			// Force calculations
			// 
			// Here each edge is assigned one thread
			//
			///////////////////////////////////////////////
	
			// zero all forces
			if(tid<num_nodes*DIM){
				forces[tid] = 0.0;
			}
			__syncthreads();

			force = 0.0;
					
			if(tid < num_edges){
			// Check if connection exists 
			if(n1 != SPCL_NUM && n2 != SPCL_NUM){

				// read the positions of the crosslinkers
				x1 = R[n1*DIM];
				y1 = R[n1*DIM + 1];
				x2 = R[n2*DIM];
				y2 = R[n2*DIM + 1];
				//if(step==0 && iter==0){printf("(%f, %f); (%f, %f)\n", x1, y1, x2, y2);}
				
				// Calculate distance, unit vector and force
				// Shared memory is per block. If num_edges*DIM is too large each block
				// can be held responsible for separate pairs and then atomic adds can 
				// be done. Another approach could be to have each thread implement force
				// calc for one node to avoid atomic adds but that requires each thread to 
				// run through all edges and figure out which ones to add. That will be order
				// n whereas atomic adds should be order z extra work

				// check for PBC;
				if(pbc==true){
					dist = hypotf(x1-x2-pbc_v[0], y1-y2-pbc_v[1]);
					// add unitvector 3 for DIM 3. __in future use for loop here
					unitvector[0] = (x1 - x2 - pbc_v[0])/dist;
					unitvector[1] = (y1 - y2 - pbc_v[1])/dist;
				}
				else{
					dist = hypotf(x1-x2, y1-y2);
					// add unitvector 3 for DIM 3. __in future use for loop here
					unitvector[0] = (x1 - x2)/dist;
					unitvector[1] = (y1 - y2)/dist;
				}
			
				
				// calculate force
				force = force_wlc_cuda(dist, L);

				// Break if force too high
				if(force==999999){
					//printf("Breaking bond between %d and %d at iter %d, step %d\n", n1, n2, iter, step);
					n1 = SPCL_NUM;
					n2 = SPCL_NUM;
					edges[pair] = n1;
					edges[pair + 1] = n2;
					force = 0.0;
					damage_integral[tid] = 1.1;
				}
				else{
					// add the forces calculated to the nodes (atomic add)
					atomicAdd(&forces[n1*DIM], -1.0*force*unitvector[0]);
					atomicAdd(&forces[n1*DIM+1], -1.0*force*unitvector[1]);

					atomicAdd(&forces[n2*DIM], force*unitvector[0]);
					atomicAdd(&forces[n2*DIM+1], force*unitvector[1]);
				}
				
			}}
			__syncthreads();
			
			/////////////////////////////////////////////////
			//
			// Optimization step
			//
			/////////////////////////////////////////////////

			// Assign each thread to nodes*DIM
			if(tid<n_moving){
				n_t = moving_nodes[tid];
				grad[0] = forces[n_t*DIM];
				grad[1] = forces[n_t*DIM + 1];

				rms_history[0] = alpha*rms_history[0] + (1-alpha)*grad[0]*grad[0];
				rms_history[1] = alpha*rms_history[1] + (1-alpha)*grad[1]*grad[1];

				delR[0] = eta/__frsqrt_rn((rms_history[0] + 1.0e-6)) * grad[0];
				delR[1] = eta/__frsqrt_rn((rms_history[1] + 1.0e-6)) * grad[1];
				
				R[n_t*DIM] += delR[0];
				//if(fabs(delR[0])>10.0){printf("For node %d we have forces %0.3f\n",n_t, grad[0] );}
				R[n_t*DIM + 1] += delR[1];
				//if(fabs(delR[1])>10.0){printf("For node %d we have forces %0.3f\n",n_t, grad[1] );}
			}
			__syncthreads();
		}


		// Check if connection exists 
		if(tid<num_edges){
			if(damage_integral[tid]<1.0){
				// Update damage integral
				damage_integral[tid] += kfe_cuda(force)*TIME_STEP;
			}
			// Update edges acc. to damage
			if(damage_integral[tid] > 1.0){
				damage_integral[tid] = 1.1;
				n1 = SPCL_NUM;
				n2 = SPCL_NUM;
				edges[pair] = SPCL_NUM;
				edges[pair + 1] = SPCL_NUM;
			}
		}

		// update the force in the array
		if(tid<n_tnodes){
			n_t = tnodes[tid];
			atomicAdd(&plate_force[iter*DIM], forces[DIM*n_t]);
			atomicAdd(&plate_force[iter*DIM + 1], forces[DIM*n_t + 1]);
		}

		// move top nodes acc. to velocity
		if(tid < 2*n_tnodes && tid >= n_tnodes){
			n_t = tnodes[tid - n_tnodes];
			R[DIM*n_t] += vel[0]*TIME_STEP;
			R[DIM*n_t + 1] += vel[1]*TIME_STEP;
		}
		__syncthreads();
	
	}
}


void sanity_check(hostvars* vars){
	int n1, n2,c=0;
	float dist, x1, x2, y1, y2, L;
	float pbc_v[2] = {vars->PBC_vector[0], vars->PBC_vector[1]};
	for(int i=0; i<vars->n_elems; i++){
		n1 = vars->edges[2*i];
		n2 = vars->edges[2*i+1];
		x1 = vars->R[n1*DIM];
		y1 = vars->R[n1*DIM + 1];
		x2 = vars->R[n2*DIM];
		y2 = vars->R[n2*DIM + 1];
		if(vars->PBC[i] == true){x2 += pbc_v[0]; y2 += pbc_v[1]; }
		dist = sqrtf((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
		L = vars->L[i];
		if (dist>=L){
			printf("(%0.2f, %0.2f); (%0.2f, %0.2f) and L = %0.2f\n", x1, y1, x2, y2, L);
			c += 1;
		}
		}
	printf("%d of %d elements make no sense!", c, vars->n_elems);
}

void pull_CUDA(hostvars* vars, int n_steps){
	// Pass all host variables in a struct
	
	// Initialize nodes and edges
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
	
	//sanity_check(vars);

	// Check if the transfers are ok
	printf("n_nodes: %d \t", n_nodes);
	printf("n_elems: %d \t", n_elems);
	printf("n_tnodes: %d \t", n_tside);
	printf("n_moving: %d \n", n_moving);

	// GPU allocations
	cudaMalloc((void**)&R_d, n_nodes*DIM*sizeof(float));
	cudaMalloc((void**)&forces_d, n_nodes*DIM*sizeof(float));
	cudaMalloc((void**)&edges_d, 2*n_elems*sizeof(int));
	cudaMalloc((void**)&tsideNodes_d, n_tside*sizeof(int));
	cudaMalloc((void**)&moving_nodes_d, n_moving*sizeof(int));
	cudaMalloc((void**)&PBC_d, n_elems*sizeof(bool));
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
	cudaMemcpy(L_d, vars->L, n_elems*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(damage_d, vars->damage, n_elems*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(pull_forces_d, vars->pull_forces, n_steps*DIM*sizeof(float), cudaMemcpyHostToDevice);
	printf("Transfer successful! \n");

	// Define grid and block size
	// Launch atleast as many threads as edges
	dim3 gridsize((n_elems-1)/BLOCK_SIZE + 1);
	dim3 blocksize(BLOCK_SIZE);

	// Launch timer code
	clock_t t = clock();

	int n_e = 0;
	for(int i=0; i< n_elems; i++){
		if(vars->edges[2*i]!=SPCL_NUM){
			n_e++;
		}
	}
	printf("We have %d edges\n",n_e);

	printf("Launching kernel...\n");
	optimize_cuda<<< gridsize, blocksize >>>(
		R_d, edges_d, damage_d, forces_d, \
		L_d, n_nodes, n_elems, PBC_d, \
		PBC_vector_d, tsideNodes_d, n_tside, \
		moving_nodes_d, n_moving, \
		pull_forces_d, n_steps);
	cudaDeviceSynchronize();
	//gpuErrchk(cudaDeviceSynchronize());
	printf("%d steps took %0.5f s\n", n_steps, float(clock()-t)/CLOCKS_PER_SEC);

	// Copy device to host
	cudaMemcpy(vars->R, R_d,  n_nodes*DIM*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(vars->forces, forces_d,  n_nodes*DIM*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(vars->damage, damage_d,  n_elems*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(vars->pull_forces, pull_forces_d,  n_steps*DIM*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(vars->edges, edges_d,  2*n_elems*sizeof(int), cudaMemcpyDeviceToHost);

	n_e = 0;
	for(int i=0; i< n_elems; i++){
		if(vars->edges[2*i]!=SPCL_NUM || vars->edges[2*i+1]!=SPCL_NUM){
			n_e++;
		}
	}
	printf("We have %d edges\n",n_e);

	// for(int i=0; i<n_elems; i++){
	// 	if (vars->damage[i]>=1.0){
	// 		printf("damage[%d]\t%0.5f\n",i, vars->damage[i]);}
	// 	// for(int d=0; d<DIM; d++){
	// 	// 	printf("%0.3f\t",vars->R[i*DIM + d]);
	// 	// }
	// 	// printf("\n");
	// }

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
