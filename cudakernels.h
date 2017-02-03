#define SPCL_NUM -1.0
#define BLOCK_SIZE 1024


struct hostvars{
	// all the host variables that need to memcpy'd
	// to device need to be placed here
	const float* PBC_vector;
<<<<<<< HEAD
	float* R_d; int* edges_d;
	float* forces_d;
	int* bsideNodes_d; int* tsideNodes_d;
	int* lsideNodes_d; int* rsideNodes_d;
	bool* PBC_d;
	float* L_d; float* damage_d;
	float* pull_forces_d;
	const float* PBC_vector_d;
=======
	float* R; int* edges;
	float* forces;
	int* bsideNodes; int* tsideNodes;
	int* lsideNodes; int* rsideNodes;
	bool* PBC;
	float* L; float* damage;
	float* pull_forces;
	int n_nodes, n_elems;
	int n_tnodes, n_bnodes;
	int n_side_nodes;
>>>>>>> 3c983a41043bd95e04788ba87070f7e85d7aa2f1
};

__global__ void optimize_cuda(float*, int*, float*, float*, \
	const float*, const int, const int, \
	const bool*, const float*, const int*, int, const int*, int, \
	float* , int , float , float, int);

__device__ float kfe_cuda(float);

__device__ float force_wlc_cuda(float, float);

__device__ bool notmember(int, const int*, int, const int*, int);

void pull_CUDA(hostvars*, int);