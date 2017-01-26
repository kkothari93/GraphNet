#define SPCL_NUM -1.0
#define BLOCK_SIZE 1024


struct hostvars{
	// all the host variables that need to memcpy'd
	// to device need to be placed here
	const float* PBC_vector;
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
};

__global__ void optimize_cuda(float*, int*, float*, float*, \
	const float*, const int, const int, \
	const bool*, const float*, const int*, int, const int*, int, \
	float* , int , float , float, int);

__device__ float kfe_cuda(float);

__device__ float force_wlc_cuda(float, float);

__device__ bool notmember(int, const int*, int, const int*, int);

void pull_CUDA(hostvars*, int);