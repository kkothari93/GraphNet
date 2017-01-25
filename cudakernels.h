#define SPCL_NUM -1.0
#define BLOCK_SIZE 1024

struct hostvars{
	// all the host variables that need to memcpy'd
	// to device need to be placed here
	const float* PBC_vector;
	float* R_d; int* edges_d;
	float* forces_d;
	int* bsideNodes_d; int* tsideNodes_d;
	int* lsideNodes_d; int* rsideNodes_d;
	bool* PBC_d;
	float* L_d; float* damage_d;
	float* pull_forces_d;
	const float* PBC_vector_d;
}

void optimize(float*, int*, float*, float*, \
	const float*, const int, const int, \
	const bool*, const float*, const int*, int, const int*, int, \
	float* , int , float , float, int);

float kfe_cuda(float);

float force_wlc_cuda(float, float);

bool notmember(int, const int*, int, const int*, int);