#ifndef __helper__
#define __helper__
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
#include "vel.h"
//#include "gnuplot_i.hpp"
using namespace std;

#define __params__	
#define DIM 2								// Number of dimensions
#define TIME_STEP 1e-4						// Time step
#define SIM_TIME 9.0						// Simulation time
#define TOL 1e-6							// Tolerance
#define STEPS int(SIM_TIME/TIME_STEP)		// Number of time steps
#define L_MEAN 250.0f						// Average for contour length
#define L_STD 150.0f						// Std. deviation for contour lengths
#define MAXBOUND 500.0f
#define PAD MAXBOUND*1.03
#define SACBONDS false
#define IMPLEMENT_PBC true
#define FNAME_STRING "high_disorder_high_L_"
#define CRACKED false
//#define GENERATOR std::uniform_real_distribution<float>(L_MEAN - L_STD, L_MEAN + L_STD
#if CRACKED
#define PROB_REMOVAL 0.8
#else
#define PROB_REMOVAL 0.0
#endif


// Define constants
#define __constants__

#define kB 1.38064852e-5					// Boltzmann constant
#define b_poly 0.1								// Persistence length
#define T 300 								// Temperature
#define ae 0.1 								// Strength of bond - includes activation energy
#define delxe 0.15 							// parameter for breaking crosslink connection
#define af 0.3
#define delxf 0.25

void side_nodes(float* R,\
	int* lnodes, int* rnodes, int* tnodes, int* bnodes,\
	int& n_lside, int& n_rside, int& n_tside, int& n_bside, int n);

inline int get_num_vertices(int elem_type){
	switch(elem_type){
		//case 1: return 2;
		case 2: return 3;
		case 3: return 4;
		//case 4: return 4;
		//case 5: return 8;
		//case 7: return 5;
		default: return -1;
	}
}

inline bool contains(vector<int>& vec, int elem){
	return (std::find(vec.begin(), vec.end(), elem) != vec.end());
}

inline float getnorm(const float* vec, const int dim = DIM){
	float s = 0,t;
	#pragma unroll
	for (int j = 0; j<DIM; j++){
		t = vec[j];
		s += t*t;
	}
	return sqrt(s);
}

void mapping(int& edge_counter, int elem_type, int n_vertices, stringstream& input, int* edges);

void mapping(int& edge_counter, int elem_type);

void read_n(int& n_nodes, int& n_elems, string& fname);

void take_input(float* R, int* edges, int n_nodes, int n_elems, string& fname);

void __init__(float* L, float* damage, bool* PBC, int n_elems);

inline void normalize_vector(float* result, const float* vec){
	float norm = getnorm(vec);
	#pragma unroll
	for (int i = 0; i<DIM; i++){
		result[i] = vec[i] / norm;
	}
}

inline bool does_file_exist(string& fname){
	ifstream infile(fname);
	return infile.good();
}

inline void normalize_vector(float* vec){
	float norm = getnorm(vec);
	#pragma unroll
	for (int i = 0; i<DIM; i++){
		vec[i] = vec[i] / norm;
	}
}

inline void unitvector(float* result, float* r1, float* r2){
	#pragma unroll
	for (int j = 0; j<DIM; j++){
		result[j] = r1[j] - r2[j];
	}
	normalize_vector(result);
}

inline float force_wlc(float x, float L){
	float t = x / L;
	if (t < 0.99){ return kB*T / b_poly * (t + 1.0 / 4.0 / pow((1 - t), 2) - 1.0 / 4.0); }
	else { return 999999.0; }
}

inline void convert_to_vector(float* result, const float mag, const float* direction){
	#pragma unroll
	for (int i = 0; i<DIM; i++){
		result[i] = mag*direction[i];
	}
}

inline float dist(const float* r1, const float* r2){
	float s = 0.0, t;
	#pragma unroll
	for (int j = 0; j<DIM; j++){
		t = r1[j] - r2[j];
		s += t * t;
	}
	return sqrt(s);
}

void forcevector(float* result, float* r1, float* r2, float L);


inline float kfe(float force_mag){
	return ae * exp(force_mag * delxe / kB / T);
}

inline bool ismember(int item, int* array,size_t size){
	bool is_inside = false; 
	for(int k=0; k<size; k++){if(array[k]==item){is_inside = true; break;}}
	return is_inside;
}

float getabsmax(float* arr, size_t sizeofarr);

template <typename t>
void write_to_file(string& fname, t* arr, int rows, int cols){

	ofstream logger;
	std::time_t result = std::time(nullptr);
	logger.open(fname, ios::trunc|ios_base::out);


	logger<<"1D pulling of a 2D gel"<<"\n";
	logger<<"File created at "<<std::asctime(std::localtime(&result));
	logger<<"\n";

	logger<<"Sim dimension : "<<DIM<<"\n";
	logger<<"Simulation time : "<<SIM_TIME<<"\n";
	logger<<"Velocity : "<<vel_x<<"\t"<<vel_y<<"\n";
	logger<<"Maxbound : "<<MAXBOUND<<"\n";

	logger<<"Disorder characteristics : "<<"\n";
	logger<<" -- L_MEAN : "<<L_MEAN<<"\n";
	logger<<" -- L_STD : "<<L_STD<<"\n";
	logger<<"Others : SACBONDS = "<<SACBONDS<<"; IMPLEMENT_PBC = "<<IMPLEMENT_PBC<<"\n";
	
	logger<<"Cracked? : "<<CRACKED<<"\n";

    // logger.open(fname, ios::trunc|ios_base::out);
	for(int i =0; i < rows; i++){
		for(int j = 0; j< cols; j++){
			logger<<arr[i*cols + j]<<"\t";
		}
		logger<<"\n";
	}
    logger.close();
	cout<<"Stored everything in "<<fname<<"!\n";
}


inline float kf(float force){
	return af*expf(force*delxf/kB/T);
}

void __init__(float* L, int* m, float* damage, float* sacdamage, bool* PBC, int n_elems);

#endif