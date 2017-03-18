#include <iostream>
#include <cmath>
#include "Crack.cpp"
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
#include "gnuplot_i.hpp"
using namespace std;

void side_nodes(float* R,\
	int* lnodes, int* rnodes, int* tnodes, int* bnodes,\
	int& n_lside, int& n_rside, int& n_tside, int& n_bside, int n){
	// make the list of side nodes
	for(int i=0; i<n; i++){
		if(fabs(R[i*DIM]-0.0)<TOL){lnodes[n_lside]=i; n_lside++;}
		if(fabs(R[i*DIM]-MAXBOUND)<TOL){rnodes[n_rside]=i; n_rside++;}
		if(fabs(R[i*DIM + 1]-0.0)<TOL){bnodes[n_bside]=i; n_bside++;}
		if(fabs(R[i*DIM + 1]-MAXBOUND)<TOL){tnodes[n_tside]=i; n_tside++;}
	}
	//cout << __LINE__ << endl;
}

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

void mapping(int& edge_counter, int elem_type, int n_vertices, stringstream& input, int* edges){
	int* local_nodes = new int[n_vertices];
	int garbage[3];
	
	input>>garbage[0]>>garbage[1]>>garbage[2];
	for(int n = 0; n<n_vertices; n++){
		input>>local_nodes[n];
	}
	switch(elem_type){
		case 2: 
		// [0,1], [0,2], [1, 2]
		edges[edge_counter*2] = local_nodes[0] -1 ;
		edges[edge_counter*2+1] = local_nodes[1] - 1 ;
		edges[edge_counter*2+2] = local_nodes[0] - 1 ;
		edges[edge_counter*2+3] = local_nodes[2] - 1 ;
		edges[edge_counter*2+4] = local_nodes[1] - 1 ;
		edges[edge_counter*2+5] = local_nodes[2] - 1 ;
		edge_counter += 3;
		break;
		case 3:
		// [0,1],[0,3],[1,2],[2,3]
		edges[edge_counter*2] = local_nodes[0] -1 ;
		edges[edge_counter*2+1] = local_nodes[1] - 1 ;
		edges[edge_counter*2+2] = local_nodes[0] - 1 ;
		edges[edge_counter*2+3] = local_nodes[3] - 1 ;
		edges[edge_counter*2+4] = local_nodes[1] - 1 ;
		edges[edge_counter*2+5] = local_nodes[2] - 1 ;
		edges[edge_counter*2+6] = local_nodes[2] - 1 ;
		edges[edge_counter*2+7] = local_nodes[3] - 1 ;
		edge_counter += 4;
		break;	 
		//TODO: write mapping code
	}
	delete[] local_nodes;
}

void mapping(int& edge_counter, int elem_type){
	switch(elem_type){
		case 2:
			edge_counter += 3;
			break;
		case 3:
			edge_counter += 4;
			break;
		//TODO: add other mappings
	}

}

void read_n(int& n_nodes, int& n_elems, string& fname){

	string line;
	ifstream source;
	stringstream in_pos;
	source.open(fname);
	bool read_num_nodes = false, read_n_elems = false;
	
	n_nodes = 0;
	n_elems = 0;
	do{
		getline(source, line);
		if(line.find("$Nodes") != string::npos){
			getline(source, line);
			n_nodes = stoi(line);
			read_num_nodes = true;
		}
	}while(!read_num_nodes);
	do{
		getline(source, line);
		if(line.find("$Elements") != string::npos){
			int id, elem_type;
			int num_vertices, num_lines;
			getline(source, line);
			num_lines = stoi(line);
			int c = 0;
			for(int i=0; i<num_lines; i++){
				getline(source, line);
				//cout<<line<<endl;
				in_pos<<line;
				in_pos>>id;
				in_pos>>elem_type;
				//cout<<"element type is "<<elem_type<<endl;
				num_vertices = get_num_vertices(elem_type);
				//if (num_vertices > 3){cout<<"Found elem_type 3 at line "<<i+1<<endl;}
				if(num_vertices>0){
					mapping(n_elems, elem_type);
				}
				in_pos.str(std::string());
				in_pos.clear();
			}
			read_n_elems = true;
			}
			//TODO: get elements
	}while(!read_n_elems);
	source.close();
}

void take_input(float* R, int* edges, int n_nodes, int n_elems, string& fname) {
	
	string line;
	ifstream source;
	stringstream in_pos;
	source.open(fname);

	bool read_num_nodes = false, read_n_elems = false;

	do{
		getline(source, line);
		//cout<<line<<"\n";
		//cout<<read_n_elems;
		if(line.find("$Nodes") != string::npos){
			getline(source, line);
			float r[3]; int id;
			for(int i=0; i<n_nodes; i++){
				getline(source, line);
				in_pos<<line;
				//cout<<in_pos.str()<<endl;
				in_pos>>id>>r[0]>>r[1]>>r[2];
				for(int d=0; d<DIM; d++){
					R[(id-1)*DIM + d] = r[d];
				}
				in_pos.str(std::string());
				in_pos.clear();
			}
			read_num_nodes = true;
		}

	}while(!read_num_nodes);
	// for(int i=0; i<num_nodes; i++){
	// 	cout<<endl;
	// 	for(int d =0; d<DIM; d++){
	// 		cout<<R[i*DIM + d]<<"\t";
	// 	}
	// }

	do{
		getline(source, line);
		//cout<<line<<"\n";
		if(line.find("$Elements") != string::npos){
			int id, elem_type;
			int num_vertices, num_lines;
			getline(source, line);
			num_lines = stoi(line);
			int c = 0;
			for(int i=0; i<num_lines; i++){
				getline(source, line);
				//cout<<line<<endl;
				in_pos<<line;
				in_pos>>id;
				in_pos>>elem_type;
				//cout<<"element type is "<<elem_type<<endl;
				num_vertices = get_num_vertices(elem_type);
				if(num_vertices>0){
					mapping(c, elem_type, num_vertices, in_pos, edges);
				}
				in_pos.str(std::string());
				in_pos.clear();
			}
			read_n_elems = true;
			if(c==n_elems){
				cout<<"Everything's alright here!\n";
			}
			else{
				cout<<"Damn! Something's not right here...\n";
			}
			}
			//TODO: get elements
	}while(!read_n_elems);
	source.close();
}

void __init__(float* L, float* damage, bool* PBC, int n_elems){
	std::default_random_engine seed;
	std::uniform_real_distribution<float> generator(L_MEAN - L_STD, L_MEAN + L_STD);
	// std::normal_distribution<float> generator(L_MEAN, L_STD);
	#pragma unroll
	for(int i=0; i<n_elems; i++){
		L[i] = generator(seed);
		damage[i] = 0.0;
		PBC[i] = false;
	}	
}

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
	if (t < 0.99){ return kB*T / b * (t + 1.0 / 4.0 / pow((1 - t), 2) - 1.0 / 4.0); }
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

void forcevector(float* result, float* r1, float* r2, float L){
	float rhat[DIM];
	float s = dist(r1, r2);
	unitvector(rhat, r1, r2);
	float force = force_wlc(s, L);
	convert_to_vector(result, force, rhat);
}


inline float kfe(float force_mag){
	return ae * exp(force_mag * delxe / kB / T);
}

inline bool ismember(int item, int* array,size_t size){
	bool is_inside = false; 
	for(int k=0; k<size; k++){if(array[k]==item){is_inside = true; break;}}
	return is_inside;
}

float getabsmax(float* arr, size_t sizeofarr){
	float max_elem = 0.0;
	for(int i = 0; i<sizeofarr; i++){
		if(fabs(arr[i])>max_elem){max_elem = fabs(arr[i]);}
	}
	return max_elem;
}

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

void __init__(float* L, int* m, float* damage, float* sacdamage, bool* PBC, int n_elems){
	std::default_random_engine seed;
	// std::normal_distribution<float> generator(L_MEAN, L_STD);
	std::uniform_real_distribution<float> generator(L_MEAN - L_STD, L_MEAN + L_STD);
	std::uniform_real_distribution<float> hidden(0.1*L_MEAN, 0.2*L_MEAN);
	std::uniform_int_distribution<int> number(2, 4);
	#pragma unroll
	for(int i=0; i<n_elems; i++){
		L[i] = generator(seed);
		m[i] = number(seed);
		L[i] -= m[i]*hidden(seed); 
		sacdamage[i] = 0.0;
		damage[i] = 0.0;
		PBC[i] = false;
	}	
}