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

<<<<<<< HEAD
#define __params__	
#define DIM 2								// Number of dimensions
#define TIME_STEP 1e-4						// Time step
#define SIM_TIME 8.0						// Simulation time
#define TOL 1e-6							// Tolerance
#define STEPS int(SIM_TIME/TIME_STEP)		// Number of time steps
#define L_MEAN 240.0f						// Average for contour length
#define L_STD 0.001f							// Std. deviation for contour lengths
#define Z_MAX 10							// Max coordination number
#define MAXBOUND 500.0f
#define SACBONDS false
#define IMPLEMENT_PBC true
#define FNAME_STRING "zero_disorder_high_L_"
#define CRACKED true
=======

#define __params__	
#define DIM 2								// Number of dimensions
#define TIME_STEP 1e-2						// Time step
#define SIM_TIME 100.0						// Simulation time
#define TOL 1e-6							// Tolerance
#define STEPS int(SIM_TIME/TIME_STEP)		// Number of time steps
#define L_MEAN 250.0f						// Average for contour length
#define L_STD 75.0f						// Std. deviation for contour lengths
#define MAXBOUND 500.0f
#define SACBONDS false
#define IMPLEMENT_PBC true
#define FNAME_STRING "damage_sbyL_reducedtsby10_"
#define CRACKED false
//#define GENERATOR std::uniform_real_distribution<float>(L_MEAN - L_STD, L_MEAN + L_STD
>>>>>>> BW
#if CRACKED
#define PROB_REMOVAL 0.8
#else
#define PROB_REMOVAL 0.0
#endif


// Define constants
#define __constants__

#define kB 1.38064852e-5					// Boltzmann constant
#define b 0.1								// Persistence length
#define T 300 								// Temperature
#define ae 0.1 								// Strength of bond - includes activation energy
#define delxe 0.15 							// parameter for breaking crosslink connection
#define af 0.3
#define delxf 0.25

using namespace std;

const float PBC_vector[DIM] = {MAXBOUND*1.1, 0};
const float vel[DIM] = {vel_x, vel_y};

class Network {

public:
	Network();
	Network(Network const & source);
	Network(string& fname);
	virtual ~Network();
	Network const & operator=(Network const & other);
	void build_network();
	void apply_crack(Cracklist &);
	void load_network(string&);
	void malloc_network(string&);
	void make_edge_connections(float dely_allowed = 10.0);
	void get_forces(bool);
	void move_top_plate();
	void get_plate_forces(float*, int);
	void optimize(float, float, int);
	float get_weight();
	float set_weight(float);
	bool get_stats();
	int get_current_edges();
	void plotNetwork(int, bool);

protected:

	void clear();
	void copy(Network const & source);
	Gnuplot gnu;
	bool cracked;
	//int DIM;
	int n_nodes;
	int n_elems;
	int n_rside, n_lside, n_bside, n_tside;
	//int num_edges;
	float * R; //n_nodes*DIM
	int * edges; //n_elems*2
	float * forces; //n_nodes*DIM
	float * damage;// n_nodes
	float * L; //n_elems
	bool* PBC; //n_elems
	int* lsideNodes; //max_nodes_on_a_side*2
	int* rsideNodes;
	int* tsideNodes;
	int* bsideNodes;
	bool initialized;
	//add moving nodes to speed up force
	int* moving_nodes;
	int n_moving;

	
};

class sacNetwork : public Network{
	private:
		int* m; //n_elems
		float* sacdamage; // n_elems
	public:
		sacNetwork();
		~sacNetwork();
		sacNetwork(string& fname);
		void clear();
		void malloc_network(string&);
		void load_network(string&);
		void get_forces(bool);


};

/******************************************************************************************************************************************************************************
******************************************************************************************************************************************************************************/
						/**HELPER FUNCTIONS**/

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

/*****Actual Class Functions ******/
Network::Network() {
	initialized = false;

}

Network::Network(Network const & source) {

	initialized = false;
	copy(source);
	initialized = true;

}

Network::Network(string& fname) {

	initialized = false;
	load_network(fname);
	initialized = true;

}

Network::~Network() {

	clear();

}

Network const & Network::operator=(Network const & other) {

	if (this != &other) {
		clear();
		copy(other);
	}

	return *this;
	
}

void Network::clear() {

	free(R);
	free(edges);
	free(forces);
	free(damage);
	free(L);
	free(PBC);
	free(lsideNodes);
	free(rsideNodes);
	free(tsideNodes);
	free(bsideNodes);
	free(moving_nodes);

}

void Network::build_network() {
	//TODO: Add small-network code
	
}


void Network::malloc_network(string& fname){
	
	read_n(n_nodes, n_elems, fname);


	int max_nodes_on_a_side = int(sqrt(n_nodes)*2.0);
	// add memory for side connections
	n_elems += 3*max_nodes_on_a_side;

	size_t sf = sizeof(float);
	size_t si = sizeof(int);
	size_t sb = sizeof(bool);

	R = (float*)malloc(n_nodes*DIM*sf);
	edges = (int*)malloc(n_elems*2*si);
	forces = (float*)malloc(n_nodes*DIM*sf);
	damage = (float*)malloc(n_elems*sf);
	L = (float* )malloc(n_elems*sf);
	PBC = (bool* )malloc(n_elems*sb);
	lsideNodes = (int* )malloc(max_nodes_on_a_side*si);
	rsideNodes = (int* )malloc(max_nodes_on_a_side*si);
	bsideNodes = (int* )malloc(max_nodes_on_a_side*si);
	tsideNodes = (int* )malloc(max_nodes_on_a_side*si);

	// 	initialise n_xside for side nodes
	n_rside = 0;
	n_lside = 0;
	n_bside = 0;
	n_tside = 0;

	// adjust n_elems back to actual
	n_elems -= 3*max_nodes_on_a_side;
}

void Network::load_network(string& fname) {

	if (initialized) {
		clear();
	}
	//Check if file exists
	bool exists = does_file_exist(fname);
	if(!exists){
		cout<<"File does not exist!\n";
		return;
	}
	
	//Malloc all variables
	malloc_network(fname);

	cout<<"Malloc was successful!\n";

	cout<<"Reading the mesh...\n";
	take_input(R, edges, n_nodes, n_elems, fname);
	cout<<"Mesh read successfully!\n";
	cout<<"Number of nodes are: "<<n_nodes<<endl;
	cout<<"Number of elements are: "<<n_elems<<endl;

	side_nodes(R, lsideNodes, rsideNodes, tsideNodes, bsideNodes, \
		n_lside, n_rside, n_tside, n_bside, n_nodes);


	// moving nodes
	n_moving = n_nodes - n_tside - n_bside;
	moving_nodes = (int*)malloc(sizeof(int)*n_moving);
	
	int c = 0;
	for(int i =0; i<n_nodes; i++){
		if(!ismember(i, tsideNodes, n_tside) && !ismember(i, bsideNodes, n_bside)){
			moving_nodes[c] = i;
			c++;
		}
	}
	cout<<"We have "<<c<<" moving nodes\n";


	cout<<"Side nodes written successfully! \n";
	__init__(L, damage, PBC, n_elems);

	if(IMPLEMENT_PBC){
		this->make_edge_connections(15.0);
		cout<<"Number after new connections made: "<<n_elems<<endl;}

}

void Network::copy(Network const & source) {

	cracked = source.cracked;
	n_nodes = source.n_nodes;
	n_elems = source.n_elems;
	n_rside = source.n_rside;
	n_lside = source.n_lside;
	n_bside = source.n_bside;
	n_tside = source.n_tside;
	n_moving = source.n_moving;

	size_t sf = sizeof(float);
	size_t si = sizeof(int);
	size_t sb = sizeof(bool);

	R = (float*)malloc(n_nodes*DIM*sf);
	edges = (int*)malloc(n_elems*2*si);
	forces = (float*)malloc(n_nodes*DIM*sf);
	damage = (float*)malloc(n_elems*sf);
	L = (float* )malloc(n_elems*sf);
	PBC = (bool* )malloc(n_elems*sb);
	lsideNodes = (int* )malloc(n_lside*si);
	rsideNodes = (int* )malloc(n_rside*si);
	bsideNodes = (int* )malloc(n_bside*si);
	tsideNodes = (int* )malloc(n_tside*si);
	moving_nodes = (int *)malloc(n_moving*si);

	for (int i = 0; i < n_nodes * DIM; i++) {
		R[i] = source.R[i];
		forces[i] = source.forces[i];
	}
	
	for (int i = 0; i < n_elems * 2; i++) {
		edges[i] = source.edges[i];
	}
	
	for (int i = 0; i < n_lside; i++) {
		lsideNodes[i] = source.lsideNodes[i];
	}
	for (int i = 0; i < n_rside; i++) {
		rsideNodes[i] = source.rsideNodes[i];
	}
	for (int i = 0; i < n_tside; i++) {
		tsideNodes[i] = source.tsideNodes[i];
	}
	for (int i = 0; i < n_bside; i++) {
		bsideNodes[i] = source.bsideNodes[i];
	}
	for (int i = 0; i < n_moving; i++) {
		moving_nodes[i] = source.moving_nodes[i];
	}

	for (int i = 0; i < n_elems; i++) {
		damage[i] = source.damage[i];
		L[i] = source.L[i];
		PBC[i] = source.PBC[i];
	}

}

void Network::get_forces(bool update_damage = false) {

	int node1, node2;
	int j, k, id; // loop variables
	float r1[DIM]; float r2[DIM] ;
	float edge_force[DIM];
	float rhat[DIM];
	float s;
	float force;


<<<<<<< HEAD
	memset(forces, 0.0, n_nodes*DIM*sizeof(*forces));
=======
	memset(forces, 0.0, n_nodes*DIM*sizeof(float));
>>>>>>> BW

	for (j = 0; j < n_elems; j++){
		// read the two points that form the edge // 2 because 2 points make an edge! Duh.
		node1 = edges[j * 2]; 
		node2 = edges[j * 2 + 1];
		
		// check if pair exists
		if(node1 == -1 || node2 == -1) {
			continue;
		}

		// read the positions
		#pragma unroll
		for(k = 0; k<DIM; k++){
			r1[k] = R[node1*DIM + k]; 
			r2[k] = R[node2*DIM + k];
		}
		
		// check PBC_STATUS
		if (PBC[j]) {
			// add PBC_vector to get new node position
			#pragma unroll
			for (k = 0; k < DIM; k++){
				r2[k] += PBC_vector[k];
			}
			// get force on node1 due to node2
			s = dist(r1, r2);
			unitvector(rhat, r1, r2);
			force = force_wlc(s, L[j]);
			if(force == 999999){edges[j*2] = -1; edges[j*2 +1] = -1; force =0.0;}
			convert_to_vector(edge_force, force, rhat);
			// subtract back the PBC_vector to get original node position
			// #pragma unroll
			// for (k = 0; k < DIM; k++){
			// 	r2[k] -= PBC_vector[k];
			// }
		}
		else{
			s = dist(r1, r2);
			unitvector(rhat, r1, r2);
			force = force_wlc(s, L[j]);
			if(force == 999999){edges[j*2] = -1; edges[j*2 +1] = -1; force =0.0;}
			convert_to_vector(edge_force, force, rhat);
		}
		#pragma unroll
		for (k = 0; k < DIM; k++){
			forces[node1*DIM + k] -= edge_force[k];
			forces[node2*DIM + k] += edge_force[k];
		}
		//update damage if needed
		if (update_damage){
<<<<<<< HEAD
			damage[j] += kfe(force)*TIME_STEP;
			//remove edge ... set to special value
			if(damage[j] > 1.0){
				cout<<"Breaking bond between "
				<<edges[j*2]<<" and "<<edges[2*j +1]<<" F,s = "<<force \
				<<", "<<s<<endl;
=======
			damage[j] = s/L[j]/0.90;
			//remove edge ... set to special value
			if(damage[j] > 1.0){
				cout<<"Breaking bond between "
				<<edges[j*2]<<" and "<<edges[2*j +1]<<" F, s/L, PBC = "<<force \
				<<", "<<s/L[j]<<", "<<PBC[j]<<endl;
>>>>>>> BW
				for(int d = 0; d<DIM; d++){
					cout<<r1[d]<<"\t "<<r2[d]<<"\t";
				}
				cout<<endl;
			edges[j*2] = -1; edges[j*2+1] = -1;}
		}
	}

}

void Network::make_edge_connections(float dely_allowed) {

	std::default_random_engine seed;
<<<<<<< HEAD
	std::normal_distribution<float> generator(L_MEAN, L_STD);
=======
	std::uniform_real_distribution<float> generator(L_MEAN - L_STD, L_MEAN + L_STD);
>>>>>>> BW
	int nl, nr, lnode, rnode;
	for(nl= 0; nl < n_lside; nl++){
		lnode = lsideNodes[nl];
		for(nr= 0; nr < n_rside; nr++){
			rnode = rsideNodes[nr];
			if (fabs(R[lnode*DIM + 1] - R[rnode*DIM + 1]) < dely_allowed){
				edges[n_elems*2] = rnode;
				edges[n_elems*2 + 1] = lnode;
				//cout<<"Connected node "<<lnode<<" and "<<rnode<<"\n";
				L[n_elems] = generator(seed);
				damage[n_elems] = 0.0;
				PBC[n_elems] = true;
				n_elems += 1;
			}
		}
	}

}

void Network::apply_crack(Cracklist & alist) {
	std::default_random_engine seed;
	std::uniform_real_distribution<float> generator(0, 1);
	float equation = 0;
	//make clusters
	vector<int> within_circle;
	vector<int> nodes_to_remove;
	int edges_removed = 0, node1, node2;
	Crack crack;
	float x;
	float y;
	float si,co;
	for(int ci = 0; ci < alist.n_cracks; ci++){
		crack.setter(alist[ci]);
		
		for(int i=0; i<n_nodes; i++){
			x = R[i*DIM];
			y = R[i*DIM+1];
			x -= crack.c[0];
			y -= crack.c[1];
			si = crack.trig[0];
			co = crack.trig[0];
			equation = pow(x*co + y*si, 2)/pow(crack.a[0],2) + \
						pow(x*si - y*co, 2)/pow(crack.a[1],2);

			equation -= 1.0;

			if(equation<=0.0){
				within_circle.push_back(i);
			}
		}
	}

	for(int i= 0; i < n_nodes; i++){
		if(contains(within_circle, i)){
			continue;
		}
		else{
			nodes_to_remove.push_back(i);
		}
	}

	for(int i = 0; i<n_elems; i++){
		node1 = edges[2*i];
		node2 = edges[2*i+ 1];
		if(contains(nodes_to_remove, node1) && contains(nodes_to_remove, node2)){
			if(generator(seed) < PROB_REMOVAL){
				edges[2*i] = -1;
				edges[2*i + 1] = -1;
				edges_removed += 1;
			}
		}
	}
	cout<<"Edges removed : "<<edges_removed<<endl;

}



void Network::plotNetwork(int iter_step, bool first_time){
	ofstream f;
	float c;
	f.open("data.txt",std::ofstream::out | std::ofstream::trunc);
	// std::default_random_engine seed;
	// std::uniform_real_distribution<float> arbitcolor(0.0,1.0);
	int node1, node2;
	if(first_time){
		float r1[DIM];
		float r2[DIM];
		float s;
		int j, k;
		for (j = 0; j < n_elems; j++){
		// read the two points that form the edge // 2 because 2 points make an edge! Duh.
			node1 = edges[j * 2]; 
			node2 = edges[j * 2 + 1];
			
			// check if pair exists
			if(node1 == -1 || node2 == -1) {
				continue;
			}

			// read the positions
			#pragma unroll
			for(k = 0; k<DIM; k++){
				r1[k] = R[node1*DIM + k]; 
				r2[k] = R[node2*DIM + k];
			}
			
			// check PBC_STATUS
			if (PBC[j]) {
				// add PBC_vector to get new node position
				#pragma unroll
				for (k = 0; k < DIM; k++){
					r2[k] += PBC_vector[k];
				}
				// get force on node1 due to node2
				s = dist(r1, r2);
			}
			else{
				s = dist(r1, r2);
			}
			
			c = s/L[j];
			for(int d = 0; d<DIM; d++){
					f<<R[node1*DIM+d]<<"\t";
				}
				f<<c<<endl;
			if(!PBC[j]){
				for(int d = 0; d<DIM; d++){
					f<<R[node2*DIM+d]<<"\t";
				}
			}
			else{
				for(int d = 0; d<DIM; d++){
					f<<(R[node1*DIM+d]+10)<<"\t";
				}				
			}
			f<<c<<endl<<endl;
		}
		f.close();
	}
	else{
	for(int i = 0; i<n_elems; i++){
		c = damage[i];
		node1 = edges[2*i];
		node2 = edges[2*i+1];

		if(node1!=-1 && node2!=-1){
			for(int d = 0; d<DIM; d++){
					f<<R[node1*DIM+d]<<"\t";
				}
				f<<c<<endl;
			if(!PBC[i]){
				for(int d = 0; d<DIM; d++){
					f<<R[node2*DIM+d]<<"\t";
				}
			}
			else{
				for(int d = 0; d<DIM; d++){
					f<<(R[node1*DIM+d]+10)<<"\t";
				}				
			}
		f<<c<<endl<<endl;
		}
	}
	f.close();
	}
	//Plot to h
	stringstream convert;
	convert<<"[-100:"<<MAXBOUND+100<<"]";
	string xrange = convert.str();
	gnu.cmd("set xrange " + xrange);
	gnu.cmd("load 'viridis.pal'");
	gnu.cmd("set key off");
	gnu.cmd("set colorbox default vertical");
	gnu.cmd("set cbrange [0:1]");
	gnu.cmd("set cbtics ('0.0' 0.0,'1.0' 1.0) offset 0,0.5 scale 0");
	gnu.cmd("plot 'data.txt' every ::0::1 u 1:2:3 w l lc palette lw 2 title '"+\
		std::to_string(iter_step)+"'");
	gnu.cmd("set term png size 1200,900");
	gnu.cmd("set output '"+std::to_string(iter_step)+".png'");
	gnu.cmd("replot");
	gnu.cmd("set term x11");
	gnu.reset_plot();
}

int Network::get_current_edges(){
	int n_elems_current = 0;
	for(int i=0; i<n_elems; i++){
		if(edges[i*2] != -1 && edges[i*2 + 1]!=-1){
			n_elems_current += 1;
		}
	}
	cout<<"Number of edges are "<<n_elems_current<<"\n";
	return n_elems_current;
}

bool Network::get_stats(){
	int node1, node2;
	int j, k, id; // loop variables
	float r1[DIM]; float r2[DIM] ;
	float mean_x = 0.0;
	float mean_t = 0.0;
	float var_x = 0.0;
	float var_t = 0.0;
	float max_x = 0.0;
	float max_t = 0.0;
	float s;
	int short_L = 0;
	float shortage = 0.0;
	int c = 0;

	for (j = 0; j < n_elems; j++){
		// read the two points that form the edge // 2 because 2 points make an edge! Duh.
		node1 = edges[j * 2]; 
		node2 = edges[j * 2 + 1];
		
		// check if pair exists
		if(node1 == -1 || node2 == -1) {
			continue;
		}
		c++;


		// read the positions
		#pragma unroll
		for(k = 0; k<DIM; k++){
			r1[k] = R[node1*DIM + k]; 
			r2[k] = R[node2*DIM + k];
		}
		
		// check PBC_STATUS
		if (PBC[j]) {
			// add PBC_vector to get new node position
			#pragma unroll
			for (k = 0; k < DIM; k++){
				r2[k] += PBC_vector[k];
			}
			// check if L too short
			s = dist(r1, r2);
			mean_x += s;
			mean_t += s/L[j];
			if(s>max_x){max_x =s;}
			if(s/L[j]>max_t){max_t=s/L[j];}
			var_x += s*s;
			var_t += s*s/L[j]/L[j];
			if(s >= 0.99*L[j]){
				short_L++;
				shortage = s - L[j];
			}
		}
		else{
			s = dist(r1, r2);
			mean_x += s;
			mean_t += s/L[j];
			if(s>max_x){max_x =s;}
			if(s/L[j]>max_t){max_t=s/L[j];}
			var_x += s*s;
			var_t += s*s/L[j]/L[j];
			if(s >= 0.99*L[j]){
				short_L++;
				shortage = s - L[j];
			}
		}
	}
	mean_x /= c;
	var_x /= c;
	var_x -= mean_x*mean_x;

	mean_t /= c;
	var_t /= c;
	var_t -= mean_t*mean_t;
	cout<<"Number of edges: "<<c<<endl;
	cout<<"Number of edges with x>0.99L: "<<short_L<<endl;
	cout<<"node-node distance mean: "<<mean_x<<endl;
	cout<<"node-node distance std_dev: "<<sqrt(var_x)<<endl;
	cout<<"node-node distance max: "<<max_x<<endl;
	cout<<"\n";
	cout<<"node-node x/L mean: "<<mean_t<<endl;
	cout<<"node-node x/L std_dev: "<<sqrt(var_t)<<endl;	
	cout<<"node-node x/L max: "<<max_t<<endl;

<<<<<<< HEAD
	if(shortage/get_weight() > 0.1){
		cout<<"Disorder too high for given mesh! Exiting...\n";
		return true;
	}
=======
	// if(shortage/get_weight() > 0.1){
	// 	cout<<"Disorder too high for given mesh! Exiting...\n";
	// 	return true;
	// }
>>>>>>> BW
	if((float)c/(float)n_elems < 0.02){
		cout<<"Too few edges remain in given mesh! Exiting...\n";
		return true;
	}
	return false;
}

float Network::get_weight(){
	float sum_length = 0.0;
	for(int i = 0; i < n_elems; i++){
		if(edges[i*2] != -1 && edges[i*2 + 1]!=-1){
			sum_length += L[i];
		}
	}
	cout<<"\nNetwork weight is "<<sum_length<<"\n";
	return sum_length;
}

float Network::set_weight(float weight){
	float curr_weight = get_weight();
	float alpha = weight/curr_weight;
	for(int i = 0; i < n_elems; i++){
		if(edges[i*2] != -1 && edges[i*2 + 1]!=-1){
			L[i] *= alpha;
		}
	}
	cout<<"\n Network weight reset to "<<weight<<"\n";
	cout<<"Average L now is : "<<weight/n_elems<<"\n";
	return alpha;
}

void Network::move_top_plate(){
	int node;
	for(int i = 0; i<n_tside; i++){
		node = tsideNodes[i];
		#pragma unroll
		for(int d=0; d<DIM; d++){
			R[node*DIM + d] += TIME_STEP*vel[d]; 
		}
	}
}

float getabsmax(float* arr, size_t sizeofarr){
	float max_elem = 0.0;
	for(int i = 0; i<sizeofarr; i++){
		if(fabs(arr[i])>max_elem){max_elem = fabs(arr[i]);}
	}
	return max_elem;
}

void Network::optimize(float eta = 0.1, float alpha = 0.9, int max_iter = 800){
	float* rms_history = new float[n_moving*DIM](); // () allows 0.0 initialization
<<<<<<< HEAD
	float* delR = new float[n_moving*DIM]();
	float g;
=======
	float g, delR;
>>>>>>> BW
	char p;
	int id, d, node;
	for(int step = 0; step < max_iter; step++){
		get_forces(false);
<<<<<<< HEAD
=======
		// plotNetwork(step, true);
		// cin>>p;
>>>>>>> BW
		if(getabsmax(forces,n_nodes*DIM)>TOL){
			for(id = 0; id < n_moving; id++){
				node = moving_nodes[id];
				#pragma unroll
				for(d = 0; d<DIM; d++){
					g = forces[DIM*node+d];
<<<<<<< HEAD
					rms_history[id*DIM + d] = alpha*rms_history[id] + (1-alpha)*g*g;
					delR[id*DIM + d] = sqrt(1.0/(rms_history[id] + TOL))*eta*g;
					R[node*DIM + d] += delR[id*DIM + d];
=======
					rms_history[id*DIM + d] = alpha*rms_history[id*DIM + d] + (1-alpha)*g*g;
					delR = sqrt(1.0/(rms_history[id*DIM + d] + TOL))*eta*g;
					R[node*DIM + d] += delR;
>>>>>>> BW
				}		
			}
		}
		else{
			break;
		}
	}
	//
	get_forces(true);
	delete[] rms_history;
<<<<<<< HEAD
	delete[] delR;
=======
>>>>>>> BW
}

void Network::get_plate_forces(float* plate_forces, int iter){
	int node;
	for(int i=0; i<n_tside; i++){
		node = tsideNodes[i];
		plate_forces[iter*DIM+0] += forces[node*DIM + 0];
		plate_forces[iter*DIM+1] += forces[node*DIM + 1];
		if(DIM>2){
			plate_forces[iter*DIM+2] = forces[node*DIM + 2];
		}
	}
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

sacNetwork::sacNetwork(){
	initialized = false;
}

sacNetwork::sacNetwork(string& fname){
	initialized = false;
	load_network(fname);
	initialized = true;
}

sacNetwork::~sacNetwork(){
	clear();
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

void sacNetwork::clear() {
	// DO NOT CLEAR baseclass' (Network) variables from here
	// free(R);
	// free(edges);
	// free(forces);
	// free(damage);
	free(sacdamage);
	free(m);
	// free(L);
	// free(PBC);
	// free(lsideNodes);
	// free(rsideNodes);
	// free(tsideNodes);
	// free(bsideNodes);
	// free(moving_nodes);

}

void sacNetwork::get_forces(bool update_damage = false) {

	int node1, node2;
	int j, k, id; // loop variables
	float r1[DIM]; float r2[DIM] ;
	float edge_force[DIM];
	float rhat[DIM];
	float s;
	float force;


	memset(forces, 0.0, n_nodes*DIM*sizeof(*forces));

	for (j = 0; j < n_elems; j++){
		// read the two points that form the edge // 2 because 2 points make an edge! Duh.
		node1 = edges[j * 2]; 
		node2 = edges[j * 2 + 1];
		
		// check if pair exists
		if(node1 == -1 || node2 == -1) {
			continue;
		}

		// read the positions
		#pragma unroll
		for(k = 0; k<DIM; k++){
			r1[k] = R[node1*DIM + k]; 
			r2[k] = R[node2*DIM + k];
		}
		
		// check PBC_STATUS
		if (PBC[j]) {
			// add PBC_vector to get new node position
			#pragma unroll
			for (k = 0; k < DIM; k++){
				r2[k] += PBC_vector[k];
			}
			// get force on node1 due to node2
			s = dist(r1, r2);
			unitvector(rhat, r1, r2);
			force = force_wlc(s, L[j]);
			if(force == 999999){edges[j*2] = -1; edges[j*2 +1] = -1; 
				 force =0.0;}
			convert_to_vector(edge_force, force, rhat);
			// subtract back the PBC_vector to get original node position
			// #pragma unroll
			// for (k = 0; k < DIM; k++){
			// 	r2[k] -= PBC_vector[k];
			// }
		}
		else{
			s = dist(r1, r2);
			unitvector(rhat, r1, r2);
			force = force_wlc(s, L[j]);
			if(force == 999999){edges[j*2] = -1; edges[j*2 +1] = -1; 
				force =0.0;}
			convert_to_vector(edge_force, force, rhat);
		}
		#pragma unroll
		for (k = 0; k < DIM; k++){
			forces[node1*DIM + k] -= edge_force[k];
			forces[node2*DIM + k] += edge_force[k];
		}


		//update damage if needed
		if (update_damage){
			damage[j] += kfe(force)*TIME_STEP;
			if(m[j] > 0){
				sacdamage[j] += kf(force)*TIME_STEP;
				if(sacdamage[j] > 1.0){
					L[j] += (L_MEAN - L[j])/m[j] ; 
					m[j] -= 1;
					sacdamage[j] = 0.0;
				}
			}
			
			//remove edge ... set to special value
			if(damage[j] > 1.0){
				cout<<"Breaking bond between "
				<<edges[j*2]<<" and "<<edges[2*j +1]<<" F,s = "<<force \
				<<", "<<s<<endl;
				for(int d = 0; d<DIM; d++){
					cout<<r1[d]<<"\t "<<r2[d]<<"\t";
				}
				cout<<endl;
			edges[j*2] = -1; edges[j*2+1] = -1;}
		}
	}

}

void sacNetwork::malloc_network(string& fname){
	
	read_n(n_nodes, n_elems, fname);


	int max_nodes_on_a_side = int(sqrt(n_nodes)*2.0);
	// add memory for side connections
	n_elems += 3*max_nodes_on_a_side;

	size_t sf = sizeof(float);
	size_t si = sizeof(int);
	size_t sb = sizeof(bool);

	R = (float*)malloc(n_nodes*DIM*sf);
	edges = (int*)malloc(n_elems*2*si);
	forces = (float*)malloc(n_nodes*DIM*sf);
	damage = (float*)malloc(n_elems*sf);
	sacdamage = (float*)malloc(n_elems*sf);
	m = (int*)malloc(n_elems*si);
	L = (float* )malloc(n_elems*sf);
	PBC = (bool* )malloc(n_elems*sb);
	lsideNodes = (int* )malloc(max_nodes_on_a_side*si);
	rsideNodes = (int* )malloc(max_nodes_on_a_side*si);
	bsideNodes = (int* )malloc(max_nodes_on_a_side*si);
	tsideNodes = (int* )malloc(max_nodes_on_a_side*si);


	// 	initialise n_xside for side nodes
	n_rside = 0;
	n_lside = 0;
	n_bside = 0;
	n_tside = 0;

	// adjust n_elems back to actual
	n_elems -= 3*max_nodes_on_a_side;
}

void sacNetwork::load_network(string& fname) {

	if (initialized) {
		clear();
	}
	//Check if file exists
	bool exists = does_file_exist(fname);
	if(!exists){
		cout<<"File does not exist!\n";
		return;
	}
	
	//Malloc all variables
	malloc_network(fname);

	cout<<"Malloc was successful!\n";

	cout<<"Reading the mesh...\n";
	take_input(R, edges, n_nodes, n_elems, fname);
	cout<<"Mesh read successfully!\n";
	cout<<"Number of nodes are: "<<n_nodes<<endl;
	cout<<"Number of elements are: "<<n_elems<<endl;

	side_nodes(R, lsideNodes, rsideNodes, tsideNodes, bsideNodes, \
		n_lside, n_rside, n_tside, n_bside, n_nodes);


	// moving nodes
	n_moving = n_nodes - n_tside - n_bside;
	moving_nodes = (int*)malloc(sizeof(int)*n_moving);
	
	int c = 0;
	for(int i =0; i<n_nodes; i++){
		if(!ismember(i, tsideNodes, n_tside) && !ismember(i, bsideNodes, n_bside)){
			moving_nodes[c] = i;
			c++;
		}
	}
	cout<<"We have "<<c<<" moving nodes\n";


	cout<<"Side nodes written successfully! \n";
	__init__(L, m, damage, sacdamage, PBC, n_elems);

	if(IMPLEMENT_PBC){	
		this->make_edge_connections(15.0);
		cout<<"Number after new connections made: "<<n_elems<<endl;
	}
}

int main() {

	//string path = "/media/konik/Research/2D sacrificial bonds polymers/cpp11_code_with_cuda/template2d.msh";
	string path = "./template2d_z4.msh";
	#if SACBONDS
	#define DECL_NET sacNetwork test_network(path)
	#else
	#define DECL_NET Network test_network(path)
	#endif
	DECL_NET;

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
<<<<<<< HEAD
		if((i+1)%10 == 0){
=======
		if((i+1)%100 == 0){
>>>>>>> BW
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
	return 0;
}
