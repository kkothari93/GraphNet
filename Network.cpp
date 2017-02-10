// #define __params__	
// #define DIM 2								// Number of dimensions
// #define TIME_STEP 1e-4						// Time step
// #define SIM_TIME 10.0						// Simulation time
// #define TOL 1e-6							// Tolerance
// #define STEPS int(SIM_TIME/TIME_STEP)		// Number of time steps
// #define L_MEAN 120.0f						// Average for contour length
// #define L_STD 4.0f							// Std. deviation for contour lengths
// #define Z_MAX 10							// Max coordination number
// #define MAXBOUND 500.0f

// // Define constants
// #define __constants__

// #define kB 1.38064852e-5					// Boltzmann constant
// #define b 0.1								// Persistence length
// #define T 300 								// Temperature
// #define ae 0.1 								// Strength of bond - includes activation energy
// #define delxe 0.15 							// parameter for breaking crosslink connection
// #define BLOCK_SIZE 1024

#include "Network.h"
using namespace std;


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
		//case 3: return 4;
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
		//TODO: write mapping code

	}
	delete[] local_nodes;
}

void mapping(int& edge_counter, int elem_type){
	switch(elem_type){
		case 2:
			edge_counter += 3;
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
	std::normal_distribution<float> generator(L_MEAN, L_STD);
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
	float toRet =  ae * exp(force_mag * delxe / kB / T);
	return toRet;
}

inline bool ismember(int item, int* array,size_t size){
	bool is_inside = false; 
	for(int k=0; k<size; k++){if(array[k]==item){is_inside = true; break;}}
	return is_inside;
}

/*****Actual Class Functions ******/
Network::Network() {

	/**
	DIM = 2;
	cracked = false;
	n_nodes = 1600;
	n_elems = 11000;
	n_rside = 0;
	n_lside = 0;
	n_bside = 0;
	n_tside = 0;
	R = new float[n_nodes * DIM];
	forces = new float[n_nodes * DIM];
	edges = new int[Z_MAX * n_nodes * 2];
	int max_nodes_on_a_side = int(sqrt(n_nodes))*2;
	lsideNodes = new int[max_nodes_on_a_side];
	rsideNodes = new int[max_nodes_on_a_side];
	tsideNodes = new int[max_nodes_on_a_side];
	bsideNodes = new int[max_nodes_on_a_side];
	damage = new float[2 * n_elems];
	L = new float[2 * n_elems];
	PBC = new bool[2 * n_elems];
	edge_matrix = new bool*[n_elems];
	for (int i = 0; i < n_elems; i++) {
		edge_matrix[i] = new bool[i+1];
		for (int j = 0; j <= i; i++) {
			edge_matrix[i][j] = false;
		}
	}
	**/
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

	// for (int i = 0; i < n_nodes; i++) {
	// 	free(edge_matrix[i]);
	// }
	free(edge_matrix);
	free(L);
	free(PBC);
	free(lsideNodes);
	free(rsideNodes);
	free(tsideNodes);
	free(bsideNodes);
	free(moving_nodes);

}


// void Network::clear() {

// 	delete[] R;
// 	R = NULL;
// 	delete[] edges;
// 	edges = NULL;
// 	delete[] forces;
// 	forces = NULL;
// 	delete[] damage;
// 	damage = NULL;

// 	for (int i = 0; i < n_nodes; i++) {
// 		delete[] edge_matrix[i];
// 		edge_matrix[i] = NULL;
// 	}
// 	delete[] edge_matrix;
// 	edge_matrix = NULL;
// 	delete[] L;
// 	L = NULL;
// 	delete[] PBC;
// 	PBC = NULL;
// 	delete[] lsideNodes;
// 	lsideNodes = NULL;
// 	delete[] rsideNodes;
// 	rsideNodes = NULL;
// 	delete[] tsideNodes;
// 	tsideNodes = NULL;
// 	delete[] bsideNodes;
// 	bsideNodes = NULL;

// }

void Network::build_network() {

	
	//num_edges = Z_MAX * n_nodes * 2;
	/** TO-DO **/
	
}

// TODO: Check what's wrong here
// template <typename t>
// void malloc_2d(t** array_2d, int r, int c){
// 	array_2d = (t**)malloc(r*c*sizeof(t*));
// 	for(int i=0; i< r; i++){
// 		array_2d[i] = (t*)malloc(c*sizeof(t));
// 	}
// }

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

	// malloc_2d<bool>(edge_matrix, n_nodes, n_nodes);
	edge_matrix = (bool**)malloc(n_nodes*n_nodes*sizeof(bool*));
	for(int i = 0; i<n_nodes; i++){
		edge_matrix[i] = (bool*)malloc((i+1)*sizeof(bool));
	}

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

	this->make_edge_connections(15.0);
	cout<<"Number after new connections made: "<<n_elems<<endl;

	
	for (int i = 0; i < n_elems; i++) {
		int node1 = edges[2*i];
		int node2 = edges[2*i + 1];
		if(node1 == -1 || node2 == -1){
			continue;
		}
		else if (node1 < node2) {
			edge_matrix[node1][node2] = true;
		}

	}

}

void Network::copy(Network const & source) {

	cracked = source.cracked;
	//DIM = source.DIM;
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

	// malloc_2d<bool>(edge_matrix, n_nodes, n_nodes);
	edge_matrix = (bool**)malloc(n_nodes*sizeof(bool*));
	for(int i = 0; i<n_nodes; i++){
		edge_matrix[i] = (bool*)malloc((i+1)*sizeof(bool));
	}

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

	for (int i = 0; i < n_elems; i++) {
		int node1 = edges[2*i];
		int node2 = edges[2*i + 1];
		if(node1 == -1 || node2 == -1){
			continue;
		}
		else if (node1 < node2) {
			edge_matrix[node1][node2] = true;
		}

	}

}

bool Network::get_forces(bool update_damage = false, int lo, int hi) {

	bool BROKEN = false;
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
		if ((node1 <= lo || node1 >= hi || node2 <= lo || node2 >= hi)) {
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
		if (update_damage && (node1 < node2)){
			damage[j] += kfe(force)*TIME_STEP;
			//remove edge ... set to special value
			if(damage[j] > 1.0){
				cout<<"Breaking bond between "
				<<edges[j*2]<<" and "<<edges[2*j +1]<<" F,s = "<<force \
				<<", "<<s<<endl;
				for(int d = 0; d<DIM; d++){
					cout<<r1[d]<<"\t "<<r2[d]<<"\t";
				}
				cout<<endl;
				edges[j*2] = -1; edges[j*2+1] = -1;
				BROKEN = true;
			}

		}
	}
	return BROKEN;
}

void Network::make_edge_connections(float dely_allowed) {

	std::default_random_engine seed;
	std::normal_distribution<float> generator(L_MEAN, L_STD);
	int nl, nr, lnode, rnode;
	for(nl= 0; nl < n_lside; nl++){
		lnode = lsideNodes[nl];
		for(nr= 0; nr < n_rside; nr++){
			rnode = rsideNodes[nr];
			if (fabs(R[lnode*DIM + 1] - R[rnode*DIM + 1]) < dely_allowed){
				edges[n_elems*2] = rnode;
				edges[n_elems*2 + 1] = lnode;
				cout<<"Connected node "<<lnode<<" and "<<rnode<<"\n";
				L[n_elems] = generator(seed);
				damage[n_elems] = 0.0;
				PBC[n_elems] = true;
				n_elems += 1;
			}
		}
	}

}

void Network::apply_crack(Crack const & crack) {

	float equation = 0;
	vector<int> nodes_to_remove;
	for(int i=0; i<n_nodes; i++){
		equation = 0;
		for(int d=0; d<DIM; d++){
			equation += pow(R[i*DIM+d]-crack.c[d],2)/pow(crack.a[d],2);
		}
		equation -= 1.0;

		if(equation<=0.0){
			nodes_to_remove.push_back(i);
		}
	}
	int edges_removed = 0, node1, node2;
	for(int i = 0; i<n_elems; i++){
		node1 = edges[2*i];
		node2 = edges[2*i+ 1];
		if(contains(nodes_to_remove, node1) || contains(nodes_to_remove, node2)){
			edges[2*i] = -1;
			edges[2*i + 1] = -1;
			edges_removed += 1;
		}
	}
	cout<<"Edges removed : "<<edges_removed<<endl;

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

	cout<<"Number of edges with x>0.99L: "<<short_L<<endl;
	cout<<"node-node distance mean: "<<mean_x<<endl;
	cout<<"node-node distance std_dev: "<<sqrt(var_x)<<endl;
	cout<<"node-node distance max: "<<max_x<<endl;
	cout<<"\n";
	cout<<"node-node x/L mean: "<<mean_t<<endl;
	cout<<"node-node x/L std_dev: "<<sqrt(var_t)<<endl;	
	cout<<"node-node x/L max: "<<max_t<<endl;

	if(shortage/get_weight() > 0.1){
		cout<<"Disorder too high for given mesh! Exiting...\n";
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

void Network::set_weight(float weight){
// TODO
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

void Network::optimize(float eta = 0.1, float alpha = 0.9, int max_iter = 800, int lo, int hi, bool& BROKEN){
	float* rms_history = new float[n_moving*DIM](); // () allows 0.0 initialization
	float* delR = new float[n_moving*DIM]();
	float g;
	int id, d, node;
	for(int step = 0; step < max_iter; step++){
		get_forces(false, lo, hi);
		if(getabsmax(forces,n_nodes*DIM)>TOL){
			for(id = 0; id < n_moving; id++){
				node = moving_nodes[id];
				#pragma unroll
				for(d = 0; d<DIM; d++){
					g = forces[DIM*node+d];
					rms_history[id*DIM + d] = alpha*rms_history[id] + (1-alpha)*g*g;
					delR[id*DIM + d] = sqrt(1.0/(rms_history[id] + TOL))*eta*g;
					R[node*DIM + d] += delR[id*DIM + d];
				}		
			}
		}
		else{
			break;
		}
	}
	BROKEN = get_forces(true);
	delete[] rms_history;
	delete[] delR;
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

// void Network::split_for_MPI(float * R_split, int * edges_split, float * forces, int number_of_procs, int curr_proc_rank) {

// 	int R_total = n_nodes * DIM;
// 	int R_split_size = ceil(R_total/number_of_procs);

// 	if (curr_proc_rank < number_of_procs-1) {

// 		for (int i = curr_proc_rank*R_split_size, int j = 0; i < (curr_proc_rank+1)*(R_split_size); i++, j++) {
// 			R_split[j] = R[i];
// 		}

// 	}
// 	else if (curr_proc_rank == number_of_procs-1) {

// 		for (int i = 0; )
// 	}


// }
template <typename t>
void write_to_file(string& fname, t* arr, int rows, int cols){

	ofstream logger;
	std::time_t result = std::time(nullptr);


	cout<<"1D pulling of a 2D gel"<<endl;
	cout<<"File created at "<<std::asctime(std::localtime(&result));
	cout<<endl;

	cout<<"Sim dimension : "<<DIM<<endl;
	cout<<"Simulation time : "<<SIM_TIME<<endl;
	cout<<"Velocity : "<<vel_x<<"\t"<<vel_y<<endl;
	cout<<endl;

	cout<<"Disorder characteristics : "<<endl;
	cout<<" -- L_MEAN : "<<L_MEAN<<endl;
	cout<<" -- L_STD : "<<L_STD<<endl;
	cout<<endl;
	
	cout<<"Cracked? : "<<CRACKED<<endl;

	logger.open(fname, ios::trunc|ios_base::out);
	for(int i =0; i < rows; i++){
		for(int j = 0; j< cols; j++){
			logger<<arr[i*cols + j]<<"\t";
		}
		logger<<"\n";
	}
}

int main() {

	//string path = "/media/konik/Research/2D sacrificial bonds polymers/cpp11_code_with_cuda/template2d.msh";
	string path = "./template2d.msh";
	Network test_network(path);
	float* plate_forces;

	test_network.get_weight();
	bool should_stop = test_network.get_stats();

	if(should_stop){return 0;}

	plate_forces = (float*)malloc(sizeof(float)*DIM*STEPS);
	memset(plate_forces, 0.0, STEPS*DIM*sizeof(*plate_forces));

	clock_t t = clock();
	cout<<"\n Will run for "<<STEPS<<":\n";
	for(int i = 0; i<STEPS; i++ ){
		test_network.optimize();
		test_network.move_top_plate();
		test_network.get_plate_forces(plate_forces, i);
		if((i+1)%100 == 0){
			test_network.get_stats();
			test_network.get_current_edges();
			cout<<"Step "<<(i+1)<<" took "<<float(clock()-t)/CLOCKS_PER_SEC<<" s\n";
			t = clock();  // reset clock
		}
	}

	string fname = "forces_disorder_0.1.txt";
	write_to_file<float>(fname, plate_forces, STEPS, DIM);

	free(plate_forces);
	return 0;
}





