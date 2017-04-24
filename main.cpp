#include <iostream>
#include <cmath>
#include <cstdlib>
#include <random>
#include <fstream>
#include <cmath>
#include <time.h>
#include <vector>
#include <algorithm>
#include <ctime>
#include <chrono>


#include "gnuplot_i.hpp"
#include "input.h"


#define DIM 2								// Number of dimensions
#define TIME_STEP 1e-4						// Time step
#define SIM_TIME 10.0						// Simulation time
#define TOL 1e-6							// Tolerance
#define STEPS int(SIM_TIME/TIME_STEP)		// Number of time steps
#define L_MEAN 120.0f						// Average for contour length
#define L_STD 4.0f							// Std. deviation for contour lengths
#define Z_MAX 10							// Max coordination number
#define MAXBOUND 500.0f
#define CRACKED true

// Define constants
#define kB 1.38064852e-5					// Boltzmann constant
#define b 0.1								// Persistence length
#define T 300 								// Temperature
#define ae 0.1 								// Strength of bond - includes activation energy
#define delxe 0.15 							// parameter for breaking crosslink connection

static double vel[DIM] = {0.0, 100.0};

inline double getnorm(const double* vec, const int dim = DIM){
	double s = 0;
	for (int j = 0; j<DIM; j++){
		s += pow(vec[j], 2);
	}
	return sqrt(s);
}

inline void convert_to_vector(double* result, const double mag, const double* direction){
	for (int i = 0; i<DIM; i++){
		result[i] = mag*direction[i];
	}
}

inline void normalize_vector(double* result, const double* vec){
	double norm = getnorm(vec);
	for (int i = 0; i<DIM; i++){
		result[i] = vec[i] / norm;
	}
}

inline void normalize_vector(double* vec){
	double norm = getnorm(vec);
	for (int i = 0; i<DIM; i++){
		vec[i] = vec[i] / norm;
	}
}

inline double dist(const double* r1, const double* r2){
	double s = 0.0;
	for (int j = 0; j<DIM; j++){
		s += pow(r1[j] - r2[j], 2.0);
	}
	return sqrt(s);
}

inline void unitvector(double* result, double* r1, double* r2){
	for (int j = 0; j<DIM; j++){
		result[j] = r1[j] - r2[j];
	}
	normalize_vector(result);
}

inline double force_wlc(double x, double L){
	double t = x / L;
	if (t < 0.99){ return kB*T / b * (t + 1.0 / 4.0 / pow((1 - t), 2) - 1.0 / 4.0); }
	else { return 999999.0; }
}

inline double kfe(double force_mag){
	return ae * exp(force_mag * delxe / kB / T);
}

void forcevector(double* result, double* r1, double* r2, double L){
	double rhat[DIM];
	double s = dist(r1, r2);
	unitvector(rhat, r1, r2);
	double force = force_wlc(s, L);
	convert_to_vector(result, force, rhat);
}

double psi_wlc(double stretch, double L){
	return 3.0 / 4.0*kB*T / b*L*pow(L, 2)*(1.0 + 1.0 / 3.0 * stretch / (1 - stretch));
}

void delpsi_wlc_delrvec(double* result, double* r1, double* r2, double L){
	double stretch_vector[DIM];
	for (int i = 0; i < DIM; i++){
		stretch_vector[i] = (r1[i] - r2[i]) / L;
	}
	double stretch = getnorm(stretch_vector);
	double delpsi_delr_mag = 3.0*(1.0 + stretch / (1.0 - stretch) / 3.0) + stretch / 2.0 * pow(1.0 / (1.0 - stretch), 2);
	double normed_stretch_vector[DIM];
	normalize_vector(normed_stretch_vector, stretch_vector);
	convert_to_vector(result, delpsi_delr_mag, normed_stretch_vector);
}

template <typename t>
void print_array(t* arr, int n, int dim){
	for (int i = 0; i<n; i++){
		std::cout<<i<<"\t";
		for (int j = 0; j<dim; j++){
			std::cout << arr[i*dim + j] << "\t";
		}
		std::cout << std::endl;
	}
}

void get_pull_forces(double* , double* , int , const int* , int);

void get_forces(double* forces, double* R, int* edges, double* damage,\
	const double* chain_len, const int num_edges, const bool* PBC_STATUS,\
	const double* PBC_vector, bool update_damage = false){
	// if PBC_STATUS == True, assumes that the 
	// second node position = second_node_position + PBC_vector
	
	int node1, node2;
	int j, k, id; // loop variables
	double r1[DIM]; double r2[DIM] ;
	double edge_force[DIM];

	for (j = 0; j < num_edges; j++){
		// read the two points that form the edge // 2 because 2 points make an edge! Duh.
		node1 = edges[j * 2]; node2 = edges[j * 2 + 1];
		
		// check if pair exists
		if(node1 == -1 || node2 == -1){continue;}

		// read the positions
		for(k = 0; k<DIM; k++){
			r1[k] = R[node1*DIM + k]; 
			r2[k] = R[node2*DIM + k];
		}
		
		// check PBC_STATUS
		if (PBC_STATUS[j]) {
			// add PBC_vector to get new node position
			for (k = 0; k < DIM; k++){
				r2[k] += PBC_vector[k];
			}
			// get force on node1 due to node2
			forcevector(edge_force, r1, r2, chain_len[j]);
			// subtract back the PBC_vector to get original node position
			for (k = 0; k < DIM; k++){
				r2[k] -= PBC_vector[k];
			}
		}
		else{
			forcevector(edge_force, r1, r2, chain_len[j]);
		}
			for (k = 0; k < DIM; k++){
				// if (edge_force[k] > 99.0) {
				// 	cout<<"Node "<<node1<<" and node "<<node2<<endl;
				// 	cout<<"Distance is "<<dist(r1, r2)<<endl;
				// 	cout<<"L is "<<chain_len[j]<<endl;
				// 	cout<<endl;
				// }
			forces[node1*DIM + k] -= edge_force[k];
			forces[node2*DIM + k] += edge_force[k];
		}
		//update damage if needed
		if (update_damage){
			damage[j] += kfe(getnorm(edge_force))*TIME_STEP;
			//remove edge ... set to special value
			if(damage[j] > 1.0){cout<<"Breaking bond between "
				<<edges[j*2]<<" and "<<edges[2*j +1]<<endl;
			edges[j*2] = -1; edges[j*2+1] = -1;}
		}

	}
	//cout<<"Got here!\n";	

}



template <typename t>
inline bool ismember(t item, const t* array, const int size){
	bool is_inside = false; 
	for(int k=0; k<size; k++){if(array[k]==item){is_inside = true; break;}}
	return is_inside;
}

template <typename t>
t inline t_abs(t var){
	if(var<0.0){var *= -1;}
	return var;
}

template <typename t>
t inline getabsmax(t* arr, int sizeofarr){
	t max_arr = TOL;
	t current_item;
	for(int k=0; k<sizeofarr; k++){
		current_item = t_abs<t>(arr[k]);
		if(current_item > max_arr){max_arr=current_item;}
	}
	return max_arr;
}
template <typename t>
inline void zero_arr(t* arr, int n){
	for(int i=0; i< n; i++){
		arr[i] = 0;
	}
}

void optimize(double*R, int* edges, double* damage_integral, \
	const double* chain_len, const int num_nodes, const int num_edges, \
	const bool* PBC_STATUS, const double* PBC_vector, \
	const int* tnodes, int n_tnodes, const int* bnodes, int n_bnodes, \
	double* plate_force, int iter,\
	double eta = 0.1, double alpha = 0.9, int max_iter = 1000){
	// Momentum based gradient method
	/*
	Passive nodes do not participate in optimization
	*/
	double* forces = new double[num_nodes*DIM];
	double* delR = new double[num_nodes*DIM];
	double* rms_history = new double[num_nodes*DIM]();
	double* oneby_sqrt_rms_history = new double[num_nodes*DIM]();
	double g;
	for(int step = 0; step < max_iter; step++){
		int id = 0;
		zero_arr<double>(forces, num_nodes*DIM);
		get_forces(forces, R, edges, damage_integral, chain_len, \
		num_edges, PBC_STATUS, PBC_vector, false);
		// cout<<"Got upto this point\n";
		if (getabsmax<double>(forces, num_nodes*DIM) > TOL){
		// get update step
			for(id = 0; id < num_nodes*DIM; id++){
				if (ismember<int>(id/DIM, tnodes, n_tnodes) || ismember<int>(id/DIM, bnodes, n_bnodes)){
					id++;
					continue;
				}
				g = forces[id];
				rms_history[id] = alpha*rms_history[id] + (1.0-alpha)*g*g;
				oneby_sqrt_rms_history[id] = sqrt(1.0/(rms_history[id] + TOL));
				delR[id] = eta*oneby_sqrt_rms_history[id] * g;
				R[id] += delR[id]; 
			}
		}
		else{
			break;
		}
	}
	zero_arr<double>(forces, num_nodes*DIM);
	get_forces(forces, R, edges, damage_integral, chain_len, \
		num_edges, PBC_STATUS, PBC_vector, true);
	get_pull_forces(forces, plate_force, iter, tnodes, n_tnodes);

	delete[] forces;
	delete[] delR;
	delete[] rms_history;
	delete[] oneby_sqrt_rms_history;
}

void inline force_on_plate(double* plate_force, double* forces,\
	int* mNodes, int num_mnodes){
	int node, d = 0;
	for(d=0;d<DIM;d++){
		plate_force[d] =  0.0;
	}		
	for(d=0; d<num_mnodes; d++){
		node = mNodes[d];
		for(d=0;d<DIM;d++){
			plate_force[d] += forces[node*DIM + d];
		}
	}
}

void __init__(double* L, double* damage, bool* PBC, int num_elems){
	std::default_random_engine seed;
	std::normal_distribution<double> generator(L_MEAN, L_STD);
	for(int i=0; i<num_elems; i++){
		L[i] = generator(seed);
		damage[i] = 0.0;
		PBC[i] = false;
	}	
}

void make_edge_connections(double* R, int* edges, int& num_edges, \
	int* lsideNodes, int* rsideNodes, int num_lnodes, int num_rnodes, \
	bool* PBC_STATUS, double* chainlength, double* damage_integral, \
	double dely_allowed = 10.0){
	// Checks y coordinates of nodes and connects 
	// if nodes are close-enough
	std::default_random_engine seed;
	std::normal_distribution<double> generator(L_MEAN, L_STD);
	int nl, nr, lnode, rnode;
	for(nl= 0; nl < num_lnodes; nl++){
		lnode = lsideNodes[nl];
		for(nr= 0; nr < num_rnodes; nr++){
			rnode = rsideNodes[nr];
			if (fabs(R[lnode*DIM + 1] - R[rnode*DIM + 1]) < dely_allowed){
				edges[num_edges*2] = rnode;
				edges[num_edges*2 + 1] = lnode;
				cout<<"Connected node "<<lnode<<" and "<<rnode<<"\n";
				chainlength[num_edges] = generator(seed);
				damage_integral[num_edges] = 0.0;
				PBC_STATUS[num_edges] = true;
				num_edges += 1;
			}
		}
	}

}

void side_nodes(double* R,\
	int* lnodes, int* rnodes, int* tnodes, int* bnodes,\
	int& n_lside, int& n_rside, int& n_tside, int& n_bside, int n){
	// make the list of side nodes
	for(int i=0; i<n; i++){
		if(fabs(R[i*DIM]-0.0)<TOL){lnodes[n_lside]=i; n_lside++;}
		if(fabs(R[i*DIM]-MAXBOUND)<TOL){rnodes[n_rside]=i; n_rside++;}
		if(fabs(R[i*DIM + 1]-0.0)<TOL){bnodes[n_bside]=i; n_bside++;}
		if(fabs(R[i*DIM + 1]-MAXBOUND)<TOL){tnodes[n_tside]=i; n_tside++;}
	}
}

void move_top_nodes(double* R, int* tnodes, int n_tnodes){
	for(int i=0; i < n_tnodes; i++){
		R[tnodes[i]*DIM + 1] += vel[1] * TIME_STEP;
	}
}

// void convert_neighbors_for_cuda(double* cuda_props, int* neighbors, double* R,  const int n){
	
// 	std::default_random_engine generator;
//   	std::normal_distribution<double> distribution(L_MEAN,L_STD);
// 	int i,j,k = 0;
// 	int neighbor_id;
// 	int source_id, dest_id;
// 	for(i; i<n; i++){
// 		for(j; j<Z_MAX; j++){
// 			source_id = i*Z_MAX + j;
// 			if (neighbors[source_id] != -1){
// 				neighbor_id = neighbors[source_id];
// 				dest_id = i*Z_MAX*(DIM+1);
// 				for(k;k<DIM;k++){
// 					cuda_props[dest_id + k] = R[neighbor_id*DIM+k];
// 				}
// 				//Update L
// 				cuda_props[dest_id + DIM] = distribution(generator);
// 			}
// 		}
// 	}
// }

void get_file(const string& filename, double* arr, \
	int n_dim, int n_nodes){
	// Write simulation parameters
	std::time_t result = std::time(nullptr);

	cout<<"1D pulling of a 2D gel"<<endl;
	cout<<"File created at "<<std::asctime(std::localtime(&result));
	cout<<endl;

	cout<<"Sim dimension : "<<DIM<<endl;
	cout<<"Simulation time : "<<SIM_TIME<<endl;
	cout<<"Velocity : "<<vel[1]<<endl;
	cout<<endl;

	cout<<"Disorder characteristics : "<<endl;
	cout<<" -- L_MEAN : "<<L_MEAN<<endl;
	cout<<" -- L_STD : "<<L_STD<<endl;
	cout<<endl;
	
	cout<<"Cracked? : "<<CRACKED<<endl;

	ofstream fout;
	fout.open(filename, ios::trunc|ios_base::out);
	for(int i = 0; i < n_nodes; i++){
		for(int d = 0; d<n_dim; d++){
			fout<<arr[i*n_dim + d]<<"\t";
		}
		fout<<endl;
	}
	fout.close();
}

void plot_network(Gnuplot& h, double* R, int* edges, bool* PBC,\
	int n_nodes, int n_elems, int iter_step){

	// write to file
	ofstream f;
	f.open("data.txt");
	int node1, node2;
	for(int i = 0; i<n_elems; i++){
		node1 = edges[2*i];
		node2 = edges[2*i+1];
		if(node1!=-1 && node2!=-1){
			for(int d = 0; d<DIM; d++){
					f<<R[node1*DIM+d]<<"\t";
				}
				f<<endl;
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
		f<<endl<<endl;
		}
	}
	f.close();

	//Plot to h
	h.cmd("plot 'data.txt' every ::0::1 with lp title '"+\
		std::to_string(iter_step)+"'");
	h.cmd("set term png");
	h.cmd("set output '"+std::to_string(iter_step)+".png'");
	h.cmd("replot");
	h.cmd("set term x11");
	h.reset_plot();

	//delete temporary file
	//std::remove("data.txt");
}

void plot_forces(Gnuplot& h, \
	const vector<double>& p_x, const vector<double>& p_y, \
	int iter_step){
	
	h.remove_tmpfiles();
	h.reset_plot();
	h.plot_x(p_x, "fx");
	h.plot_x(p_y, "fy");
	h.set_xrange(0, iter_step+1200);
	h.set_yrange(-0.05,\
	 2*(*max_element(p_y.begin(), p_y.end())));	
}


void get_pull_forces(double* forces, double* plate_force,\
 int iter_step, const int* tside, int tnodes){
	// Zero the forces first
	plate_force[2*iter_step] = 0.0;
	plate_force[2*iter_step + 1] = 0.0;
	for(int i=0; i<tnodes; i++){
		plate_force[2*iter_step] += forces[2*tside[i]];
		plate_force[2*iter_step+1] += forces[2*tside[i]+1]; 
	}
}

void get_components(vector<double>& vec_x , vector<double>& vec_y,\
	 double* arr, int index){
	vec_x.push_back(-arr[2*index]);
	vec_y.push_back(-arr[2*index+1]);
}

inline bool contains(vector<int>& vec, int elem){
	return (std::find(vec.begin(), vec.end(), elem) != vec.end());
}

void crack(double* c, double* a, double* R, int n_nodes,
	int* edges, int n_edges){
	
	double equation = 0;
	vector<int> nodes_to_remove;
	
	for(int i=0; i<n_nodes; i++){
		equation = 0;
		for(int d=0; d<DIM; d++){
			equation += pow(R[i*DIM+d]-c[d],2)/pow(a[d],2);
		}
		equation -= 1.0;

		if(equation<=0.0){
			nodes_to_remove.push_back(i);
		}
	}
	int edges_removed = 0, node1, node2;
	for(int i = 0; i<n_edges; i++){
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

void gen_random(char* s, const int len){
	static const char alpha[]="ABCDEFGHIJKLMNOPQRSTUVWXYZ";

	for(int i = 6; i<len+6; i++){
		s[i] = alpha[rand()%(sizeof(alpha)-1)];
	}
	s[len] = 0;
}

int main(){
	int n_nodes=1600, n_elems = 11000;

	// Initialize nodes and edges
	double* R = (double*)malloc(n_nodes*DIM*sizeof(double));
	int* edges = (int*)malloc(Z_MAX*n_nodes*2*sizeof(int));
	double* pull_forces = (double*)malloc(STEPS*DIM*sizeof(double));

	vector<double> p_x;
	vector<double> p_y;
	// file to read and write position data
	const string fname = "coordinates.txt";

	// Read in mesh
	cout<<"Reading the mesh...\n";	
	take_input(R, edges, n_nodes, n_elems);
	cout<<"Mesh read successfully!\n";
	cout<<"Number of nodes are: "<<n_nodes<<endl;
	cout<<"Number of elements are: "<<n_elems<<endl;

	int max_nodes_on_a_side = int(sqrt(n_nodes))*2;
	int n_rside = 0, n_lside = 0, n_bside = 0, n_tside = 0;
	// get left and right nodes
	int* lsideNodes = (int*)malloc(max_nodes_on_a_side*sizeof(double));
	int* rsideNodes = (int*)malloc(max_nodes_on_a_side*sizeof(double));
	// get top and bottom nodes
	int* tsideNodes = (int*)malloc(max_nodes_on_a_side*sizeof(double));
	int* bsideNodes = (int*)malloc(max_nodes_on_a_side*sizeof(double));

	side_nodes(R, lsideNodes, rsideNodes, tsideNodes, bsideNodes,\
		n_lside, n_rside, n_tside, n_bside, n_nodes);	

	// Initialize edge properties
	double* damage = (double*)malloc(2*n_elems*sizeof(double));
	double* L = (double*)malloc(2*n_elems*sizeof(double));
	bool* PBC = (bool*)malloc(2*n_elems*sizeof(bool));
	__init__(L, damage, PBC, n_elems);

	// Make PBC connections
	const double PBC_vector[DIM] = {MAXBOUND*1.2, 0};
	// TODO: get lside and rside nodes
	make_edge_connections(R, edges, n_elems, \
		lsideNodes, rsideNodes, n_lside, n_rside, \
		PBC, L, damage, 15.0);
	cout<<"Number after new connections made: "<<n_elems<<endl;

	if(CRACKED){
		double c[] = {MAXBOUND/2.0, MAXBOUND/2.0};
		double a[] = {MAXBOUND/20.0, MAXBOUND/25.0};
		crack(c, a, R, n_nodes, edges, n_elems);
	}
	
	// gnuplots
	//Gnuplot gforces("lines lw 2");
	//Gnuplot gnetwork;

	//plot_network(gnetwork, R, edges, PBC, \
	//			n_nodes, n_elems, 0);
	// Track time for 1000 iterations
	clock_t t = clock(); 
	
	int iter = 0; // needed to write forces later

	for(iter = 0; iter<STEPS; iter++){
		if((iter+1)%1000 == 0){ // +1 required to have values in p_x, p_y
			cout<<(iter+1)<<endl; 
			cout<<"That took "<<(clock()-t)/CLOCKS_PER_SEC<<" s\n";
			t = clock();  // reset clock

			// plot network
			//plot_network(gnetwork, R, edges, PBC, \
			 	n_nodes, n_elems, iter);
			
			// Plot forces
			//plot_forces(gforces, p_x, p_y, iter);


		}
		optimize(R, edges, damage, L, n_nodes, n_elems,\
			PBC, PBC_vector, tsideNodes, n_tside, bsideNodes, n_bside,\
			pull_forces, iter);
	
		get_components(p_x, p_y, pull_forces, iter);

		move_top_nodes(R, tsideNodes, n_tside);

	}

	
	char rand_string[] = "forcesxxxxxx.txt";
	gen_random(rand_string, 6);
	
	get_file(rand_string, pull_forces, DIM, iter);

	// wait to check out the plots	
	char p;
	cin>>p;

	// free all variables
	free(pull_forces);
	free(bsideNodes); free(tsideNodes);
	free(lsideNodes); free(rsideNodes);
	free(PBC);
	free(R); free(edges);
	free(L); free(damage);	
	
	return 0;
}


	///////////////////////////////////////////////////////////
	// unit tests for functions
	///////////////////////////////////////////////////////////


	// double r1[] = {0.0,0.0};
	// double r2[] = {1.0, 1.0};
	// double result[] = {0.0,0.0};
	// double ch_l = 45;

	// unitvector(result,r1,r2);
	// cout<<"Unit vector is : \n";	
	// print_array<double>(result, 1, DIM);
	// cout<<"Distance is : "<<dist(r1, r2)<<endl;
	// forcevector(result, r1, r2 , ch_l);
	// print_array<double>(result, 1, DIM);

	//print_array<double>(forces, n_nodes, DIM);

	// cout<<"Pair\tchain length\t distance"<<endl;
	// int node1, node2;
	// for(int i=0; i<n_elems; i++){
	// 	node1 = edges[2*i]; node2 = edges[2*i+1];
	//  	for(int k = 0; k<DIM; k++){
	//  		r1[k] = R[node1*DIM + k];
	//  		r2[k] = R[node2*DIM + k];
	//  	}
	//  	cout<<"("<<node1<<", "<<node2<<")\n";
	//  	print_array<double>(r1, 1, DIM);
	//  	print_array<double>(r2, 1, DIM);
	//  	unitvector(result, r1, r2);
	//  	print_array<double>(result, 1, DIM);
	//  	cout<<"Distance, L, PBC =  "<<dist(r1, r2)<<", "<<L[i]<<", "<<PBC[i]<<endl;
	//  	cout<<"Force = "<<force_wlc(dist(r1,r2), L[i])<<endl;
	//  // 	cout<<"("<<node1<<", "<<node2<<")\t"\
	//  // 	<<forces[i*DIM]<<"\t"<<forces[i*DIM + 1]<<endl;
	// 	// //cout<<PBC[i]<<" \t";
	//  // 	cout<<L[i]<<" \t";
	//  // 	cout<<dist(r1,r2)<<endl;
	//  	cout<<endl;
	//  }

// Herbert Hui -- Nature Materials