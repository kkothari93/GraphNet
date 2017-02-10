#include <iostream>
#include <cmath>
#include <cstdlib>
#include <random>
#include <fstream>
#include <cmath>
#include <time.h>
#include <vector>
#include <cstring>
#include <algorithm>
#include <ctime>
#include <chrono>

#include "gnuplot_i.hpp"
#include "input.h"
#include "vel.h"

#define __params__

#define DIM 2								// Number of dimensions
#define TIME_STEP 1e-4						// Time step
#define SIM_TIME 10.0						// Simulation time
#define TOL 1e-6							// Tolerance
#define STEPS int(SIM_TIME/TIME_STEP)		// Number of time steps
#define L_MEAN 120.0f						// Average for contour length
#define L_STD 4.0f							// Std. deviation for contour lengths
#define Z_MAX 10							// Max coordination number
#define MAXBOUND 500.0f
#define CRACKED false
#define USE_CUDA true

// Define constants
#define __constants__

#define kB 1.38064852e-5					// Boltzmann constant
#define b 0.1								// Persistence length
#define T 300 								// Temperature
#define ae 0.1 								// Strength of bond - includes activation energy
#define delxe 0.15 							// parameter for breaking crosslink connection

#include <cuda.h>
#include <cuda_runtime_api.h>
#include "cudakernels.h"

static float vel[2] = {vel_x, vel_y};
// Make PBC connections
const float PBC_vector[DIM] = {MAXBOUND*1.1, 0};

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


void __init__(float* L, float* damage, bool* PBC, int num_elems){
	std::default_random_engine seed;
	std::normal_distribution<float> generator(L_MEAN, L_STD);
	for(int i=0; i<num_elems; i++){
		L[i] = generator(seed);
		damage[i] = 0.0;
		PBC[i] = false;
	}	
}

void make_edge_connections(float* R, int* edges, int& num_edges, \
	int* lsideNodes, int* rsideNodes, int num_lnodes, int num_rnodes, \
	bool* PBC_STATUS, float* chainlength, float* damage_integral, \
	float dely_allowed = 10.0){
	// Checks y coordinates of nodes and connects 
	// if nodes are close-enough
	std::default_random_engine seed;
	std::normal_distribution<float> generator(L_MEAN, L_STD);
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
}




void get_file(const string& filename, float* arr, \
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
	cout<<"Written to file!\n"<<filename<<"\n";
	fout.close();
}

void plot_network(Gnuplot& h, float* R, int* edges, bool* PBC,\
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
	const vector<float>& p_x, const vector<float>& p_y, \
	int iter_step){
	
	h.remove_tmpfiles();
	h.reset_plot();
	h.plot_x(p_x, "fx");
	h.plot_x(p_y, "fy");
	h.set_xrange(0, iter_step+1200);
	h.set_yrange(-0.05,\
	 2*(*max_element(p_y.begin(), p_y.end())));	
}


inline bool contains(vector<int>& vec, int elem){
	return (std::find(vec.begin(), vec.end(), elem) != vec.end());
}

void crack(float* c, float* a, float* R, int n_nodes,
	int* edges, int n_edges){
	
	float equation = 0;
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

void get_moving_nodes(int* moving_nodes, int& n_moving, \
	int* tnodes, int n_tside, int* bnodes, int n_bside,\
	int n_nodes){
	int c = 0, nt,nb;
	bool found = false;
	for(int i=0; i<n_nodes; i++){
		found =false;
		for(nt= 0; nt< n_tside; nt++){
			if(tnodes[nt]==i){
				found=true; 
				break;
			}
		}
		for(nb= 0; nb< n_bside; nb++){
			if(bnodes[nb]==i){
				found=true;
				break;
			}
		}
		if(!found){
			moving_nodes[c] = i;
			c++;
		}
	}
	cout<<"We have "<<c<<" moving nodes!\n";

}
int main(){

	int n_nodes, n_elems;
	int steps = 100;

	string fname = "./template2d.msh";
	//Check if file exists
	// bool exists = does_file_exist(fname);
	// if(!exists){
	// 	cout<<"File does not exist!\n";
	// 	return 0;
	// }
	read_n(n_nodes, n_elems, fname);

	int max_nodes_on_a_side = int(sqrt(n_nodes)*2.0);
	// add memory for side connections
	n_elems += 3*max_nodes_on_a_side;

	size_t sf = sizeof(float);
	size_t si = sizeof(int);
	size_t sb = sizeof(bool);

	float* R; int* edges;
	float* forces; float* damage;
	int* tsideNodes;
	int* bsideNodes;
	int* rsideNodes;
	int* lsideNodes;
	int* moving_nodes;
	bool* PBC;
	float* L; 
	float* pull_forces;

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
	pull_forces = (float* )malloc(steps*DIM*sf);

	std::memset(pull_forces, 0.0, steps*DIM*sf);

	// 	initialise n_xside for side nodes
	int n_rside = 0;
	int n_lside = 0;
	int n_bside = 0;
	int n_tside = 0;

	// adjust n_elems back to actual
	n_elems -= 3*max_nodes_on_a_side;


	// Read in mesh
	cout<<"Reading the mesh...\n";	
	take_input(R, edges, n_nodes, n_elems, fname);
	cout<<"Mesh read successfully!\n";
	cout<<"Number of nodes are: "<<n_nodes<<endl;
	cout<<"Number of elements are: "<<n_elems<<endl;

	side_nodes(R, lsideNodes, rsideNodes, tsideNodes, bsideNodes,\
		n_lside, n_rside, n_tside, n_bside, n_nodes);

	// get free (no BC) nodes 
	int n_moving = n_nodes - n_tside - n_bside;
	moving_nodes = (int*)malloc(n_moving*sizeof(int));
	get_moving_nodes(moving_nodes, n_moving, tsideNodes, n_tside,\
		bsideNodes, n_bside, n_nodes);


	// Initialize edge properties
	__init__(L, damage, PBC, n_elems);
	
	make_edge_connections(R, edges, n_elems, \
		lsideNodes, rsideNodes, n_lside, n_rside, \
		PBC, L, damage, 10.0);
	cout<<"Number after new connections made: "<<n_elems<<endl;

	if(CRACKED){
		float c[] = {MAXBOUND/2.0, MAXBOUND/2.0};
		float a[] = {MAXBOUND/20.0, MAXBOUND/25.0};
		crack(c, a, R, n_nodes, edges, n_elems);
	}
	
	// GPU copy packing
	hostvars vars;
	vars.PBC_vector = PBC_vector;
	vars.R = R; 
	vars.edges = edges;
	vars.forces = forces; 
	vars.tsideNodes = tsideNodes;
	vars.moving_nodes = moving_nodes;
	vars.PBC = PBC;
	vars.L = L; 
	vars.damage = damage;
	vars.pull_forces = pull_forces;
	vars.n_nodes = n_nodes; 
	vars.n_elems = n_elems;
	vars.n_tnodes = n_tside; 
	vars.n_moving = n_moving;

	// gnuplots
	// Gnuplot gforces("lines lw 2");
	// Gnuplot gnetwork;

	// plot_network(gnetwork, R, edges, PBC, \
	// 			n_nodes, n_elems, 0);
	// Track time for 1000 iterations
	clock_t t = clock(); 



	pull_CUDA(&vars, steps);

	cout<<"Simulation DONE! Storing data...\n";

	string filename = "forces.txt";
	
	get_file(filename, pull_forces, DIM, steps);

	// wait to check out the plots	

	// free all variables on CPU
	free(pull_forces);
	free(bsideNodes); free(tsideNodes);
	free(lsideNodes); free(rsideNodes);
	free(moving_nodes);
	free(PBC);
	free(R); free(edges);
	free(L); free(damage);	

	return 0;
}