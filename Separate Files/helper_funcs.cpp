/**
@file helper_functions.cpp
\brief Implements functions prototyped in helper_functions.h. Functions
implemented here are not expected to be called by user. Hence, only brief 
descriptions without parameter descriptions is attached.
*/

// Contains all the helper functions needed

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

using namespace std;

#include "vel.h"
#include "helper_funcs.h"


// ----------------------------------------------------------------------- 
/// \brief Converts mesh elements from GMSH file to edges and nodes. Only 
/// implemented elements of <element-id> = 2 or 3. Need to extend for 3D.
///
// -----------------------------------------------------------------------
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

// ----------------------------------------------------------------------- 
/// \brief Reads the number of edges in elements. Only implemented elements
/// of <element-id> = 2 or 3.
///
// -----------------------------------------------------------------------
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


// ----------------------------------------------------------------------- 
/// \brief Reads the number of nodes and poly chains left at the end of 
/// previous simulation. The code does not take care of parameters. Assumes
/// user will change params.h file while restarting the previous simulation.
///
// -----------------------------------------------------------------------
string read_dump(int& n_nodes, int& n_elems, string& dname){

	// reset n_nodes and n_elems
	n_nodes = 0;
	n_elems = 0;

	// this will return the last iter stored in the dump
	string ret_last_iter;

	bool read=false;
	string line;
	ifstream dump;
	dump.open(dname);

	getline(dump, line);

	while(!dump.eof()){
		getline(dump, line);
		if(line.find("INFO FOR ITER") != string::npos){
			getline(dump, line);
			ret_last_iter = line;

		}
	}



	// Start at the beginning of the file
	dump.close();
	dump.open(dname);


	// Get to the last iter
	while(!dump.eof()){
		getline(dump, line); 

		if(line == ret_last_iter){ // get to last iter

			getline(dump,line); // skip line START NODE POSITIONS
			getline(dump,line); // get to the first node
			if(dump.eof()){
				break;
			}

			while(line.find("END NODE POSITIONS") == string::npos){
				getline(dump, line);
				n_nodes++;
			}
			// tell user how many nodes were found -- debug
			cout<<"Found "<<n_nodes<<" nodes in the dump!"<<endl;

			getline(dump,line); // skip line START ACTIVE EDGES
			getline(dump,line); // get to the first edge
			while(line.find("END ACTIVE EDGES") == string::npos){
				getline(dump, line);
				n_elems++;
			}
			// tell user how many nodes were found -- debug
			cout<<"Found "<<n_elems<<" elems in the dump!"<<endl;

			// exit loop
			break;
		}
	}
	dump.close();
	return ret_last_iter;
}

// ----------------------------------------------------------------------- 
/// \brief Reads the number of nodes and edges according to the .msh file.
/// Useful for allocating memory in the class' constructors.
///
// -----------------------------------------------------------------------
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

// ----------------------------------------------------------------------- 
/// \brief Reads the input and converts to appropriate network data
///
// -----------------------------------------------------------------------
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


// ----------------------------------------------------------------------- 
/// \brief Initializes various properties for each edge in the graph
///
// -----------------------------------------------------------------------
void __init__(float* L, float* damage, bool* PBC, int n_elems){
	std::default_random_engine seed;
	seed.seed(std::chrono::system_clock::now().time_since_epoch().count());
	// std::uniform_real_distribution<float> generator(L_MEAN - L_STD, L_MEAN + L_STD);
	std::normal_distribution<float> generator(L_MEAN, L_STD);
	#pragma unroll
	for(int i=0; i<n_elems; i++){
		L[i] = generator(seed);
		damage[i] = 0.0;
		PBC[i] = false;
	}	
}

// ----------------------------------------------------------------------- 
/// \brief Gives a force vector given the nodes and the contour length of the 
/// polymer connecting these.
///
// -----------------------------------------------------------------------
void forcevector(float* result, float* r1, float* r2, float L){
	float rhat[DIM];
	float s = dist(r1, r2);
	unitvector(rhat, r1, r2);
	float force = force_wlc(s, L);
	convert_to_vector(result, force, rhat);
}



// ----------------------------------------------------------------------- 
/// \brief Gets the absolute maximum of an array
///
// -----------------------------------------------------------------------
float getabsmax(float* arr, size_t sizeofarr){
	float max_elem = 0.0;
	for(int i = 0; i<sizeofarr; i++){
		if(fabs(arr[i])>max_elem){max_elem = fabs(arr[i]);}
	}
	return max_elem;
}



// ----------------------------------------------------------------------- 
/// \brief Overrides __init__ for sacNetworks.
///
// -----------------------------------------------------------------------
void __init__(float* L, int* m, float* damage, float* sacdamage, bool* PBC, int n_elems){
	std::default_random_engine seed;
	seed.seed(std::chrono::system_clock::now().time_since_epoch().count());
	std::normal_distribution<float> generator(L_MEAN, L_STD);
	// std::uniform_real_distribution<float> generator(L_MEAN - L_STD, L_MEAN + L_STD);
	// std::uniform_real_distribution<float> hidden(0.10*L_MEAN, 0.20*L_MEAN);
	std::uniform_real_distribution<float> sacornosac(0.0,1.0);
	std::uniform_real_distribution<float> hidden(0.08*L_MEAN, 0.10*L_MEAN);
	std::uniform_int_distribution<int> number(4, 6);
	#pragma unroll
	float p = 1.0; // percentage of bonds that will be having sacrificial bonds
	for(int i=0; i<n_elems; i++){
		L[i] = generator(seed);
		if (sacornosac(seed) < p){
			m[i] = number(seed);
		}
		else{
			m[i] = 0;
		}
		L[i] -= m[i]*hidden(seed); 
		sacdamage[i] = 0.0;
		damage[i] = 0.0;
		PBC[i] = false;
	}	
}


// ----------------------------------------------------------------------- 
/// \brief Gives an integer value to be appended to the network dump 
/// filename
///
// -----------------------------------------------------------------------
int filename(const string& fname){
	string temp = fname + ".txt";
	std::ifstream f(temp); // Associate stream f to file
	bool exists = f.good();
	int c = 0;
	if(!exists){
		return -1;
	}
	else{
		while(f.good()){
			c++;
			f.close(); //Dissociate stream from file
			temp = fname + "_" + std::to_string(c) + ".txt"; // Overwrite temp
			cout<<temp<<endl;
			f.open(temp);
		}
		return c;
	}
}