/**
@file sac_network.cpp
\brief Implements all functions prototyped in sac_network.h
*/

#include "sac_network.h"

// ----------------------------------------------------------------------- 
/// \brief Default constructor
// ----------------------------------------------------------------------- 
sacNetwork::sacNetwork(){
	initialized = false;
}

// ----------------------------------------------------------------------- 
/// \brief Constructor that instantiates a sacNetwork object according
/// to given \emph{fname} file. Please note that there are no checks made
/// on the GMSH .msh filepath and that responsibility is left to the user.
///
/// \param fname --> takes a filename and instantiates Network object
///
// ----------------------------------------------------------------------- 
sacNetwork::sacNetwork(string& fname, bool from_dump){
	initialized = false;
	if(from_dump){
		load_from_dump(fname);
	}

	else{
		load_network(fname);
		iter_offset = 0;
	}

	initialized = true;
}


// ----------------------------------------------------------------------- 
/// \brief Default destructor
// ----------------------------------------------------------------------- 
sacNetwork::~sacNetwork(){
	clear();
}

// ----------------------------------------------------------------------- 
/// \brief Copy constructor
// ----------------------------------------------------------------------- 
sacNetwork::sacNetwork(sacNetwork const & source): Network(source){
	size_t sf = sizeof(float);
	size_t si = sizeof(int);
	m = (int*)malloc(si*n_elems);
	sacdamage = (float*)malloc(sf*n_elems);
	for(int d=0; d<n_elems; d++){
		m[d] = source.m[d];
		sacdamage[d] = 0.0f; 
	}
}


// ----------------------------------------------------------------------- 
/// \brief Clears extra variables declared by sacNetwork object. Called
/// by destructor.
// ----------------------------------------------------------------------- 
void sacNetwork::clear() {
	// DO NOT CLEAR baseclass' (Network) variables from here
	free(sacdamage);
	sacdamage = NULL;
	free(m);
	m = NULL;

}

// ----------------------------------------------------------------------- 
/// \brief Overrides parent class' function. 
/// Calculates forces along edges of the polymer network graph (i.e
/// forces in each polymer). The forces are then assigned to nodes which 
/// prepares the network for the optimization step. If update_damage is 
/// true, the appropriate damage metric is also updated. If damage exceeds
/// predefined thresholds, the edge is broken. Hidden lengths are also opened 
/// according to damage in these hidden bonds. The function does not return
/// anything as it updates the objects forces directly.
///
/// \param update_damage (bool) --> flag to update damage
// -----------------------------------------------------------------------
void sacNetwork::get_forces(bool update_damage = false) {

	int node1, node2;
	int j, k, id; // loop variables
	float r1[DIM]; float r2[DIM] ;
	float edge_force[DIM];
	float rhat[DIM];
	float s;
	float force;


	memset(forces, 0.0, n_nodes*DIM*sizeof(forces));

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
			convert_to_vector(edge_force, force, rhat);
		}
		#pragma unroll
		for (k = 0; k < DIM; k++){
			forces[node1*DIM + k] -= edge_force[k];
			forces[node2*DIM + k] += edge_force[k];
		}


		//update damage if needed
		if (update_damage){
			if(m[j] > 0){
				if(RATE_DAMAGE){
					sacdamage[j] += kf(force)*TIME_STEP;
				}
				else {
					sacdamage[j] = s/L[j]/.72;
				}

				if(sacdamage[j] > 1.0){
					if(weight_multiplier*L_MEAN > L[j]){
						L[j] += (weight_multiplier*L_MEAN - L[j])/m[j];
					}  
					m[j] -= 1;
					sacdamage[j] = 0.0;
					force = force_wlc(s, L[j]);
				}
			}
			if(RATE_DAMAGE){
				damage[j] += kfe(force)*TIME_STEP;
			}
			else {
				damage[j] = s/L[j]/0.9;
			}
			
			//remove edge ... set to special value
			if(damage[j] > 1.0){
				cout<<"Breaking bond between "
				<<edges[j*2]<<" and "<<edges[2*j +1]<<" F,s/L = "<<force \
				<<", "<<s/L[j]<<endl;
				edges[j*2] = -1; edges[j*2+1] = -1;
			}
		}
	}

}

// ----------------------------------------------------------------------- 
/// \brief Extends parent class' function to add memory allocation for extra
/// class variables.
///
// -----------------------------------------------------------------------
void sacNetwork::malloc_network(){
	
	Network::malloc_network();
	size_t sf = sizeof(float);
	size_t si = sizeof(int);
	m = (int*)malloc(n_elems*si);
	sacdamage = (float* )malloc(n_elems*sf);
}


// ----------------------------------------------------------------------- 
/// \brief Overrides parent class' function to initialize extra variables.
/// Called by sacNetwork(string& fname).
///
// -----------------------------------------------------------------------
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

	// read n_nodes, n_elems
	read_n(n_nodes, n_elems, fname);
	
	//Malloc all variables
	malloc_network();

	cout<<"Malloc was successful!\n";

	cout<<"Reading the mesh...\n";
	take_input(R, edges, n_nodes, n_elems, fname);

	// remove duplicates
	remove_duplicates(n_elems);

	cout<<"Mesh read successfully!\n";
	cout<<"Number of nodes are: "<<n_nodes<<endl;
	cout<<"Number of elements are: "<<n_elems<<endl;

	side_nodes();


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