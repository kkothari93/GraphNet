#include "sac_network.h"

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

//TODO: Copy constructor for sac network
sacNetwork(sacNetwork const & source): Network(source){
}


void sacNetwork::clear() {
	// DO NOT CLEAR baseclass' (Network) variables from here
	// free(R);
	// free(edges);
	// free(forces);
	// free(damage);
	free(sacdamage);
	sacdamage = NULL;
	free(m);
	m = NULL;
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
				<<edges[j*2]<<" and "<<edges[2*j +1]<<" F,s/L = "<<force \
				<<", "<<s/L[j]<<endl;
				edges[j*2] = -1; edges[j*2+1] = -1;
		}
		}
	}

}

void sacNetwork::malloc_network(string& fname){
	
	Network::malloc_network(fname);
	size_t sf = sizeof(float);
	size_t si = sizeof(int);
	m = (int*)malloc(n_elems*si);
	sacdamage = (float* )malloc(n_elems*sf);
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