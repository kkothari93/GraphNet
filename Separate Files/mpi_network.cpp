#include "mpi_network.h"

MPI_Network::MPI_Network(): Network() {

}

MPI_Network::MPI_Network(MPI_Network const & source) {

	initialized = false;
	copy(source);
	initialized = true;

}

MPI_Network::MPI_Network(Network const & source): Network(source) {

}

MPI_Network::MPI_Network(sacNetwork const & source): Network(source) {


}

MPI_Network::~MPI_Network() {
	clear();
}

void MPI_Network::clear() {
	Network::clear();
	free(not_moving_nodes);
	not_moving_nodes = NULL;
	//cout << __LINE__ << endl;
	delete[] chunk_nodes;
	chunk_nodes = NULL;
	delete[] chunk_edges;
	chunk_edges = NULL;
}

void MPI_Network::copy(MPI_Network const & source) {

	Network::copy(source);

	n_not_moving = source.n_not_moving;
	not_moving_nodes = (int*)malloc(sizeof(int)*n_not_moving);
	chunk_edges_len = source.chunk_edges_len;
	chunk_edges = new int[chunk_edges_len];
	chunk_nodes_len = source.chunk_nodes_len;
	chunk_nodes = new int[chunk_nodes_len];

	for (int i = 0; i < n_not_moving; i++) {
		not_moving_nodes[i] = source.not_moving_nodes[i];
	}
	for (int i = 0; i < chunk_edges_len; i++) {
		chunk_edges[i] = source.chunk_edges[i];
	}
	for (int i = 0; i < chunk_nodes_len; i++) {
		chunk_nodes[i] = source.chunk_nodes[i];
	}

}

void MPI_Network::load_network(string& fname){


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

	// cout<<"Malloc was successful!\n";

	// cout<<"Reading the mesh...\n";
	take_input(R, edges, n_nodes, n_elems, fname);
	// cout<<"Mesh read successfully!\n";
	// cout<<"Number of nodes are: "<<n_nodes<<endl;
	// cout<<"Number of elements are: "<<n_elems<<endl;

	side_nodes(R, lsideNodes, rsideNodes, tsideNodes, bsideNodes, \
		n_lside, n_rside, n_tside, n_bside, n_nodes);


	// NOT moving nodes
	n_not_moving = n_tside + n_bside;
	not_moving_nodes = (int*)malloc(sizeof(int)*n_not_moving);
	
	int c = 0;
	for(int i =0; i<n_not_moving; i++){
		if(i>n_tside){
			not_moving_nodes[i] = bsideNodes[i-n_tside];
		}
		else{
			not_moving_nodes[i] = tsideNodes[i];
		}
	}
	//cout<<"We have "<<c<<" moving nodes\n";


	//cout<<"Side nodes written successfully! \n";
	__init__(L, damage, PBC, n_elems);

	this->make_edge_connections(15.0);
	cout<<"Number after new connections made: "<<n_elems<<endl;


}

inline bool Rinset(float* R, float x_lo, float x_hi, float y_lo, float y_hi){
	return (R[0] >= x_lo && R[0] < x_hi && R[1] >= y_lo && R[1] < y_hi);
}

void MPI_Network::init_MPI(int world_rank, int world_size) {

	int y_level = world_rank % 2 == 0? world_rank/2 : (world_rank-1)/2;
	float x_lo = world_rank%2 == 0? 0 : (PAD/2);
	float x_hi = world_rank%2 == 0? PAD/2 : PAD;
	float y_lo = ((PAD*2.0/world_size) * y_level);
	float y_hi = (world_rank == world_size - 1) || (world_rank == world_size - 2)? PAD : ((PAD*2.0/world_size) * (y_level+1));

	chunk_nodes_len = n_nodes/world_size + sqrt(n_nodes)*2;
	chunk_nodes = new int[chunk_nodes_len];
	int k = 0;
	for (int i = 0; i < n_nodes; i+=1) {
		if (Rinset(&(R[i*DIM]),x_lo, x_hi, y_lo, y_hi)) {
			chunk_nodes[k] = i;
			k++;
		}
	}
	for (; k < chunk_nodes_len; k++) {
		chunk_nodes[k] = -1;
	}
	chunk_edges_len = n_elems/world_size + sqrt(n_elems)*5;
	chunk_edges = new int[chunk_edges_len];
	k = 0;
	for (int i = 0; i < n_elems; i+=1) {

		// check this condition for both edges, use Rinset
		if (!Rinset(&(R[edges[i*2]*DIM]),x_lo, x_hi, y_lo, y_hi) && !Rinset(&(R[edges[i*2+1]*DIM]),x_lo, x_hi, y_lo, y_hi)) {
			continue;
		}
		else {
			chunk_edges[k] = i;
			k++;
		}
	}
	for ( ; k < chunk_edges_len; k++) {
		chunk_edges[k] = -1;
	}

}

void MPI_Network::get_forces(bool update_damage = false) {

	int node1, node2;
	int j, k, l, id; // loop variables
	float r1[DIM]; float r2[DIM] ;
	float edge_force[DIM];
	float rhat[DIM];
	float s;
	float force;


	memset(forces, 0, n_nodes*DIM*sizeof(float));

	
	for (l = 0; l < chunk_edges_len; l++){
		// read the two points that form the edge // 2 because 2 points make an edge! Duh.
		//node1 = edges[j * 2]; 
		// chunk edges is the index of the edge array and not the edge itself
		j = chunk_edges[l];
		if(j==-1){
			break;
		}
		node1 = edges[j * 2];
		node2 = edges[(j * 2) + 1];
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
			// if (s>L[j]){
			// 	cout<<"Something's wrong for edge "\
			// 	<<j<< ", nodes "<<node1 << " - " <<\
			// 	node2<<": s = "<<s<<", L = "<<L[j]<<"\n";
			// }
			force = force_wlc(s, L[j]);
			if(force == 999999){edges[j*2] = -1; edges[j*2 +1] = -1; force =0.0; damage[j] = 0.0;}
			convert_to_vector(edge_force, force, rhat);
		}
		else{
			s = dist(r1, r2);
			// if (s>L[j]){
			// 	cout<<"Something's wrong for edge "<<j<< ", nodes "<<node1 << " - " <<node2<<": s = "<<s<<", L = "<<L[j]<<"\n";
			// }
			unitvector(rhat, r1, r2);
			force = force_wlc(s, L[j]);
			if(force == 999999){edges[j*2] = -1; edges[j*2 +1] = -1; force =0.0; damage[j] = 0.0;}
			convert_to_vector(edge_force, force, rhat);
		}
		#pragma unroll
		for (k = 0; k < DIM; k++){
			forces[node1*DIM + k] -= edge_force[k];
			forces[node2*DIM + k] += edge_force[k];
		}
		//update damage if needed
		if (update_damage) { //&& (node1 < node2)){
			damage[j] += kfe(force)*TIME_STEP;
			//remove edge ... set to special value
			if(damage[j] > 1.0){
				cout<<"Breaking bond between "
				<<node1<<" and "<<node2<<" F,s = "<<force \
				<<", "<<s<<endl;
				for(int d = 0; d<DIM; d++){
					cout<<r1[d]<<"\t "<<r2[d]<<"\t";
				}
				cout<<endl;
				edges[j*2] = -1; edges[j*2+1] = -1;
			}

		}
	}
	//printf("forces: %f\t%f\n", forces[4], forces[5]);
	return;

}

bool MPI_Network::notmoving(int nodeid){
	bool init = false;
	for(int i=0;i<n_not_moving;i++){
		if(nodeid == not_moving_nodes[i]){
			init = true;
			break;
		}
	}
	return init;
}

void MPI_Network::optimize(float eta, float alpha, int max_iter) {

	float* rms_history = new float[n_nodes*DIM](); // () allows 0.0 initialization
	float g, delR;
	int id, d, node;
	for(int step = 0; step < max_iter; step++){
		get_forces(false);
	
		if(getabsmax(forces,n_nodes*DIM)>TOL){
			for(id = 0; id < chunk_nodes_len; id++){
				node = chunk_nodes[id];
				if (chunk_nodes[id] ==-1 || notmoving(chunk_nodes[id])) {
						continue;
				}				
				#pragma unroll
					for(d = 0; d<DIM; d++){
						g = forces[DIM*node+d];
						rms_history[id*DIM + d] = alpha*rms_history[id*DIM + d] + (1-alpha)*g*g;
						delR = sqrt(1.0/(rms_history[id*DIM + d] + TOL))*eta*g;
						R[node*DIM + d] += delR;
					}				
			}
		}
		else{
			break;
		}
	}
	delete[] rms_history;
	rms_history = NULL;

}



