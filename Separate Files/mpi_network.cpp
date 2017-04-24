/**
@file network.cpp
\brief This file implements all functions defined in header mpi_network.h
*/
#include "mpi_network.h"

// ----------------------------------------------------------------------- 
/// \brief Default constructor
// ----------------------------------------------------------------------- 
MPI_Network::MPI_Network(): Network() {

}

// ----------------------------------------------------------------------- 
/// \brief Copy constructor
// ----------------------------------------------------------------------- 
MPI_Network::MPI_Network(MPI_Network const & source) {

	initialized = false;
	copy(source);
	initialized = true;

}

// ----------------------------------------------------------------------- 
/// \brief Constructor to initialize from a Network object
/// Note: This constructor uses initialization lists for initializing the 
/// member data.
// ----------------------------------------------------------------------- 
MPI_Network::MPI_Network(Network const & source): Network(source) {
}

// ----------------------------------------------------------------------- 
/// \brief Constructor to initialize from a sacNetwork object
/// Note: Since MPI_Network virtually (public) inherits from Network and 
/// publically inherits from sacNetwork (which itself virtual public 
/// inherits from Network), we need to initialize the base class (Network)
/// object too. The constructors for virtual base classes are always called
/// from the most derived class, using any arguments it might pass in. So,
/// sacNetwork copy constructor will not call the base Network copy constructor.
/// Therefore, we need to construct this object using base Network's copy 
/// constructor explicitly.
// -----------------------------------------------------------------------
MPI_Network::MPI_Network(sacNetwork const & source): Network(source), sacNetwork(source) {
}

// ----------------------------------------------------------------------- 
/// \brief Destructor
///
// -----------------------------------------------------------------------
MPI_Network::~MPI_Network() {
	clear();
}

// ----------------------------------------------------------------------- 
/// \brief Clears the extra variables instantiated in the MPI_Network object.
/// Called by the destructor of the class.
///
// -----------------------------------------------------------------------
void MPI_Network::clear() {
	// free(not_moving_nodes);
	// not_moving_nodes = NULL;
	//cout << __LINE__ << endl;
	delete[] chunk_nodes;
	chunk_nodes = NULL;
	delete[] chunk_edges;
	chunk_edges = NULL;

}

// ----------------------------------------------------------------------- 
/// \brief Copies source's data members into current object. Called by copy
/// constructor. 
///
// -----------------------------------------------------------------------
void MPI_Network::copy(MPI_Network const & source) {

	if(SACBONDS){sacNetwork::copy(source);}
	else{Network::copy(source);}

	// n_not_moving = source.n_not_moving;
	// not_moving_nodes = (int*)malloc(sizeof(int)*n_not_moving);
	chunk_edges_len = source.chunk_edges_len;
	chunk_edges = new int[chunk_edges_len];
	chunk_nodes_len = source.chunk_nodes_len;
	chunk_nodes = new int[chunk_nodes_len];

	// for (int i = 0; i < n_not_moving; i++) {
	// 	not_moving_nodes[i] = source.not_moving_nodes[i];
	// }
	for (int i = 0; i < chunk_edges_len; i++) {
		chunk_edges[i] = source.chunk_edges[i];
	}
	for (int i = 0; i < chunk_nodes_len; i++) {
		chunk_nodes[i] = source.chunk_nodes[i];
	}

}

// ----------------------------------------------------------------------- 
/// \brief Checks if a node lies in a particular partitioning of the graph
///
/// \param R --> A float array for a nodes position
/// \param x_lo --> Takes the x-dimension of the left side of the partition
/// \param x_hi
/// \param y_lo
/// \param_y_hi
// -----------------------------------------------------------------------
inline bool Rinset(float* R, float x_lo, float x_hi, float y_lo, float y_hi){
	return (R[0] >= x_lo && R[0] < x_hi && R[1] >= y_lo && R[1] < y_hi);
}

// ----------------------------------------------------------------------- 
/// \brief Rewrites moving nodes to include only nodes included in object's
/// partition
///
// -----------------------------------------------------------------------
void MPI_Network::rewrite_moving_nodes(){
	// Takes intersection of base class' moving nodes 
	// with chunk_nodes of MPI_class
	int chunk_n_moving = 0;
	int chunk_moving[n_moving];
	bool found;
	for(int i=0; i<chunk_nodes_len; i++){
		for(int d=0; d<n_moving; d++){
			if(chunk_nodes[i]==moving_nodes[d]){
				chunk_moving[chunk_n_moving] = chunk_nodes[i]; 
				chunk_n_moving++;
			}
		}
	}
	cout<<"Found "<<chunk_n_moving<<" nodes!\n";
	n_moving = chunk_n_moving;
	for(int d= 0; d<n_moving; d++){
		moving_nodes[d] = chunk_moving[d];
	}
}

// ----------------------------------------------------------------------- 
/// \brief 	Chunks the graph using squares and ranking the squares using 
/// the process' rank in the MPI execution
///
/// \param world_rank --> rank of the MPI process
/// \param world_size --> size of the MPI_Comm
// -----------------------------------------------------------------------
void MPI_Network::init_MPI(int world_rank, int world_size) {

	int y_level = world_rank % 2 == 0? world_rank/2 : (world_rank-1)/2;
	float x_lo = world_rank%2 == 0? 0 : (PAD_X/2);
	float x_hi = world_rank%2 == 0? PAD_X/2 : PAD_X;
	float y_lo = ((PAD_Y*2.0/world_size) * y_level);
	float y_hi = (world_rank == world_size - 1) || (world_rank == world_size - 2)? PAD_Y : ((PAD_Y*2.0/world_size) * (y_level+1));
	
	
	chunk_nodes_len = int(n_nodes/world_size*1.3);
	chunk_nodes = new int[chunk_nodes_len];
	
	// initialization can be inside for loop for k, increments need to be inside if loop
	int k = 0;
	for (int i = 0; i < n_nodes; i+=1) {
		if (Rinset(&(R[i*DIM]),x_lo, x_hi, y_lo, y_hi)) {
			chunk_nodes[k] = i;
			k++;
		}
	}

	for(; k< chunk_nodes_len; k++){
		chunk_nodes[k] = -1;
	}

	chunk_edges_len = int(n_elems/world_size*2);
	chunk_edges = new int[chunk_edges_len];
	int node1 , node2;
	k = 0;
	for (int i = 0; i < n_elems; i+=1) {
		node1 = edges[i*2];
		node2 = edges[i*2 + 1];

		// check this condition for both nodes, use Rinset
		if(node1 != -1 && node2 != -1){
			bool node1inset = Rinset(&(R[node1*DIM]),x_lo, x_hi, y_lo, y_hi);
			bool node2inset = Rinset(&(R[node2*DIM]),x_lo, x_hi, y_lo, y_hi);
			if (!node1inset && !node2inset) {
				continue;
			}
			else {
				chunk_edges[k] = i;
				k++;
			}
		}
	}
	for(;k<chunk_edges_len;k++){
		chunk_edges[k] = -1;
	}


	rewrite_moving_nodes();

}

// ----------------------------------------------------------------------- 
/// \brief 	Overrides the get_forces of parent classes (Network or sacNetwork)
/// to calculate forces in its own chunk
///
/// \param update_damage --> described in Network::get_forces
// -----------------------------------------------------------------------
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
			force = force_wlc(s, L[j]);
			convert_to_vector(edge_force, force, rhat);
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
		//update damage if needed
		if (update_damage){
			if(RATE_DAMAGE){
				if(SACBONDS){
					if(m[j] > 0){
						sacdamage[j] += kf(force)*TIME_STEP;
					if(sacdamage[j] > 1.0){
						L[j] += (L_MEAN - L[j])/m[j] ; 
						m[j] -= 1;
						sacdamage[j] = 0.0;
						force = force_wlc(s, L[j]);
						}
					}
				}

				damage[j] += kfe(force)*TIME_STEP;
				//remove edge ... set to special value
				if(damage[j] > 1.0){
					cout<<"Breaking bond between "
					<<edges[j*2]<<" and "<<edges[2*j +1]<<" F, s/L = "<<force \
					<<", "<<s/L[j]<<endl;
					edges[j*2] = -1; edges[j*2+1] = -1;
				}
			}
			else{
				damage[j] = s/L[j];
				if(damage[j] > 0.9){
					cout<<"Breaking bond between "
					<<edges[j*2]<<" and "<<edges[2*j +1]<<" F, s/L = "<<force \
					<<", "<<s/L[j]<<endl;
					edges[j*2] = -1; edges[j*2+1] = -1;
				}
			}
		}
	}
	return;

}