/**
@file network.cpp
This file implements all functions defined in header network.h
*/
#include "network.h"
#include "gnuplot_i.hpp"

// ----------------------------------------------------------------------- 
/// \brief Default constructor
// ----------------------------------------------------------------------- 
Network::Network() {

	initialized = false;

}

// ----------------------------------------------------------------------- 
/// \brief Copy constuctor for class Network
///
/// \param source --> copies from a Network object source
// ----------------------------------------------------------------------- 
Network::Network(Network const & source) {

	initialized = false;
	copy(source);
	initialized = true;

}

// ----------------------------------------------------------------------- 
/// \brief Constructor that instantiates a Network object according
/// to given \emph{fname} file. Please note that there are no checks made
/// on the GMSH .msh filepath and that responsibility is left to the user.
///
/// \param fname --> takes a filename and instantiates Network object
///
// ----------------------------------------------------------------------- 
Network::Network(string& fname, bool from_dump) {

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
/// \brief Destructor for class network
///
// ----------------------------------------------------------------------- 
Network::~Network() {

	clear();

}


// ----------------------------------------------------------------------- 
/// \brief Overloads assignment operator (=) to implement copy constructor
///
// ----------------------------------------------------------------------- 
Network const & Network::operator=(Network const & other) {

	if (this != &other) {
		clear();
		copy(other);
	}

	return *this;
	
}

// ----------------------------------------------------------------------- 
/// \brief Called by destructor of Network objects to clear variables 
/// from memory
///
// ----------------------------------------------------------------------- 
void Network::clear() {

	free(R);
	R = NULL;
	free(edges);
	edges = NULL;
	free(forces);
	forces = NULL;
	free(damage);
	damage = NULL;
	free(L);
	L = NULL;
	free(PBC);
	PBC = NULL;
	free(lsideNodes);
	lsideNodes = NULL;
	free(rsideNodes);
	rsideNodes = NULL;
	free(tsideNodes);
	tsideNodes = NULL;
	free(bsideNodes);
	bsideNodes = NULL;
	free(moving_nodes);
	moving_nodes = NULL;

}

// ----------------------------------------------------------------------- 
/// \brief This function is for users to build custom networks. It has no 
/// body and implements nothing. 
///
// ----------------------------------------------------------------------- 
void Network::build_network() {
	//TODO: Add small-network code
	
}

// ----------------------------------------------------------------------- 
/// \brief This function allocates memory for all class data members. 
/// This is called by the Network constructor if the same is called using
/// a filepath string. User is discouraged from directly calling this function.
/// 
///  
// -----------------------------------------------------------------------
void Network::malloc_network(){

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

	// 	initialise n_xside for side nodes
	n_rside = 0;
	n_lside = 0;
	n_bside = 0;
	n_tside = 0;

	// adjust n_elems back to actual
	n_elems -= 3*max_nodes_on_a_side;
}



// ----------------------------------------------------------------------- 
/// \brief This function loads the last network state from a dump file.
/// This is called by the Network constructor load_from_dump is set to true.
/// Do not load from dump if network has already broken.
/// 
/// \param dname (string) --> dump filename
/// 
// -----------------------------------------------------------------------
void Network::load_from_dump(string& dname){
	
	string last_iter = read_dump(n_nodes, n_elems, dname);
	iter_offset = stoi(last_iter);

	cout<<"Restarting from iter "<<iter_offset<<endl;


	malloc_network();

	//read file
	ifstream dump;
	string line;
	dump.open(dname);

	stringstream out;
	//read R positions
	while(!dump.eof()){
		getline(dump, line);
		if(line==last_iter){
			break;	
		}
	}

	// added to skip over line "START NODE POSITIONS"
	getline(dump, line);

	// read R positions
	int index;
	float x, y;
	for(int i=0; i<n_nodes; i++){
		getline(dump, line);
		out<<line;
		out>>index;
		out>>x>>y;
		out.str(std::string());
		out.clear();
		//cout<<index<<", "<<x<<", "<<y<<endl;
		R[DIM*index] = x;
		R[DIM*index+ 1] = y;
	}
	// --finished reading R

	while(!dump.eof()){
		getline(dump, line);
		if(line.find("START ACTIVE EDGES")!=string::npos){
			break;	
		}
	}


	// read edges
	int node1, node2, p;
	float d, s;
	for(int i=0; i<n_elems; i++){
		getline(dump, line);
		
		out<<line;
		out>>index>>node1>>node2>>d>>p;
		out.str(std::string());
		out.clear();

		edges[2*i] = node1;
		edges[2*i + 1] = node2;
		if(node1 >= 0 && node2 >= 0){	
			s = dist(&R[node1*DIM], &R[node2*DIM]);

			if(d!=0){
				L[i] = s/d;
			}
			else{
				cout<<"could not find L for "<<node1<<", "<<node2<<endl;
				L[i] = L_MEAN; // should not be the case!
			}
		}		
	
		damage[i] = d;

		PBC[i] = (p==0) ? false : true;
	}
	// --finished reading edges

	dump.close();

	float max_y = 0;

	// This will find the max_y for the loaded iter from dump file
	// The max_y is the position for the top plate
	for(int i = 0; i<n_nodes; i++){
		y = R[DIM*i + 1];
		if(y>max_y){
			max_y = y;
		}
	}

	// side_nodes
	side_nodes(0.0, max_y);

	// find moving nodes
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
}



// ----------------------------------------------------------------------- 
/// \brief Gives the boundary nodes for a network.
///
// -----------------------------------------------------------------------
void Network::side_nodes(float max_x, float max_y){

	// make the list of side nodes
	for(int i=0; i<n_nodes; i++){
		if(fabs(R[i*DIM]-0.0)<TOL){n_lside++;}
		if(fabs(R[i*DIM]-max_x)<TOL){n_rside++;}
		if(fabs(R[i*DIM + 1]-0.0)<TOL){n_bside++;}
		if(fabs(R[i*DIM + 1]-max_y)<TOL){n_tside++;}
	}

	size_t si = sizeof(int);
	lsideNodes = (int* )malloc(n_lside*si);
	rsideNodes = (int* )malloc(n_rside*si);
	bsideNodes = (int* )malloc(n_bside*si);
	tsideNodes = (int* )malloc(n_tside*si);
	
	n_lside = 0; n_rside = 0; n_bside = 0; n_tside=0;

	for(int i=0; i<n_nodes; i++){
		if(fabs(R[i*DIM]-0.0)<TOL){lsideNodes[n_lside]=i; n_lside++;}
		if(fabs(R[i*DIM]-max_x)<TOL){rsideNodes[n_rside]=i; n_rside++;}
		if(fabs(R[i*DIM + 1]-0.0)<TOL){bsideNodes[n_bside]=i; n_bside++;}
		if(fabs(R[i*DIM + 1]-max_y)<TOL){tsideNodes[n_tside]=i; n_tside++;}
	}
}



// ----------------------------------------------------------------------- 
/// \brief Adds y-directional long range edges to the network. Chooses two nodes
/// uniformly at random from the n_add divisions in x_direction such that the 
/// distance between them is at least 0.25*L_MEAN and at max L_MEAN
///
/// \param n_add number of long range y-connections to add
/// \param prestrech float between 0.01,0.99 to give initial x/L of these long range bonds
/// If value is not between (0.01,0.99), 0.25 is assumed 
///
// -----------------------------------------------------------------------
void Network::add_long_range_egdes_y(int n_add, float prestretch){
	// Add structured edges in y direction. To do so, first create 2*n_add circles
	// then find a node in the circle radius (chosen here to be L_STD) and add an edge
	// between these two nodes
	float* centers = new float[n_add*2*2]; // 2 times n_add circles, each taking 
										   // two floats to store centers
	
	cout<<"Asked to add "<<n_add<<" more edges!\n";

	if(n_add >= n_elems/2){
		cout<<"You're asking for too many long range bonds,"
				" can't do it! I will add none and continue..."<<endl;
	}
	
	if(prestretch < 0.01 || prestretch > 0.99){
		cout<< "Value for prestretch needs to be between 0.01 and 0.99\n";
		prestretch = 0.25;
	}

	// fill the centers
	for(int i = 0 ; i <n_add*2 ; i++){
		if(i<n_add){
			centers[2*i] = MAXBOUND_X/2/n_add * (2*i+1);
			centers[2*i + 1] = MAXBOUND_Y/4;
		}
		if(i>=n_add){
			centers[2*i] = MAXBOUND_X/2/n_add * (2*(i - n_add)+1);
			centers[2*i + 1] = 3*MAXBOUND_Y/4;	
		}
		cout<<centers[2*i]<<"\t"<<centers[2*i+1]<<endl;
	}


	// create circles for these centers
	Crack* circles = new Crack[2*n_add]();
	for(int i=0; i<2*n_add; i++){
		circles[i].setter(centers[2*i], centers[2*i+1], MAXBOUND_Y/10, MAXBOUND_Y/10, 1, 0);
	}

	// fill one id per circle
	int* node_per_circle = new int[2*n_add](); // value initialises to zero: ISO C++03 5.3.4[expr.new]/15
	for(int i = 0; i < 2*n_add; i++){
		node_per_circle[i] = 0;
	}
	

	for(int i = 4; i<n_nodes; i++){
		for(int k =0; k<2*n_add; k++){
			if(node_per_circle[k] != 0){
				continue;
			}
			if(circles[k].inside(&R[DIM*i])<=0){
				node_per_circle[k] = i;
			}
		}
	}

	float s; // will store distances
	//Add edges between selected nodes
	for(int i = 0; i < n_add; i++){
	
		edges[2*n_elems] = node_per_circle[i];
		edges[2*n_elems + 1] = node_per_circle[i + n_add];

		// Calculate distance
		s = dist(&R[node_per_circle[i]*DIM], &R[node_per_circle[i + n_add]*DIM]);

		// update the contour lengths according to position
		L[n_elems] = s/prestretch;

		// update PBC;
		PBC[n_elems] =  false;

		// update damage
		damage[n_elems] = 0;

		// increment n_elems
		n_elems++;
	}

	delete[] centers;
	centers = NULL;
	delete[] node_per_circle;
	node_per_circle = NULL;
	delete[] circles;
	circles = NULL;
}


// ----------------------------------------------------------------------- 
/// \brief Adds random long range edges to the network. Chooses two nodes
/// uniformly at random and adds an edge with x/L = 0.25 (i.e. prestressed)
///
// -----------------------------------------------------------------------
void Network::add_long_range_egdes_random(int n_add, float prestretch){
	
	cout<<"Asked to add "<<n_add<<" more edges!\n";
	int node1, node2;
	float s = 0;

	if(n_add >= n_elems/2){
		cout<<"You're asking for too many long range bonds,"
				" can't do it! I will add none and continue..."<<endl;
	}

	if(prestretch < 0.01 || prestretch > 0.99){
		cout<< "Value for prestretch needs to be between 0.01 and 0.99\n";
		prestretch = 0.25;
	}

	
	srand(time(NULL));

	for(int i = 0; i < n_add; i++){
		while(s<L_MEAN/4 || s>L_MEAN){ 
			// L_MEAN/4, L_MEAN is arbitrary threshold to keep the edge long range

			// random number between [min; max] 
			// min = 4 (all corner nodes); max = n_nodes - 1
			node1 = rand()%(n_nodes - 4) + 4;
			node2 = rand()%(n_nodes - 4) + 4;
		
			s = dist(&R[node1*DIM], &R[node2*DIM]);
		}


		cout<<"Adding edge between "<<node1<<", "<<node2<<endl;

		// add nodes to edges array
		edges[n_elems*2] = node1;
		edges[n_elems*2 + 1] = node2;

		// update the contour lengths according to position
		L[n_elems] = s/prestretch;

		// update PBC;
		PBC[n_elems] =  false;

		// update damage
		damage[n_elems] = 0;

		// increment n_elems
		n_elems++;

		s = 0;
	}

}


// ----------------------------------------------------------------------- 
/// \brief Removes duplicate edges. A much more efficient version can be made
/// of this program. But it is called only once in the simulation, so this 
/// should be fine.
// -----------------------------------------------------------------------
void Network::remove_duplicates(int& n_elems){
	int counter = 0;
	int node1, node2;
	int d_node1, d_node2, i,k;

	for(i=0;i<n_elems;i++){
		// here n_elems is needed instead of n_elems-1 for correct counting in counter
		node1 = edges[2*i];
		node2 = edges[2*i+1];
		if(node1!=-2 && node2!=-2){
			for(int j=i+1; j<n_elems; j++){
				d_node1 = edges[2*j];
				d_node2 = edges[2*j+1];
				if((node1==d_node1 && node2==d_node2) || (node2==d_node1 && node1==d_node2)){
					edges[2*j] = -2; // -1 token reserved for broken edges. -2 used for readable code
					edges[2*j+1] = -2; // -1 token reserved for broken edges. -2 used for readable code
				}
			}
			counter++;
		}
	}
	
	int* buffer = new int[2*counter];
	for(i=0, k=0;i<n_elems;i++){
		node1 = edges[2*i];
		node2 = edges[2*i+1];
		if(node1!=-2 && node2!=-2){
			buffer[2*k] = node1;
			buffer[2*k+1] = node2;
			k++;
		}
	}
	
	for(i=0, k=0;i<n_elems;i++){
		if(k<counter){
			edges[2*i] = buffer[2*k];
			edges[2*i+1] = buffer[2*k+1];
			k++;
		}
		else{
			edges[2*i] = -1;
			edges[2*i+1] = -1;
		}
	}

	cout<<n_elems-counter<<" duplicates removed!\n";
	n_elems = counter;
	delete[] buffer;
	buffer = NULL;
}

// ----------------------------------------------------------------------- 
/// \brief This function builds the network. It includes building PBC 
/// bonds. The function initializes all properties required for running
/// experiments on the instantiated Network code. This is called by constructor
/// of a Network object. User is discouraged from directly calling this function.
///
/// \param fname (string) --> string of the .msh filepath
// -----------------------------------------------------------------------
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

	// read n_nodes, n_elems
	read_n(n_nodes, n_elems, fname);
	
	//Malloc all variables
	malloc_network();

	cout<<"Malloc was successful!\n";

	cout<<"Reading the mesh...\n";
	take_input(R, edges, n_nodes, n_elems, fname);

	// remove duplicate edges
	remove_duplicates(n_elems);

	// add long range connections
	// add_long_range_egdes_y(int(log(n_elems)), 0.25);

	cout<<"Mesh read successfully!\n";
	cout<<"Number of nodes are: "<<n_nodes<<endl;
	cout<<"Number of elements are: "<<n_elems<<endl;


	// initialize boundary node arrays
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
	__init__(L, damage, PBC, n_elems);

	if(IMPLEMENT_PBC){
		this->make_edge_connections(15.0);
		cout<<"Number after new connections made: "<<n_elems<<endl;
	}
}

// ----------------------------------------------------------------------- 
/// \brief Implements the body of the copy constuctor. Allocates memory and
/// does the appropriate assignments
///
/// \param source (Network) --> an initialized object of type Network
// -----------------------------------------------------------------------
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


// ----------------------------------------------------------------------- 
/// \brief Calculates forces along edges of the polymer network graph (i.e
/// forces in each polymer). The forces are then assigned to nodes which 
/// prepares the network for the optimization step. If update_damage is 
/// true, the appropriate damage metric is also updated. If damage exceeds
/// predefined thresholds, the edge is broken. In case of a sacNetwork object
/// hidden lengths are also opened according to damage in these hidden bonds
/// apart from end damage. The function does not return anything as it updates
/// the objects forces directly.
///
/// \param update_damage (bool) --> flag to update damage
// -----------------------------------------------------------------------
void Network::get_forces(bool update_damage = false) {

	int node1, node2;
	int j, k, id; // loop variables
	float r1[DIM]; float r2[DIM] ;
	float edge_force[DIM];
	float rhat[DIM];
	float s;
	float force;


	memset(forces, 0.0, n_nodes*DIM*sizeof(float));

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
			//if(force == 999999){edges[j*2] = -1; edges[j*2 +1] = -1; force =0.0;}
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
			//if(force == 999999){edges[j*2] = -1; edges[j*2 +1] = -1; force =0.0;}
			convert_to_vector(edge_force, force, rhat);
		}
		#pragma unroll
		for (k = 0; k < DIM; k++){
			forces[node1*DIM + k] -= edge_force[k];
			forces[node2*DIM + k] += edge_force[k];
		}

		//update damage if needed
		if (update_damage){
			if(RATE_DAMAGE){
				damage[j] += kfe(force)*TIME_STEP;
				//remove edge ... set to special value
				if(damage[j] > 1.0){
					cout<<"Breaking bond between "
					<<edges[j*2]<<" and "<<edges[2*j +1]<<" F, s/L = "<<force \
					<<", "<<s<<"/"<<L[j]<<endl;
					edges[j*2] = -1; edges[j*2+1] = -1;
				}
			}
			else{
				damage[j] = s/L[j]/0.9;
				if(damage[j] > 1.0){
					cout<<"Breaking bond between "
					<<edges[j*2]<<" and "<<edges[2*j +1]<<" F, s/L = "<<force \
					<<", "<<s<<"/"<<L[j]<<endl;
					edges[j*2] = -1; edges[j*2+1] = -1;
				}
			}
		}
	}

}


// ----------------------------------------------------------------------- 
/// \brief Adds PBC bonds to the network according to dely_allowed.
///
/// \param dely_allowed(float) --> Connects two nodes i and j on the two 
/// sides of the network respectively iff $abs(y_i - y_j) < $dely_allowed
// -----------------------------------------------------------------------
void Network::make_edge_connections(float dely_allowed) {

	std::default_random_engine seed;
	std::uniform_real_distribution<float> generator(L_MEAN - L_STD, L_MEAN + L_STD);
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


// ----------------------------------------------------------------------- 
/// \brief Adds random cracks in the Network object. This is implemented 
/// by removing all connections and nodes that lie \emph{entirely} within
/// the ellipses specified in alist. Please note this assumption may be
/// relaxed and a probabilistic measure may be introduced to remove some
/// connections so as to simulate small world networks (cluster networks)
/// This is managed by the PROB_REMOVAL parameter in params.h 
///
/// \param alist (Cracklist) --> A list of ellipses
// -----------------------------------------------------------------------
void Network::apply_crack(Cracklist & alist) {

	cracked = true;
	std::default_random_engine seed;
	std::uniform_real_distribution<float> generator(0, 1);
	float val = 0;
	//make clusters
	//vector<int> within_circle;
	vector<int> nodes_to_remove;
	int edges_removed = 0; 
	int node1, node2;
	Crack crack;
	float x;
	float y;
	float si,co;
	for(int ci = 0; ci < alist.n_cracks; ci++){
		crack.setter(alist[ci]);
		for(int i=0; i<n_nodes; i++){
			val = crack.inside(&R[i*DIM]);
			if(val<=0.0){
				nodes_to_remove.push_back(i);
			}
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


// ----------------------------------------------------------------------- 
/// \brief Plots the graph of the network. Note the edges in the graph just
/// signify connection between two nodes and do not represent the actual 
/// polymer geometry in any sense. Note that this plot will be saved as .png 
/// file in the folder specified by FLDR_STRING in params.h
///
/// \param iter_step (int) --> The iteration number of the experiment. This
/// number will also be the filename for the .png file
///
/// \param first_time (bool) --> This just effects what quantity is mapped
/// to the color of the graph. If True, colors represent magnitude of forces
/// in the polymer else, they represent the damage in a polymer connection
// -----------------------------------------------------------------------
void Network::plotNetwork(int iter_step, bool first_time){
	ofstream f;
	Gnuplot gnu;
	string fname = "data" + std::string(FLDR_STRING) + ".txt";
	float c;
	f.open(fname,std::ofstream::out | std::ofstream::trunc);
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
	float aspect_ratio = (MAXBOUND_X)/(R[tsideNodes[0]*DIM+1]); 

	stringstream convert;
	convert<<"[-100:"<<MAXBOUND_X+100<<"]";
	string xrange = convert.str();
	convert.str(std::string());
	
	int x_res = 1280;
	convert<<x_res<<", "<<int(x_res/aspect_ratio);
	string png_size = convert.str();
	convert.str(std::string());

	string bs = "/";
	string path = FLDR_STRING + bs + std::to_string(iter_step) + ".png";
	
	gnu.cmd("set xrange " + xrange);
	gnu.cmd("load 'viridis.pal'");
	gnu.cmd("set key off");
	gnu.cmd("set colorbox default vertical");
	gnu.cmd("set cbrange [0:1]");
	gnu.cmd("set cbtics ('0.0' 0.0,'1.0' 1.0) offset 0,0.5 scale 0");
	gnu.cmd("set size ratio -1"); // to have the same scale on x and y axes
	gnu.cmd("plot '" + fname + "' every ::0::1 u 1:2:3 w l lc palette lw 2 title '"+\
		std::to_string(iter_step)+"'");
	gnu.cmd("set term png size " + png_size);
	gnu.cmd("set output '"+ path + "'");
	gnu.cmd("replot");
	gnu.cmd("set term x11");
	gnu.reset_plot();
	
}

// ----------------------------------------------------------------------- 
/// \brief Gets the active (load-bearing) edges in the Network object. 
/// Also prints out that info to STDOUT or piped OUT
///
/// 
// -----------------------------------------------------------------------
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


// ----------------------------------------------------------------------- 
/// \brief Prints out certain statistics of the network object.
///
/// /param curr_force --> takes in curr_force on plate, returns True,
/// 					  if curr_force <= 0.1*__init_force
// -----------------------------------------------------------------------
bool Network::get_stats(float curr_force){
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

	if(n_elems <= 0){
		cout<<"Oops! Something's wrong! n_elems is not positive!"<<endl;
	}
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
	if(c<=0){
		cout<<"No load bearing edges remain! Stopping simulation!\n";
		return true;
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


	if((float)c/(float)n_elems < 0.05){
		cout<<"Too few edges remain in given mesh! Exiting...\n";
		return true;
	}
	if(fabs(curr_force) <= 0.1*fabs(__init_force)){
		cout<<"Force too low! Exiting...\n";
		return true;
	}
	return false;
}


// ----------------------------------------------------------------------- 
/// \brief Gets the total length of all polymer chains active in the network
/// 
// -----------------------------------------------------------------------
float Network::get_weight(){
	float sum_length = 0.0;
	for(int i = 0; i < n_elems; i++){
		if(edges[i*2] != -1 && edges[i*2 + 1]!=-1){
			sum_length += L_MEAN;
		}
	}
	cout<<"\nNetwork weight is "<<sum_length<<"\n";
	return sum_length;
}

// ----------------------------------------------------------------------- 
/// \brief Sets the total length of all polymer chains according to supplied
/// weight parameter. Also prints out the new average polymer contour lengths
/// 
/// /param weight (float) --> Weight goal required.
// -----------------------------------------------------------------------
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
	weight_multiplier = alpha;
	return alpha;
}


// ----------------------------------------------------------------------- 
/// \brief Moves the top plate according to velocity supplied in vel.h
/// 
///
// -----------------------------------------------------------------------
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


// ----------------------------------------------------------------------- 
/// \brief Equilibriates the Network object using an Implicit scheme based
/// on solving the quasi-dynamic equation
///
/// /param C --> damping coefficient
/// /param max_iter --> maximum allowable iterations in the implicit scheme
// -----------------------------------------------------------------------
void Network::qd_optimize(float C, int max_iter){
	float* Rp1 = new float[n_moving*DIM]();
	float g;
	int id, d, node;

	for(int step =0; step<max_iter; step++){
		// Predictor step
		get_forces(false);
		for(id = 0; id < n_moving; id++){
			node = moving_nodes[id];
			#pragma unroll
			for(d = 0; d<DIM; d++){
				Rp1[id*DIM + d] = R[node*DIM + d];
				g = forces[DIM*node+d];
				R[node*DIM + d] += g/C * TIME_STEP;
			}		
		}
		// Corrector step
		get_forces(false);
		for(id = 0; id < n_moving; id++){
			node = moving_nodes[id];
			#pragma unroll
			for(d = 0; d<DIM; d++){
				g = forces[DIM*node+d];
				R[node*DIM + d] = Rp1[id*DIM + d] + g/C * TIME_STEP;
			}		
		}
	}
	// update damage
	get_forces(true);
	delete[] Rp1;
}

// ----------------------------------------------------------------------- 
/// \brief Equilibriates the Network object using RMSProp algorithm due to 
/// G. Hinton et al (a gradient based force minimization over each node). 
/// Stops when maximum force on node is less than TOL specified in params.h
/// or max_iter is reached
/// 
/// /param eta (float) --> learning rate
/// /param alpha (float) --> history weighting parameter
/// /param max_iter (int) --> maximum iterations to be allowed in descent
// -----------------------------------------------------------------------
void Network::optimize(float eta, float alpha, int max_iter){
	float* rms_history = new float[n_moving*DIM](); // () allows 0.0 initialization
	float g, delR;
	char p;
	int id, d, node;
	for(int step = 0; step < max_iter; step++){
		get_forces(false);
		// plotNetwork(step, true);
		// cin>>p;
		if(getabsmax(forces,n_nodes*DIM)>TOL){
			for(id = 0; id < n_moving; id++){
				node = moving_nodes[id];
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
	//
	get_forces(true);
	delete[] rms_history;
}

// ----------------------------------------------------------------------- 
/// \brief Writes the forces on the top plate the supplied plate_forces array
/// 
/// /param plate_forces (float*) --> storage array for plate_forces
/// /param iter (int) --> force_{i} stored in plate_forces[iter*DIM + i]
/// i = 0 ... DIM - 1
/// /param max_iter (int) --> maximum iterations to be allowed in descent
// -----------------------------------------------------------------------
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
	if(iter==0){
		__init_force = plate_forces[1]; // initial force
	}
}

// ----------------------------------------------------------------------- 
/// \brief Dumps the network information for an iteration in a file. 
/// Filename would be <FLDR_STRING>.txt with <FLDR_STRING taken from params.h
/// If <FLDR_STRING>.txt exists, <FLDR_STRING>_i.txt would be created where i
/// is an integer decided based on copies present in current folder. 
/// 
/// /param iter The iteration number for the simulation
// -----------------------------------------------------------------------
void Network::dump(int iter, bool first_time){
	ofstream logger;
	std::time_t result = std::time(nullptr);
	string fname = FLDR_STRING;
	if(first_time){
		int c = filename(fname);
		if(c!=-1){
			fname = fname + "_" + std::to_string(c);
		}
	}
	fname = fname+".txt";
	logger.open(fname, ios_base::app);

	if(first_time){
		// If opened for the first time write the params
		logger<<"Network snapshots"<<iter<<"\n";
		logger<<"File created at "<<std::asctime(std::localtime(&result));
		logger<<"PARAMS:\n";

		logger<<"Sim dimension : "<<DIM<<"\n";
		logger<<"Simulation time : "<<SIM_TIME<<"\n";
		logger<<"Simulation time-step : "<<TIME_STEP<<"\n";
		logger<<"Velocity : "<<vel_x<<"\t"<<vel_y<<"\n";
		logger<<"MAXBOUND : "<<MAXBOUND_X<<"\t"<<MAXBOUND_Y<<"\n";

		logger<<"Disorder characteristics : "<<"\n";
		logger<<" -- L_MEAN : "<<L_MEAN<<"\n";
		logger<<" -- L_STD : "<<L_STD<<"\n";
		logger<<"Others : SACBONDS = "<<SACBONDS<<"; IMPLEMENT_PBC = "<<IMPLEMENT_PBC<<"\n";
		logger<<"Rate damage : "<<RATE_DAMAGE<<"\n";		
		logger<<"Cracked? : "<<CRACKED<<"; "<<PROB_REMOVAL<<"\n";
		logger<<"-------------------------------------------------\n\n";
		
		logger<<"CONSTANTS:\n";
		logger<<"Persistence length: "<<b_poly<<"\n";
		logger<<"Temperature: "<<T<<"\n";
		logger<<"End bond a: "<<ae<<"\n";
		logger<<"End bond delx: "<<delxe<<"\n";
		logger<<"Hidden bond a: "<<af<<"\n";
		logger<<"Hidden bond delx: "<<delxf<<"\n";
		logger<<"-------------------------------------------------\n\n";

	}
	if(!first_time){   
		logger<<"INFO FOR ITER\n"<<iter<<endl;
	    // Write R
	    logger<<"START NODE POSITIONS\n";
		for(int i = 0; i < n_nodes; i++){
			logger<<i<<"\t";
			for(int j = 0; j< DIM; j++){
				logger<<R[i*DIM+j]<<"\t";
			}
			logger<<"\n";
		}
		logger<<"END NODE POSITIONS\n";

		//write edges | damage | PBC_status 
		logger<<"START ACTIVE EDGES\n";
		for(int i = 0; i < n_elems; i++){
			if(edges[2*i]!=-1 && edges[2*i+1]!=-1){
				//Note i is then the original edge index
				logger<<i<<"\t"<<edges[2*i]<<"\t"<<edges[2*i+1]<<"\t"<<damage[i]<<"\t"<<PBC[i]<<"\n";
			}
		}
		logger<<"END ACTIVE EDGES\n";

		//write forces on nodes
		logger<<"START NODE FORCES\n";
		for(int i = 0; i < n_nodes; i++){
			logger<<i<<"\t";
			for(int j = 0; j< DIM; j++){
				logger<<forces[i*DIM+j]<<"\t";
			}
			logger<<"\n";
		}
		logger<<"END NODE FORCES\n\n";
	}

    logger.close();
	cout<<"Dumped everything in "<<fname<<"!\n";
}
