#define __params__	
#define TIME_STEP 1e-4						// Time step
#define SIM_TIME 10.0							// Std. deviation for contour lengths
#define Z_MAX 10							// Max coordination number
#define TOL 1e-6
#define MAXBOUND 500.0
#define L_MEAN 120.0						// Average for contour length
#define L_STD 4.0	

// Define constants
#define __constants__

#define kB 1.38064852e-5					// Boltzmann constant
#define b 0.1								// Persistence length
#define T 300 								// Temperature
#define ae 0.1 								// Strength of bond - includes activation energy
#define delxe 0.15 							// parameter for breaking crosslink connection
#define BLOCK_SIZE 1024

void side_nodes(double* R,\
	int* lnodes, int* rnodes, int* tnodes, int* bnodes,\
	int& n_lside, int& n_rside, int& n_tside, int& n_bside, int n, int DIM){
	// make the list of side nodes
	for(int i=0; i<n; i++){
		if(fabs(R[i*DIM]-0.0)<TOL){lnodes[n_lside]=i; n_lside++;}
		if(fabs(R[i*DIM]-MAXBOUND)<TOL){rnodes[n_rside]=i; n_rside++;}
		if(fabs(R[i*DIM + 1]-0.0)<TOL){bnodes[n_bside]=i; n_bside++;}
		if(fabs(R[i*DIM + 1]-MAXBOUND)<TOL){tnodes[n_tside]=i; n_tside++;}
	}
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

void take_input(double* R, int* edges, int& num_nodes, int& n_elems, int DIM) {
	
	string line;
	ifstream source;
	stringstream in_pos;
	source.open("/media/konik/New Volume/Poly Network/GMSH templates/template2d.msh");

	bool can_i_read_nodes= false,can_i_read_elems = false;
	bool read_num_nodes = false, read_n_elems = false;

	do{
		getline(source, line);
		//cout<<line<<"\n";
		//cout<<read_n_elems;
		if(line=="$Nodes"){
			getline(source, line);
			num_nodes=stoi(line);
			//cout<<"Number of nodes in the function "<<num_nodes<<"\n";
			double r[3]; int id;
			for(int i=0; i<num_nodes; i++){
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
		if(line=="$Elements"){
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
			n_elems = c;
			}
			//TODO: get elements
	}while(!read_n_elems);
	source.close();
}

void __init__(double* L, double* damage, bool* PBC, int n_elems){
	std::default_random_engine seed;
	std::normal_distribution<double> generator(L_MEAN, L_STD);
	for(int i=0; i<n_elems; i++){
		L[i] = generator(seed);
		damage[i] = 0.0;
		PBC[i] = false;
	}	
}

inline void normalize_vector(double* result, const double* vec, int DIM){
	double norm = getnorm(vec);
	for (int i = 0; i<DIM; i++){
		result[i] = vec[i] / norm;
	}
}

inline void normalize_vector(double* vec, int DIM){
	double norm = getnorm(vec);
	for (int i = 0; i<DIM; i++){
		vec[i] = vec[i] / norm;
	}
}

inline void unitvector(double* result, double* r1, double* r2, int DIM = 2){
	for (int j = 0; j<DIM; j++){
		result[j] = r1[j] - r2[j];
	}
	normalize_vector(result, DIM);
}

inline double force_wlc(double x, double L){
	double t = x / L;
	if (t < 0.99){ return kB*T / b * (t + 1.0 / 4.0 / pow((1 - t), 2) - 1.0 / 4.0); }
	else { return 999999.0; }
}

inline void convert_to_vector(double* result, const double mag, const double* direction, int DIM){
	for (int i = 0; i<DIM; i++){
		result[i] = mag*direction[i];
	}
}

void forcevector(double* result, double* r1, double* r2, double L, int DIM){
	double rhat[DIM];
	double s = dist(r1, r2);
	unitvector(rhat, r1, r2, DIM);
	double force = force_wlc(s, L);
	convert_to_vector(result, force, rhat, DIM);
}

inline double kfe(double force_mag){
	
	double toRet =  ae * exp(force_mag * delxe / kB / T);
	return toRet;
}



Network::Network() {

	DIM = 2;
	cracked = false;
	n_nodes = 1600;
	n_elems = 11000;
	n_rside = 0;
	n_lside = 0;
	n_bside = 0;
	n_tside = 0;
	R = new double[n_nodes * DIM];
	forces = new double[n_nodes * DIM];
	edges = new int[Z_MAX * n_nodes * 2];
	int max_nodes_on_a_side = int(sqrt(n_nodes))*2;
	lsideNodes = new int[max_nodes_on_a_side];
	rsideNodes = new int[max_nodes_on_a_side];
	tsideNodes = new int[max_nodes_on_a_side];
	bsideNodes = new int[max_nodes_on_a_side];
	damage = new double[2 * n_elems];
	L = new double[2 * n_elems];
	PBC = new bool[2 * n_elems];
	edge_matrix = new bool*[n_elems];
	for (int i = 0; i < n_elems; i++) {
		edge_matrix[i] = new bool[i+1];
		for (int j = 0; j <= i; i++) {
			edge_matrix[i][j] = false;
		}
	}

}

Network::Network(Network const & source) {

	copy(source);
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

	delete[] R;
	R = NULL;
	delete[] edges;
	edges = NULL;
	delete[] forces;
	forces = NULL;
	delete[] damage;
	damage = NULL;

	for (int i = 0; i < DIM; i++) {
		delete[] edge_matrix[i];
		edge_matrix[i] = NULL;
	}
	delete[] edge_matrix;
	edge_matrix = NULL;
	delete[] L;
	L = NULL;
	delete[] PBC;
	PBC = NULL;
	delete[] lsideNodes;
	lsideNodes = NULL;
	delete[] rsideNodes;
	rsideNodes = NULL;
	delete[] tsideNodes;
	tsideNodes = NULL;
	delete[] bsideNodes;
	bsideNodes = NULL;
}

void Network::build_network() {

	
	//num_edges = Z_MAX * n_nodes * 2;

	const string fname = "coordinates.txt";
	cout<<"Reading the mesh...\n";
	take_input(R, edges, n_nodes, n_elems, DIM);
	cout<<"Mesh read successfully!\n";
	cout<<"Number of nodes are: "<<n_nodes<<endl;
	cout<<"Number of elements are: "<<n_elems<<endl;
	int max_nodes_on_a_side = int(sqrt(n_nodes))*2;
	//int n_rside = 0, n_lside = 0, n_bside = 0, n_tside = 0;

	side_nodes(R, lsideNodes, rsideNodes, tsideNodes, bsideNodes, n_lside, n_rside, n_tside, n_bside, n_nodes, DIM);	

	__init__(L, damage, PBC, n_elems);

	const double PBC_vector[DIM] = {MAXBOUND*1.2, 0};

	this->make_edge_connections(15.0);
	cout<<"Number after new connections made: "<<n_elems<<endl;

	

	for (int i = 0; i < Z_MAX * n_nodes * 2; i += 2) {
		int node1 = edges[i];
		int node2 = edges[i+1];
		if(node1 == -1 || node2 == -1){
			continue;
		}
		else if (node1 < node2) {
			edge_matrix[node1][node2] = true;
		}

	}
	

}

void Network::load_network() {

}

void Network::copy(Network const & source) {

	cracked = source.cracked;
	DIM = source.DIM;
	n_nodes = source.n_nodes;
	n_elems = source.n_elems;
	n_rside = source.n_rside;
	n_lside = source.n_lside;
	n_bside = source.n_bside;
	n_tside = source.n_tside;
	R = new double[n_nodes * DIM];
	forces = new double[n_nodes * DIM];
	for (int i = 0; i < n_nodes * DIM; i++) {
		R[i] = source.R[i];
		forces[i] = source.forces[i];
	}
	edges = new int[Z_MAX * n_nodes * 2];
	for (int i = 0; i < Z_MAX * n_nodes * 2; i++) {
		edges[i] = source.edges[i];
	}
	int max_nodes_on_a_side = int(sqrt(n_nodes))*2;
	lsideNodes = new int[max_nodes_on_a_side];
	rsideNodes = new int[max_nodes_on_a_side];
	tsideNodes = new int[max_nodes_on_a_side];
	bsideNodes = new int[max_nodes_on_a_side];

	for (int i = 0; i < max_nodes_on_a_side; i++) {
		lsideNodes[i] = source.lsideNodes[i];
		rsideNodes[i] = source.rsideNodes[i];
		tsideNodes[i] = source.tsideNodes[i];
		bsideNodes[i] = source.bsideNodes[i];

	}

	damage = new double[2 * n_elems];
	L = new double[2 * n_elems];
	PBC = new bool[2 * n_elems];

	for (int i = 0; i < 2 * n_elems; i++) {
		damage[i] = source.damage[i];
		L[i] = source.L[i];
		PBC[i] = source.PBC[i];
	}

	edge_matrix = new bool*[n_elems];
	for (int i = 0; i < n_elems; i++) {
		edge_matrix[i] = new bool[i+1];
		for (int j = 0; j <= i; i++) {
			edge_matrix[i][j] = source.edge_matrix[i][j];
		}
	}



}

void Network::get_forces(const double* PBC_vector, bool update_damage = false) {

	int node1, node2;
	int j, k, id; // loop variables
	double r1[DIM]; double r2[DIM] ;
	double edge_force[DIM];

	for (j = 0; j < n_elems; j++){
		// read the two points that form the edge // 2 because 2 points make an edge! Duh.
		node1 = edges[j * 2]; 
		node2 = edges[j * 2 + 1];
		
		// check if pair exists
		if(node1 == -1 || node2 == -1) {
			continue;
		}

		// read the positions
		for(k = 0; k<DIM; k++){
			r1[k] = R[node1*DIM + k]; 
			r2[k] = R[node2*DIM + k];
		}
		
		// check PBC_STATUS
		if (PBC[j]) {
			// add PBC_vector to get new node position
			for (k = 0; k < DIM; k++){
				r2[k] += PBC_vector[k];
			}
			// get force on node1 due to node2
			forcevector(edge_force, r1, r2, L[j], DIM);
			// subtract back the PBC_vector to get original node position
			for (k = 0; k < DIM; k++){
				r2[k] -= PBC_vector[k];
			}
		}
		else{
			forcevector(edge_force, r1, r2, L[j], DIM);
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

}

void Network::make_edge_connections(double dely_allowed = 10.0) {

	std::default_random_engine seed;
	std::normal_distribution<double> generator(L_MEAN, L_STD);
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

	double equation = 0;
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


