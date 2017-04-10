#include "network.h"
#include "gnuplot_i.hpp"

Network::Network() {
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

void Network::build_network() {
	//TODO: Add small-network code
	
}


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

	if(IMPLEMENT_PBC){
		this->make_edge_connections(15.0);
		cout<<"Number after new connections made: "<<n_elems<<endl;}

}

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
		if (update_damage){
			if(RATE_DAMAGE){
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

}

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

void Network::apply_crack(Cracklist & alist) {
	std::default_random_engine seed;
	std::uniform_real_distribution<float> generator(0, 1);
	float equation = 0;
	//make clusters
	vector<int> within_circle;
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
			x = R[i*DIM];
			y = R[i*DIM+1];
			x -= crack.c[0];
			y -= crack.c[1];
			si = crack.trig[0];
			co = crack.trig[1];
			equation = pow(x*co + y*si, 2)/pow(crack.a[0],2) + \
						pow(x*si - y*co, 2)/pow(crack.a[1],2);

			equation -= 1.0;

			if(equation<=0.0){
				within_circle.push_back(i);
			}
		}
	}

	for(int i= 0; i < n_nodes; i++){
		if(contains(within_circle, i)){
			nodes_to_remove.push_back(i);
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



void Network::plotNetwork(int iter_step, bool first_time){
	ofstream f;
	Gnuplot gnu;
	float c;
	f.open("data.txt",std::ofstream::out | std::ofstream::trunc);
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
	stringstream convert;
	convert<<"[-100:"<<MAXBOUND+100<<"]";
	string xrange = convert.str();
	string bs = "/";
	string path = FLDR_STRING + bs + std::to_string(iter_step) + ".png";

	gnu.cmd("set xrange " + xrange);
	gnu.cmd("load 'viridis.pal'");
	gnu.cmd("set key off");
	gnu.cmd("set colorbox default vertical");
	gnu.cmd("set cbrange [0:1]");
	gnu.cmd("set cbtics ('0.0' 0.0,'1.0' 1.0) offset 0,0.5 scale 0");
	gnu.cmd("plot 'data.txt' every ::0::1 u 1:2:3 w l lc palette lw 2 title '"+\
		std::to_string(iter_step)+"'");
	gnu.cmd("set term png size 1200,900");
	gnu.cmd("set output '"+ path + "'");
	gnu.cmd("replot");
	gnu.cmd("set term x11");
	gnu.reset_plot();
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
	cout<<"Number of edges: "<<c<<endl;
	cout<<"Number of edges with x>0.99L: "<<short_L<<endl;
	cout<<"node-node distance mean: "<<mean_x<<endl;
	cout<<"node-node distance std_dev: "<<sqrt(var_x)<<endl;
	cout<<"node-node distance max: "<<max_x<<endl;
	cout<<"\n";
	cout<<"node-node x/L mean: "<<mean_t<<endl;
	cout<<"node-node x/L std_dev: "<<sqrt(var_t)<<endl;	
	cout<<"node-node x/L max: "<<max_t<<endl;

	// if(shortage/get_weight() > 0.1){
	// 	cout<<"Disorder too high for given mesh! Exiting...\n";
	// 	return true;
	// }
	if((float)c/(float)n_elems < 0.02){
		cout<<"Too few edges remain in given mesh! Exiting...\n";
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
	return alpha;
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

