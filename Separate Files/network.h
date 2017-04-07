#ifndef NETWORK_H
#define NETWORK_H
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
#include "vel.h"

#include "helper_funcs.h"
#include "crack.h"
using namespace std;

const float PBC_vector[DIM] = {MAXBOUND*1.1, 0};
const float vel[DIM] = {vel_x, vel_y};

class Network {

public:
	Network();
	Network(Network const & source);
	Network(string& fname);
	virtual ~Network();
	Network const & operator=(Network const & other);
	virtual void build_network();
	void apply_crack(Cracklist &);
	virtual void load_network(string&);
	virtual void malloc_network(string&);
	void make_edge_connections(float dely_allowed = 10.0);
	virtual void get_forces(bool);
	virtual void move_top_plate();
	virtual void get_plate_forces(float*, int);
	virtual void optimize(float eta = 0.1, float alpha = 0.9, int max_iter = 800);
	float get_weight();
	float set_weight(float);
	bool get_stats();
	int get_current_edges();
	virtual void plotNetwork(int, bool);

//protected:

	virtual void clear();
	void copy(Network const & source);
	bool cracked;
	//int DIM;
	int n_nodes;
	int n_elems;
	int n_rside, n_lside, n_bside, n_tside;
	float * R; //n_nodes*DIM
	int * edges; //n_elems*2
	float * forces; //n_nodes*DIM
	float * damage;// n_nodes
	float * L; //n_elems
	bool* PBC; //n_elems
	int* lsideNodes; //max_nodes_on_a_side*2
	int* rsideNodes;
	int* tsideNodes;
	int* bsideNodes;
	bool initialized;
	//add moving nodes to speed up force
	int* moving_nodes;
	int n_moving;
};

#endif