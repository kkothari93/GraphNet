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

#include "mpi_network.h"
using namespace std;



int main() {

	//string path = "/media/konik/Research/2D sacrificial bonds polymers/cpp11_code_with_cuda/template2d.msh";
	string path = "./template2d_z4.msh";
	#if SACBONDS
	#define DECL_NET sacNetwork test_network(path)
	#else
	#define DECL_NET Network test_network(path)
	#endif
	DECL_NET;

	float weight_goal = 1.03754e6; // weight of similarly sized triangular mesh network
	
	float weight_multiplier;
	float weight = test_network.get_weight();
	if (weight<weight_goal){weight_multiplier = test_network.set_weight(weight_goal);}
	bool should_stop = test_network.get_stats();
	int old_n_edges = test_network.get_current_edges();
	int curr_n_edges = old_n_edges;


	if(should_stop){return 0;}
	float* plate_forces;
	plate_forces = (float*)malloc(sizeof(float)*DIM*STEPS);
	memset(plate_forces, 0.0, STEPS*DIM*sizeof(*plate_forces));

	if(CRACKED){
		Cracklist alist(4, MAXBOUND);
		test_network.apply_crack(alist);
	}


	test_network.plotNetwork(0, true);

	clock_t t = clock();
	cout<<"\n Will run for "<<STEPS<<":\n";
	
	for(int i = 0; i<STEPS; i++ ){
		
		test_network.optimize();
		test_network.move_top_plate();
		test_network.get_plate_forces(plate_forces, i);
		if((i+1)%100 == 0){
			should_stop = test_network.get_stats();
			if(should_stop){break;}
			curr_n_edges = test_network.get_current_edges();
			if(curr_n_edges<old_n_edges){
					test_network.plotNetwork(i, false);
			}
			cout<<"Step "<<(i+1)<<" took "<<float(clock()-t)/CLOCKS_PER_SEC<<" s\n";
			t = clock();  // reset clock
		}

	}
	string sb = SACBONDS ? "true" : "false" ; 
	string fname = FNAME_STRING + std::to_string(L_STD/L_MEAN) + "_" + sb + ".txt";
	write_to_file<float>(fname, plate_forces, STEPS, DIM);

	free(plate_forces);
	return 0;
}
