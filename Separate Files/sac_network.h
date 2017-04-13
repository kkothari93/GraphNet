#ifndef SACNETWORK_H
#define SACNETWORK_H
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
#include "network.h"
using namespace std;

/**
@file sac_network.h
\brief Extends the network defintion to include hidden length
effects as described by Lieou et al., 2013 PRE 88, 012703 article

Sacrificial bonds and hidden length in biomaterials: A kinetic 
constitutive description of strength and toughness in bone
*/

class sacNetwork : virtual public Network{
	public:
		int* m; //n_elems
		float* sacdamage; // n_elems
		sacNetwork();
		~sacNetwork();
		sacNetwork(sacNetwork const &);
		sacNetwork(string& fname);
		void clear();
		void malloc_network(string&);
		void load_network(string&);
		virtual void get_forces(bool);

};

#endif