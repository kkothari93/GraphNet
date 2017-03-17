#ifndef SACNETWORK_H
#define SACNETWORK_H
#include <iostream>
#include <cmath>
#include "Crack.cpp"
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
#include "gnuplot_i.hpp"
#include "helper_funcs.cpp"
#include "network.h"
using namespace std;

class sacNetwork : public Network{
	private:
		int* m; //n_elems
		float* sacdamage; // n_elems
	public:
		sacNetwork();
		~sacNetwork();
		sacNetwork(string& fname);
		void clear();
		void malloc_network(string&);
		void load_network(string&);
		void get_forces(bool);


};

#endif