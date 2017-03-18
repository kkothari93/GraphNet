#ifndef MPI_NETWORK_H
#define MPI_NETWORK_H
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
using namespace std;

class MPI_Network : public Network {

public:

	void get_forces(bool);
	void optimize(float eta = 0.1, float alpha = 0.9, int max_iter = 800);
	


};

#endif
