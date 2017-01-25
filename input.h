#ifndef _input_included
#define _input_included

#include <fstream>
#include <sstream>
#include <string>
using namespace std;

int get_num_vertices(int);

void mapping(int&, int, int, stringstream& , int* );

void take_input(float*, int*, int&, int&);

#endif