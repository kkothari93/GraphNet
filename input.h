#ifndef _input_included
#define _input_included

#include <fstream>
#include <sstream>
#include <string>
using namespace std;

int get_num_vertices(int);

bool does_file_exist(string&);

void mapping(int&, int, int, stringstream& , int* );

void mapping(int&, int);

void take_input(float*, int*, int, int, string&);

void read_n(int&, int&, string&);

#endif