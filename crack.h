#ifndef CRACK_H
#define CRACK_H

#include <stddef.h>

using namespace std;

class Crack {

private:
	
	void clear();
	void copy(Crack const & source);

public:
	float * c;
	float * a;
	int dim;
	Crack(); //: Crack(2) {}
	Crack(int dim);
	~Crack();
	Crack(Crack const & source);
	Crack const & operator=(Crack const & other);

};

//#include "Crack.cpp"
#endif