#ifndef CRACK_H
#define CRACK_H

using namespace std;

class Crack {

private:
	
	void clear();
	void copy(Crack const & source);

public:
	double * c;
	double * a;
	int DIM;
	Crack() : Crack(2) {}
	Crack(int DIM);
	~Crack();
	Crack(Crack const & source);
	Crack const & operator=(Crack const & other);

};

#include "Crack.cpp"
#endif