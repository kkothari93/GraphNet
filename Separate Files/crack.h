#ifndef CRACK_H
#define CRACK_H

#include <stddef.h>

using namespace std;

class Crack {

protected:
	
	void clear();
	void copy(Crack const & source);

public:
	float * c;
	float * a;
	float * trig;
	int dim;
	Crack(); //: Crack(2) {}
	void setter(float, float,float, float,float, float);
	void setter(Crack &);
	virtual ~Crack();
	Crack(Crack const & source);
	Crack const & operator=(Crack const & other);
	void print_info();

};

class Cracklist : public Crack{

public:
	int n_cracks;
	Crack* listofCracks;
	Crack & operator[](int);
	Cracklist(int, float);
	Cracklist();
	~Cracklist();
} ;

//#include "Crack.cpp"
#endif