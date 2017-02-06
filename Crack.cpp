#include "crack.h"
// Crack::Crack() : Crack(2) {
	
// }

// Crack::Crack() :Crack(2) {}

Crack::Crack(float cx, float cy, float ax, float ay) {

	this->dim = 2;
	c = new float[dim]; 
	a = new float[dim];
	c[0] = cx;
	c[1] = cy;
	a[0] = ax;
	a[1] = ay;
}

Crack::Crack(float cx, float cy, float cz, float ax, float ay, float az) {
	this->dim = 3;
	c = new float[dim]; 
	a = new float[dim];
	c[0] = cx;
	c[1] = cy;
	c[2] = cz;
	a[0] = ax;
	a[1] = ay;
	a[2] = az;
}

Crack::~Crack() {
	clear();
}

Crack::Crack(Crack const & source) {
	copy(source);
}

Crack const & Crack::operator=(Crack const & other) {
	
	if (this != &other) {
		clear();
		copy(other);
	}
	return *this;	
}
	

void Crack::clear() {
	delete[] a;
	delete[] c;
	a = NULL;
	c = NULL;
}

void Crack::copy(Crack const & source) {

	dim = source.dim;
	for (int i = 0; i < dim; i++) {
		this->a[i] = source.a[i];
		this->c[i] = source.c[i];
	}
}

