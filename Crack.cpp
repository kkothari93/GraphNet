#include "crack.h"

// Crack::Crack() : Crack(2) {
	
// }

Crack::Crack() : Crack(2) {}

Crack::Crack(int DIM) {

	this->DIM = DIM;
	c = new double[DIM];
	a = new double[DIM];
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

	DIM = source.DIM;
	for (int i = 0; i < DIM; i++) {
		this->a[i] = source.a[i];
		this->c[i] = source.c[i];
	}
}

