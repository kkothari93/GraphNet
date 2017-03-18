#include "crack.h"
#include <math.h>
#include <random>
#define PI 3.14159265
using namespace std;
// Crack::Crack() : Crack(2) {
	
// }

// Crack::Crack() :Crack(2) {}

Cracklist::Cracklist(int n, float bound){
	this->dim = 2;
	n_cracks = n;
	listofCracks = new Crack[n];

	std::default_random_engine gen;
	//std::uniform_int_distribution<int> n_dist(2,4);
	std::uniform_real_distribution<float> a_dist(0.01,0.02);
	std::uniform_real_distribution<float> b_dist(0.05,0.06);
	std::uniform_real_distribution<float> angle(0.0,PI/2.0);
	std::uniform_real_distribution<float> cx_dist(0.2, 0.8);
	std::uniform_real_distribution<float> cy_dist(0.2, 0.8);
	//std::uniform_real_distribution<float> cz(-PI/2.0,PI/2.0);
	
	//int n = n_dist(gen);
	float ang = angle(gen);
	float a = a_dist(gen)*bound;
	float b = b_dist(gen)*bound;
	float cx = cx_dist(gen)*bound;
	float cy = cy_dist(gen)*bound;
	float s = sin(ang);
	float c = cos(ang);

	for(int i=0; i < n_cracks; i++){
		listofCracks[i].setter(cx,cy,a,b,s,c);
		//printf("%d\n",__LINE__);
		listofCracks[i].print_info();
		ang = angle(gen);
		a = a_dist(gen)*bound;
		b = b_dist(gen)*bound;
		cx = cx_dist(gen)*bound;
		cy = cx_dist(gen)*bound;
		s = sin(ang);
		c = cos(ang);
	}

}

Cracklist::~Cracklist(){
	delete[] listofCracks;
	listofCracks = NULL;
}

void Crack::print_info(){
	printf("\n");
	printf("Crack center at: (%0.2f, %0.2f)\n", c[0], c[1]);
	printf("Crack axes : (a: %0.2f, b: %0.2f)\n", a[0], a[1]);
	printf("Crack angle : %0.2f\n", asin(trig[0])*180.0/PI);
	printf("----------------------\n");
}


Crack::Crack(){
	this->dim = 2;
	c = new float[dim]; 
	a = new float[dim];
	trig = new float[dim];
	c[0] = 0;
	c[1] = 0;
	a[0] = 0;
	a[1] = 0;
	trig[0] = 0;
	trig[1] = 0;
}

void Crack::setter(float cx, float cy, float ax, float ay, float sine, float cosine) {
	c[0] = cx;
	c[1] = cy;
	a[0] = ax;
	a[1] = ay;
	trig[0] = sine;
	trig[1] = cosine;
}

void Crack::setter(Crack & source){
	c[0] = source.c[0];
	c[1] = source.c[1];
	a[0] = source.a[0];
	a[1] = source.a[1];
	trig[0] = source.trig[0];
	trig[1] = source.trig[1];
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


Crack & Cracklist::operator[](int b) {
	return listofCracks[b];	
}	

void Crack::clear() {
	delete[] a;
	delete[] c;
	delete[] trig;
	a = NULL;
	c = NULL;
	trig = NULL;
}

void Crack::copy(Crack const & source) {

	dim = source.dim;
	for (int i = 0; i < dim; i++) {
		this->a[i] = source.a[i];
		this->c[i] = source.c[i];
		this->trig[i] = source.trig[i];
	}
}

// int main(){
	
// 	Cracklist a(5,100);
// 	return 0;

// }