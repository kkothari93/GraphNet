/**
@file crack.cpp
\brief Implements all functions prototyped in crack.h
*/

#include "crack.h"
#include <math.h>
#include <random>
#define PI 3.14159265
using namespace std;

// ----------------------------------------------------------------------- 
//		Crack member functions 
//
// -----------------------------------------------------------------------


// ----------------------------------------------------------------------- 
/// \brief Default constructor for a Crack object. Zero-initializes all 
/// data members
///
// -----------------------------------------------------------------------
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

// ----------------------------------------------------------------------- 
/// \brief Sets the value of data members according to:
///
/// \param cx: abcissa of center of ellipse
/// \param cy: ordinate of center of ellipse
/// \param ax: length of the semi-major axis of the ellipse
/// \param ay: length of the semi-minor axis of the ellipse
/// \param sine: sine of the angle that major-axis forms with x-axis
/// \param cosine: cosine of the angle that major-axis forms with x-axis (redundant)
// -----------------------------------------------------------------------
void Crack::setter(float cx, float cy, float ax, float ay, float sine, float cosine) {
	c[0] = cx;
	c[1] = cy;
	a[0] = ax;
	a[1] = ay;
	trig[0] = sine;
	trig[1] = cosine;
}

// ----------------------------------------------------------------------- 
/// \brief Sets the value of data members according to the Crack source object
/// Effectively, overloads setter() to implement copy constructor 
///
// -----------------------------------------------------------------------
void Crack::setter(Crack & source){
	c[0] = source.c[0];
	c[1] = source.c[1];
	a[0] = source.a[0];
	a[1] = source.a[1];
	trig[0] = source.trig[0];
	trig[1] = source.trig[1];
}

// ----------------------------------------------------------------------- 
/// \brief Destructor
///
// -----------------------------------------------------------------------
Crack::~Crack() {
	clear();
}

// ----------------------------------------------------------------------- 
/// \brief Copy constructor
///
// -----------------------------------------------------------------------
Crack::Crack(Crack const & source) {
	copy(source);
}

// ----------------------------------------------------------------------- 
/// \brief Assignment(=) overloading. Calls copy constructor.
///
// -----------------------------------------------------------------------
Crack const & Crack::operator=(Crack const & other) {
	
	if (this != &other) {
		clear();
		copy(other);
	}
	return *this;	
}

// ----------------------------------------------------------------------- 
/// \brief Called by destructor. Clears all variables from memory.
///
// -----------------------------------------------------------------------
void Crack::clear() {
	delete[] a;
	delete[] c;
	delete[] trig;
	a = NULL;
	c = NULL;
	trig = NULL;
}

// ----------------------------------------------------------------------- 
/// \brief Implements copying from source. Called by copy constructor.
///
// -----------------------------------------------------------------------
void Crack::copy(Crack const & source) {

	dim = source.dim;
	for (int i = 0; i < dim; i++) {
		this->a[i] = source.a[i];
		this->c[i] = source.c[i];
		this->trig[i] = source.trig[i];
	}
}

// ----------------------------------------------------------------------- 
/// \brief Prints the data members of Crack object to STDOUT
///
// -----------------------------------------------------------------------
void Crack::print_info(){
	printf("\n");
	printf("Crack center at: (%0.2f, %0.2f)\n", c[0], c[1]);
	printf("Crack axes : (a: %0.2f, b: %0.2f)\n", a[0], a[1]);
	printf("Crack angle : %0.2f\n", asin(trig[0])*180.0/PI);
	printf("----------------------\n");
}


// ----------------------------------------------------------------------- 
//		Cracklist functions 
//
// -----------------------------------------------------------------------


// ----------------------------------------------------------------------- 
/// \brief Constructor to initialize a list of cracks with just one crack 
/// object
///
/// \param a -->  Crack object to populate the Cracklist object
// -----------------------------------------------------------------------
Cracklist::Cracklist(Crack& a){
	n_cracks = 1;
	listofCracks = new Crack[n_cracks];
	listofCracks[0].setter(a);
}


// ----------------------------------------------------------------------- 
/// \brief Constructor to randomly spawn n ellipses within the bound with 
/// random axes lengths and random orientations.
///
/// \param n -->  Number of random ellipses to be spwaned within object
// -----------------------------------------------------------------------
Cracklist::Cracklist(int n){
	this->dim = 2;
	n_cracks = n;
	listofCracks = new Crack[n];

	std::default_random_engine gen;

	std::uniform_real_distribution<float> a_dist(0.01,0.02);
	std::uniform_real_distribution<float> b_dist(0.05,0.06);
	std::uniform_real_distribution<float> angle(0.0,PI/2.0);
	std::uniform_real_distribution<float> cx_dist(0.2, 0.8);
	std::uniform_real_distribution<float> cy_dist(0.2, 0.8);
	//std::uniform_real_distribution<float> cz(-PI/2.0,PI/2.0);
	
	float bound = (MAXBOUND_X + MAXBOUND_Y)/2.0; // ~arbitrary value here. <br>
	// You think you can do better...ah, well you're right!

	float ang = angle(gen);
	float a = a_dist(gen)*bound;
	float b = b_dist(gen)*bound;
	float cx = cx_dist(gen)*bound;
	float cy = cy_dist(gen)*bound;
	float s = sin(ang);
	float c = cos(ang);

	for(int i=0; i < n_cracks; i++){
		listofCracks[i].setter(cx,cy,a,b,s,c);
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

// ----------------------------------------------------------------------- 
/// \brief Default destructor
// -----------------------------------------------------------------------
Cracklist::~Cracklist(){
	delete[] listofCracks;
	listofCracks = NULL;
}


// ----------------------------------------------------------------------- 
/// \brief Overload dereference operator([])
///
/// \param b --> Specifies which Crack object to be used from the given
/// Cracklist object.
// -----------------------------------------------------------------------
Crack & Cracklist::operator[](int b) {
	return listofCracks[b];	
}	



// int main(){
	
// 	Cracklist a(5,100);
// 	return 0;

// }