#ifndef CRACK_H
#define CRACK_H

#include <stddef.h>
#include "params.h"

using namespace std;

/**
@file crack.h
\brief Defines a crack in terms of 2D ellipses.
*/
class Crack {

protected:
	
	void clear();
	void copy(Crack const & source);

public:
	float * c; ///< Holds the coordinates of the center of the ellipse
	float * a; ///< Holds the lengths of the major and minor axes of the ellipse
	float * trig; ///< Holds the sine and cosine of the angle the major axis of the ellipse makes with the problem's x axis
	int dim; ///< Obsolete. Is generally equal to DIM from params.h
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
	int n_cracks; ///< Number of cracks in the Cracklist
	Crack* listofCracks;
	Crack & operator[](int);
	Cracklist(int);
	Cracklist();
	Cracklist(Crack&);
	~Cracklist();
} ;

#endif