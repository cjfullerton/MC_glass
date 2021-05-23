#ifndef NO_BOUNDARY_H
#define NO_BOUNDARY_H

#include "boundary.h"
#include <vector>

//interaction between planar hard boundary and spherical particle

using namespace std;

class no_boundary : public boundary_interaction {

	public:
		no_boundary() {;}
		~no_boundary() {;}

		double boundary_int_energy(vector <double> i_coords, double i_type);

};

#endif
