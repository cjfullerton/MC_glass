#ifndef BOUND_H
#define BOUND_H

#include <vector>

using namespace std;

//abstract base class for interaction energy between particles and boundary of box

class boundary_interaction {

	public:
		virtual double boundary_int_energy(vector <double> i_coords, double i_type) = 0;	

};
#endif
