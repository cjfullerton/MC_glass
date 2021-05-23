#ifndef HARD_BOUND_SPH_H
#define HARD_BOUND_SPH_H

#include "boundary.h"
#include <vector>
#include <cmath>

//interaction between planar hard boundary and spherical particle

using namespace std;

class hard_boundary_sphere : public boundary_interaction {

	protected:
		vector <int> boundary_vector;
		vector <double> box_size;


	public:
		hard_boundary_sphere(vector <int> i_boundary_vector, vector <double> i_box_size) {

			boundary_vector = i_boundary_vector;
			box_size = i_box_size;

		}
		~hard_boundary_sphere() {;}

		double boundary_int_energy(vector <double> i_coords, double i_type);

};

#endif
