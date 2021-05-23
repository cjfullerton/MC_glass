#include "hard_boundary.h"

double hard_boundary_sphere::boundary_int_energy(vector <double> i_coords, double i_type) {

	double bound_e = 0.0;

	for(int j = 0; j<i_coords.size(); j++) {

		if(boundary_vector[j] == 1) { 

			if(i_coords[j] < i_type/2.0) bound_e = bound_e + NAN;
			else if(box_size[j] - i_coords[j] < i_type/2.0) bound_e = bound_e + NAN;
		}
	}

	return bound_e;

}
