#include "grav_field.h"

double grav_field::field_energy(vector <double> i_coords, double type_in) {

	double field_e = 0.0;

	type_in = 1.0;

	for(int j = 0; j<i_coords.size(); j++) {

		if(field_direction[j] == 1) {

			field_e = field_e + type_in*type_in*type_in*mg*i_coords[j];

		}
	}

	return field_e;

}
