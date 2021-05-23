#ifndef GRAV_FIELD_H
#define GRAV_FIELD_H

#include "field.h"
#include <vector>
#include <cmath>

using namespace std;

class grav_field : public ext_field {

	protected:
		double mg;
		vector <int> field_direction;

	public:
		grav_field(vector <int> i_field_direction, double i_mg) {

			field_direction = i_field_direction; mg = i_mg;

		}
		~grav_field() {;}

		double field_energy(vector <double> i_coords, double type_in);

};


#endif
