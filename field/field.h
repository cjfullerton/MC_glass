#ifndef FIELD_H
#define FIELD_H

#include <vector>

using namespace std;

//abstract base class for energy due to uniform field applied to whole system

class ext_field {

	public:
		virtual double field_energy(vector <double> i_coords, double type_in) = 0;

};

#endif
