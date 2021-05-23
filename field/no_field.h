#ifndef NO_FIELD_H
#define NO_FIELD_H

#include "field.h"

using namespace std;

class no_field : public ext_field {

	public:
		no_field() {;}
		~no_field() {;}

		double field_energy(vector <double> i_coords, double type_in);

};


#endif

