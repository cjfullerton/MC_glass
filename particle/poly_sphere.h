#ifndef PART_POLY_SPHERE_H
#define PART_POLY_SPHERE_H

#include "particle.h"
#include <vector>
#include <cmath>

using namespace std;

class part_poly_sphere : public particle {

	protected:
		vector <double> coordinate_box;
		vector <double> coordinate_true;
		vector <double> coordinate_affine;
		int type;
		double size;
		int dim;

	public:

		part_poly_sphere(int i_dim) {

			dim = i_dim;
			type = 0;
			size = 1.0;

			coordinate_box.resize(dim, 0.0);
			coordinate_true.resize(dim, 0.0);
			coordinate_affine.resize(dim, 0.0);

		}
		part_poly_sphere(part_poly_sphere *i_part_1) {

			dim = i_part_1->dim;
			type = i_part_1->type;
			size = i_part_1->size;

			coordinate_box = i_part_1->coordinate_box;
			coordinate_true = i_part_1->coordinate_true;
			coordinate_affine = i_part_1->coordinate_affine;
		}
		~part_poly_sphere() {;}

		part_poly_sphere* clone() { return new part_poly_sphere(*this); }

		double r_coordinate_box(int j);
		vector <double> r_coordinate_box_vec();
		double r_coordinate_true(int j);
		vector <double> r_coordinate_true_vec();
		double r_coordinate_affine(int j);
		vector <double> r_coordinate_affine_vec();
		int r_type();
		double r_size();

		void c_coordinate_box(int j, double new_coordinate_box);
		void c_coordinate_true(int j, double new_coordinate_affine);
		void c_coordinate_affine(int j, double new_coordinate_true);
		void c_type(int new_type);
		void c_size(double new_size);


};

#endif
