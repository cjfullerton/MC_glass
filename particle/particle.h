#ifndef PARTICLE_H
#define PARTICLE_H

#include <cstdlib>
#include <vector>

using namespace std;

class particle {

	public:
		particle() {;}
		~particle() {;}

		virtual particle* clone() = 0;

		virtual double r_coordinate_box(int j) = 0;
		virtual vector <double> r_coordinate_box_vec() = 0;
		virtual double r_coordinate_true(int j) = 0;
		virtual vector <double> r_coordinate_true_vec() = 0;
		virtual double r_coordinate_affine(int j) = 0;
		virtual vector <double> r_coordinate_affine_vec() = 0;
		virtual int r_type() = 0;
		virtual double r_size() = 0;
		
		virtual void c_coordinate_box(int j, double new_coord_box) = 0;
		virtual void c_coordinate_true(int j, double new_coord_true) = 0;
		virtual void c_coordinate_affine(int j, double new_coord_affine) = 0;
		virtual void c_type(int new_type) = 0;
		virtual void c_size(double new_size) = 0;

};

#endif
