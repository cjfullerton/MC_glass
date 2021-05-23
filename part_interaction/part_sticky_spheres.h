#ifndef PARTICLE_INTERACTION_STICKY_SPHERES_H
#define PARTICLE_INTERACTION_STICKY_SPHERES_H

#include "part_interaction.h"
#include "particle.h"
#include "JRand_dyn.h"
#include <vector>
#include <cstdlib>
#include <cmath>

using namespace std;

class particle_interaction_sticky_spheres : public particle_interaction {

	private:
		JRand_dyn *rand_1;
		double lambda; //width of square well - should this be here, or should it go in particle?
		double well_depth; //depth of square well - should this be here, or go in particle class?
		double polydispersivity;
		double ave_diam;
		double largest_diam;
		double pi;

	public:
		particle_interaction_sticky_spheres(double i_lambda, double i_well_depth, double i_polydispersivity, JRand_dyn *i_rand_1) {
			rand_1 = i_rand_1;
			lambda = i_lambda;
			well_depth = i_well_depth; polydispersivity = i_polydispersivity; ave_diam = 1.0;
			largest_diam = 0.0;
			pi = 3.141592654;
		}
		~particle_interaction_sticky_spheres() {;}

		double interaction_energy(double r_2, double type_1, double type_2);
		double interaction_width(double type_1, double type_2);
		int bonded(double r_2, double type_1, double type_2);
		void assign_type(particle *part_1);
		double interaction_range(double size);
		double interaction_range_max();
		void c_interaction_range(double new_int_range);
		double box_muller_gaussian(double mu, double sigma);
};

#endif
