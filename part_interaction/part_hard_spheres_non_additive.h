#ifndef PARTICLE_INTERACTION_HARD_SPHERES_NON_ADD_H
#define PARTICLE_INTERACTION_HARD_SPHERES_NON_ADD_H

#include "part_interaction.h"
#include "particle.h"
#include "JRand_dyn.h"
#include <vector>
#include <cstdlib>
#include <cmath>

using namespace std;

class particle_interaction_hard_spheres_non_additive : public particle_interaction {

	private:
		JRand_dyn *rand_1;
		double polydispersity;
		double non_add_param;
		double int_range;
		double ave_diam;
		double largest_diam;
		double pi;

	public:
		particle_interaction_hard_spheres_non_additive(JRand_dyn *i_rand_1, double i_polydispersity, double i_non_add_param) {
			polydispersity = i_polydispersity;
			non_add_param = i_non_add_param;
			rand_1 = i_rand_1;
			ave_diam = 1.0;
			pi = 3.141592654;
			largest_diam = 1.0;
		}
		~particle_interaction_hard_spheres_non_additive() {;}

		double interaction_energy(double r_2, double type_1, double type_2);
		double interaction_width(double size_1, double size_2);
		int bonded(double r_2, double type_1, double type_2);
		void assign_type(particle *part_1);
		double interaction_range_max();
		double interaction_range();
		void c_interaction_range(double new_int_range);
		double box_muller_gaussian(double mu, double sigma);
};

#endif
