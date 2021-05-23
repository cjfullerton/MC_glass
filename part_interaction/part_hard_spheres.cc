#include "part_hard_spheres.h"

double particle_interaction_hard_spheres::interaction_energy(double r_2, double type_1, double type_2) {

	double hard_core, energy;

	hard_core = (type_1 + type_2)/2.0;
	
	if(r_2 < hard_core*hard_core) energy = NAN;
	else energy = 0.0;

	return energy;
}

double particle_interaction_hard_spheres::interaction_width(double type_1, double type_2) {

	return (type_1 + type_2)/2.0;

}

int particle_interaction_hard_spheres::bonded(double r_2, double type_1, double type_2) {

	double hard_core;
	int bond;

	hard_core = (type_1 + type_2)/2.0;

	if(r_2 == hard_core*hard_core) bond = 1;
	else bond = 0;

	return bond;
}

void particle_interaction_hard_spheres::assign_type(particle *part_1) {

	double x;

	if(polydispersity > 0.0) x = box_muller_gaussian(ave_diam, polydispersity);
	else x = ave_diam;

	if(x > largest_diam) largest_diam = x;

	part_1->c_type(x);

}

void particle_interaction_hard_spheres::c_interaction_range(double new_interaction_range) { largest_diam = new_interaction_range; }

double particle_interaction_hard_spheres::interaction_range(double size) { return size; }

double particle_interaction_hard_spheres::interaction_range_max() { return largest_diam; }

double particle_interaction_hard_spheres::box_muller_gaussian(double mu, double sigma) {

	double r1, r2, rand_gauss;

	r1 = rand_1->prob_gen();
	r2 = rand_1->prob_gen();

	rand_gauss = mu + sigma*sqrt(-2.0*log(r1))*cos(2.0*pi*r2);

	return rand_gauss;

}
