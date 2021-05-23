#include "part_hard_spheres_non_additive.h"

double particle_interaction_hard_spheres_non_additive::interaction_energy(double r_2, double type_1, double type_2) {

	double hard_core, energy, non_add_correction, width;

	hard_core = (type_1 + type_2)/2.0;
	non_add_correction = (1.0 + non_add_param*abs(type_1-type_2));
	width = non_add_correction*hard_core;
	
	if(r_2 < width*width) energy = NAN;
	else energy = 0.0;

	return energy;
}

double particle_interaction_hard_spheres_non_additive::interaction_width(double type_1, double type_2) {

	double hard_core, non_add_correction;

	hard_core = (type_1 + type_2)/2.0;
	non_add_correction = (1.0 + non_add_param*abs(type_1-type_2));

	return hard_core*non_add_correction;

}

int particle_interaction_hard_spheres_non_additive::bonded(double r_2, double type_1, double type_2) {

	double hard_core, non_add_correction, width;
	int bond;

	hard_core = (type_1 + type_2)/2.0;
	non_add_correction = (1.0 + non_add_param*abs(type_1-type_2));
	width = non_add_correction*hard_core;

	if(r_2 == width*width) bond = 1;
	else bond = 0;

	return bond;
}

void particle_interaction_hard_spheres_non_additive::assign_type(particle *part_1) {

	double x;

	if(polydispersity > 0.0) x = box_muller_gaussian(ave_diam, polydispersity);
	else x = ave_diam;

	if(x > largest_diam) largest_diam = x;

	part_1->c_type(x);

}

void particle_interaction_hard_spheres_non_additive::c_interaction_range(double new_interaction_range) { largest_diam = new_interaction_range; }

double particle_interaction_hard_spheres_non_additive::interaction_range(double size) { return size; }

double particle_interaction_hard_spheres_non_additive::interaction_range_max() { return largest_diam; }

double particle_interaction_hard_spheres_non_additive::box_muller_gaussian(double mu, double sigma) {

	double r1, r2, rand_gauss;

	r1 = rand_1->prob_gen();
	r2 = rand_1->prob_gen();

	rand_gauss = mu + sigma*sqrt(-2.0*log(r1))*cos(2.0*pi*r2);

	return rand_gauss;

}
