#ifndef PART_INTERACTION_H
#define PART_INTERACTION_H

#include "particle.h"
#include <cstdlib>

class particle_interaction {

	public:
		virtual double interaction_energy(double r_2, double size_1, double size_2) = 0;
		virtual double interaction_width(double size_1, double size_2) = 0;
		virtual double interaction_range_max() = 0;
		virtual double interaction_range(double size) = 0;
		virtual void c_interaction_range(double new_interaction_range) = 0;
		virtual void assign_type(particle *part_1) = 0;
		virtual int bonded(double r_2, double size_1, double size_2) = 0;
};

#endif
