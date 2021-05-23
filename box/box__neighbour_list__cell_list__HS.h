#ifndef BOX_NEIGH_CELL_HS_H
#define BOX_NEIGH_CELL_HS_H

#include "box__neighbour_list__cell_list.h"

using namespace std;

class box__neighbour_list__cell_list__HS : public box__neighbour_list__cell_list {

	public:

		double smallest_scaled_separation;
		vector <int> smallest_separated_pair;

		box__neighbour_list__cell_list__HS(int i_number_particles, double i_temperature, double i_MC_dr, vector <double> i_box_size, vector <int> i_pb_flag, particle *i_particle, particle_interaction *i_part_interaction, JRand_dyn *i_rand_1, ext_field *i_ext_field_1, boundary_interaction *i_boundary_1) : box__neighbour_list__cell_list(i_number_particles, i_temperature, i_MC_dr, i_box_size, i_pb_flag, i_particle, i_part_interaction, i_rand_1, i_ext_field_1, i_boundary_1) {;}

		~box__neighbour_list__cell_list__HS(){;}

		double r_smallest_scaled_separation();
		void c_smallest_scaled_separation(double new_scaled_separation);
		vector <int> r_smallest_separated_pair();
		void c_smallest_separated_pair(int no_in_pair, int particle);
		void c_packing_fraction(double new_packing_fraction);
		double calc_scaled_separation(int i_part_1, int i_part_2);
		double calc_scaled_separation(int i_part_1, int i_part_2, double i_rij_2);
		void find_smallest_separation();
		double interaction_energy(int i_part_1, int i_part_2);
		double particle_energy(int i_part_1);
		void rescale_all(double scale_factor);
		void generate_initial();
		void generate_initial(FILE *in_file);

};

#endif
