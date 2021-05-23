#ifndef BOX_H
#define BOX_H

#include <vector>
#include "particle.h"
#include "part_interaction.h"

using namespace std;

class box {

	public:
		virtual double particle_energy(int part_1) = 0;
		virtual double total_energy_calc() = 0;
		virtual void displace_particle(int part_1, vector <double> displacement) = 0;

		virtual void dump_configuration(char out_file_name[], int seed, char traj_start[], double i_dr, double i_dv, char type_stamp[]) = 0;

		virtual int r_number_particles() = 0;
		virtual double r_temperature() = 0;
		virtual void c_temperature(double new_temperature) = 0;
		virtual double r_pressure() = 0;
		virtual void c_pressure(double new_pressure) = 0;
		virtual double r_volume() = 0;
		virtual double r_number_density() = 0;
		virtual double r_packing_fraction() = 0;
		virtual double r_ave_sigma_3() = 0;
		virtual void c_packing_fraction(double new_packing_fraction) = 0;
		virtual void c_box_size(vector <double> new_box_size) = 0;
		virtual void rescale_all(double scale_factor) = 0;
		virtual double r_timestamp() = 0;
		virtual void c_timestamp(double new_timestamp) = 0;
		virtual int r_dim() = 0;
		virtual double r_box_size(int i_dim) = 0;
		virtual int r_pb_flag(int i_dim) = 0;
		virtual vector <particle*> r_p_vector() = 0;
		virtual particle* r_particle(int part_1) = 0;
		virtual particle_interaction* r_particle_interaction() = 0;
		virtual void copy_coordinates(vector <particle*> coords_in) = 0;

		virtual double metric(int part_1, int part_2) = 0;
		virtual double metric_true(int part_1, int part_2) = 0;
		virtual int bond(int part_1, int part_2) = 0;
		virtual vector <int> r_particle_bond_list(int part_1) = 0;
		virtual vector <int> r_particle_closer_than_list(int part_1, double range) = 0;

		virtual void calc_centre_of_mass(int version_flag) = 0;
		virtual double r_centre_of_mass(int i_dim) = 0;
		virtual void reset_displacements() = 0;

};

#endif
