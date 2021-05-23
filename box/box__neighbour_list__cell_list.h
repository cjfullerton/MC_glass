#ifndef BOX_NEIGH_CELL_H
#define BOX_NEIGH_CELL_H

#include "box.h"
#include "particle.h"
#include "part_interaction.h"
#include "field.h"
#include "boundary.h"
#include "JRand_dyn.h"

#include <cmath>
#include <cstdio>

#define SKIN 0.4

#define CELLS_FLAG 0
#define NLIST_FLAG 1

#define PI 3.141592654

using namespace std;

class box__neighbour_list__cell_list : public box {

	protected:
		vector <particle*> particle_list;
		vector <double> box_size;
		vector <int> pb_flag;
		JRand_dyn *rand_1;
		particle_interaction *part_interaction_1;
		ext_field *ext_field_1;
		boundary_interaction *boundary_1;

		int dim;
		double temperature;
		double pressure;
		double volume;
		double packing_fraction;
		double ave_sigma3;
		int number_particles;
		int cells_flag;
		int nlist_flag;
		int cells_done;
		int nlist_done;
		double MC_dr; //Used only for cell size in the case that cells and no neighbour list are used!
		double nlist_rebuild_limit;
		double nlist_cut_2;
		double cell_size_aim;
		int total_cells;
		int max_nlist_occup;
		int max_cells_occup;
		double timestamp;
		vector<int> number_cells;
		vector<double> cell_size;
		vector< vector<double> > cell_neighbours;
		vector< vector<int> > cell_contents;
		vector< vector<int> > part_cell;
		vector< vector<int> > nlist;
		vector< vector<double> > disp_since_rebuild;
		vector <double> centre_of_mass;

	public:
		box__neighbour_list__cell_list(int i_number_particles, double i_temperature, double i_MC_dr, vector <double> i_box_size, vector <int> i_pb_flag, particle *i_particle, particle_interaction *i_part_interaction, JRand_dyn *i_rand_1, ext_field *i_ext_field_1, boundary_interaction *i_boundary_1) {

			number_particles = i_number_particles;
			temperature = i_temperature;
			dim = i_box_size.size();
			MC_dr = i_MC_dr;

			cells_flag = CELLS_FLAG;
			nlist_flag = NLIST_FLAG;
			cells_done = 0;
			nlist_done = 0;
			
			cell_size_aim = 1.1 + MC_dr;

			box_size = i_box_size;
			pb_flag = i_pb_flag;
			part_interaction_1 = i_part_interaction;
			ext_field_1 = i_ext_field_1;
			boundary_1 = i_boundary_1;
			rand_1 = i_rand_1;

			for(int i = 0; i < number_particles; i++) particle_list.push_back(i_particle->clone());

			centre_of_mass.resize(dim, 0.0);

			if(cells_flag) {
				setup_cell_arrays();
				setup_cell_neighbours();
			}
		}
		~box__neighbour_list__cell_list() {;}

		int r_dim();
		double r_box_size(int i_dim);
		int r_pb_flag(int i_dim);
		double r_temperature();
		void c_temperature(double new_temperature);
		double r_pressure();
		void c_pressure(double new_pressure);
		double r_volume();
		double r_number_density();
		double r_packing_fraction();
		double r_ave_sigma_3();
		void c_packing_fraction(double new_packing_fraction);
		void c_box_size(vector <double> new_box_size);
		double r_timestamp();
		void c_timestamp(double new_timestamp);
		int r_number_particles();
		particle* r_particle(int part_i);
		particle_interaction* r_particle_interaction();
		vector <particle*> r_p_vector();
		void copy_coordinates(vector <particle*> coords_in);

		//cell list setup and management
		void setup_cell_arrays();
		void setup_cell_neighbours();
		int cell_number(int part);
		void assign_to_cell(int part);
		void assign_all_to_cells();
		void remove_from_cell(int part);
		void check_and_update_cell(int part);

		//neighbour list setup and management
		void setup_nlist_arrays();
		void setup_nlist();
		void check_particle_nlist_and_update(int part_1);

		//initial state generation, metric calculation of various flavours and MC move function
		void generate_initial();
		void generate_initial(FILE *in_file);
		double metric(int part_1, int part_2);
		double metric_true(int part_1, int part_2);
		void displace_particle(int part_1, vector <double> displacement);
		void rescale_all(double scale_factor);

		//energy, image wipe, bonds
		double total_energy_calc();
		double particle_energy(int i);
		double interaction_energy(int part_1, int part_2);
		int bond(int part_1, int part_2);
		vector <int> r_particle_bond_list(int part_1);
		vector <int> r_particle_closer_than_list(int part_1, double range);
		void calc_centre_of_mass(int version_flag);
		double r_centre_of_mass(int i_dim);
		void reset_displacements();

		//read and write configurations
		void dump_configuration(char out_file_name[], int seed, char traj_start[], double i_dr, double i_dv, char type_stamp[]);
		void dump_restart(char out_file_name[], int jobid, char traj_start[], double i_dr, double i_dv, char type_stamp[]);
		void ignore_configuration(FILE *in_file);
		void read_configuration(FILE *in_file);

};

#endif
