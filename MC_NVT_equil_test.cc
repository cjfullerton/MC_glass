#include "MC_NVT_equil_test.h"

int main(int argc, char *argv[]){

	vector <int> pb_flag;
	vector <double> box_size;
	vector <double> time_series;
	vector <int> step;
	vector <vector <double> > time_out;
	int seed, next_step, jobid;
	double L, phi_f, phi_i, time_step, log_step, time_end, count, waiting_time, simulation_time, time_current, lambda, u, dv, dr;
	FILE *config_in;
	FILE *sys_info_out;
	char file_name[1024] = {0};
	char sys_info_file_name[1024] = {0};
	char restart_file_name[1024] = {0};
	vector <std::string> file_names;

	stringstream convert_phi_f(argv[1]); convert_phi_f >> phi_f;
	stringstream convert_u(argv[2]); convert_u >> u;
	stringstream convert_lambda(argv[3]); convert_lambda >> lambda;
	stringstream convert_dr(argv[4]); convert_dr >> dr;
	stringstream convert_jobid(argv[6]); convert_jobid >> jobid;
	dv = 0.0;

	JRand_dyn *rand_1 = new JRand_dyn();
	seed = (unsigned int)(jobid+time(NULL))%1000;
	rand_1->srandom(seed);

	pb_flag.resize(DIM); box_size.resize(DIM);
	pb_flag[0] = PB_X; pb_flag[1] = PB_Y; pb_flag[2] = PB_Z;
	box_size[0] = BOX_SIZE_X; box_size[1] = BOX_SIZE_Y; box_size[2] = BOX_SIZE_Z;

	part_poly_sphere *particle_1 = new part_poly_sphere(DIM);
	particle_interaction_sticky_spheres *int_poly_sticky_1 = new particle_interaction_sticky_spheres(lambda, u, 0.0, rand_1);
	no_field *field_1 = new no_field();
	no_boundary *boundary_1 = new no_boundary();

	config_in = fopen(argv[5], "r");

	box__neighbour_list__cell_list *box_ti = new box__neighbour_list__cell_list(NUMBER_PARTICLES, TEMPERATURE, dr, box_size, pb_flag, particle_1, int_poly_sticky_1, rand_1, field_1, boundary_1);
	box_ti -> generate_initial(config_in);

	phi_i = box_ti->r_packing_fraction();

	fclose(config_in);

	integrator_MC_NPT *integrator_1 = new integrator_MC_NPT(box_ti, dr, 0.0, rand_1);

	time_current = 0.0;

	box_ti->c_timestamp(0.0);
	box_ti -> c_packing_fraction(phi_f);


	printf("%f\n", box_ti->total_energy_calc());


	return 0;
}
