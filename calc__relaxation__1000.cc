#include "calc__relaxation__1000.h"

int main(int argc, char *argv[]){

	vector <int> pb_flag;
	vector <double> box_size;
	int seed, slice, final_slice;
	double L, phi, pressure, time_step, log_step, time_end;
	FILE *config_in;
	char file_name[1024];

	JRand_dyn *rand_1 = new JRand_dyn();
	seed = (unsigned int)time(NULL)%1000;
	rand_1->srandom(seed);

	stringstream convert_final_slice(argv[2]); convert_final_slice >> final_slice;

	pb_flag.resize(DIM); box_size.resize(DIM);
	pb_flag[0] = PB_X; pb_flag[1] = PB_Y; pb_flag[2] = PB_Z;
	box_size[0] = BOX_SIZE_X; box_size[1] = BOX_SIZE_Y; box_size[2] = BOX_SIZE_Z;

	part_poly_sphere *particle_1 = new part_poly_sphere(DIM);
	particle_interaction_hard_spheres *int_poly_hard_1 = new particle_interaction_hard_spheres(rand_1, 0.0);
	no_field *field_1 = new no_field();
	no_boundary *boundary_1 = new no_boundary();

	config_in = fopen(argv[1], "r");

	box__neighbour_list__cell_list__HS *box_t0 = new box__neighbour_list__cell_list__HS(NUMBER_PARTICLES, TEMPERATURE, DR, box_size, pb_flag, particle_1, int_poly_hard_1, rand_1, field_1, boundary_1);
	box_t0 -> generate_initial(config_in);

	fclose(config_in);

	config_in = fopen(argv[1], "r");

	box__neighbour_list__cell_list__HS *box_ti = new box__neighbour_list__cell_list__HS(NUMBER_PARTICLES, TEMPERATURE, DR, box_size, pb_flag, particle_1, int_poly_hard_1, rand_1, field_1, boundary_1);
	box_ti -> generate_initial(config_in);

	sys_info__mean_squared_displacement *calc_MSD = new sys_info__mean_squared_displacement();
	calc_MSD->box_t0_in(box_t0);
	calc_MSD->box_ti_in(box_ti);

	sys_info__scattering_functions *calc_sf = new sys_info__scattering_functions();
	calc_sf->box_t0_in(box_t0);
	calc_sf->box_ti_in(box_ti);

	slice  = 0;

	while(slice < final_slice - 1) {

		calc_MSD->calc_msd(1,1);
		calc_sf->calc_fkt(7.80,1,1);
		calc_sf->calc_fkt(7.80,0,1);

		box_ti -> read_configuration(config_in);

		slice++;

	}

	calc_MSD->calc_msd(1,1);
	calc_sf->calc_fkt(7.80,1,1);
	calc_sf->calc_fkt(7.80,0,1);

	return 0;
}
