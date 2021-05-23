#include "sys_info__scattering_functions.h"

void sys_info__scattering_functions::box_t0_in(box *i_box_t0) { box_t0 = i_box_t0;}

void sys_info__scattering_functions::box_ti_in(box *i_box_ti) { box_ti = i_box_ti;}

void sys_info__scattering_functions::calc_fkt(double q, int com_flag, int version_flag) {

	double dr, df, r2_acc = 0.0, q_dot_r, sin_q_dot_r, f_qt_acc = 0.0, time;

	if(com_flag) {
		box_ti->calc_centre_of_mass(version_flag);
		box_t0->calc_centre_of_mass(version_flag);
	}

	for(int i = 0; i < box_ti->r_number_particles(); i++) {

		r2_acc = 0.0;

		for(int j = 0; j < box_ti->r_dim(); j++) {
			
			if(version_flag == 1) dr = box_ti->r_particle(i)->r_coordinate_true(j) - box_t0->r_particle(i)->r_coordinate_true(j);
			if(version_flag == 2) dr = box_ti->r_particle(i)->r_coordinate_box(j) - box_t0->r_particle(i)->r_coordinate_box(j);
			if(version_flag == 3) dr = box_ti->r_particle(i)->r_coordinate_affine(j) - box_t0->r_particle(i)->r_coordinate_affine(j);

			if(com_flag) dr = dr - (box_ti->r_centre_of_mass(j) - box_t0->r_centre_of_mass(j));

			r2_acc += dr*dr;

		}

		q_dot_r = q*sqrt(r2_acc);
		sin_q_dot_r = sin(q_dot_r);
		if(q_dot_r == 0.0) {q_dot_r = 1.0; sin_q_dot_r = 1.0;}
		f_qt_acc += sin_q_dot_r/q_dot_r;
	}

	time = box_ti->r_timestamp() - box_t0->r_timestamp();

	printf("#FQT_%d %f %f %f\n", version_flag, time, f_qt_acc/box_ti->r_number_particles(), q);

}

void sys_info__scattering_functions::calc_fkt_II(double q, int com_flag, int version_flag) {

	double dr, df, r2_acc = 0.0, q_dot_r, sin_q_dot_r, f_qt_acc = 0.0, time;

	if(com_flag) {
		box_ti->calc_centre_of_mass(version_flag);
		box_t0->calc_centre_of_mass(version_flag);
	}

	for(int i = 0; i < box_ti->r_number_particles(); i++) {

		r2_acc = 0.0;

		for(int j = 0; j < box_ti->r_dim(); j++) {
			
			if(version_flag == 1) dr = box_ti->r_particle(i)->r_coordinate_true(j) - box_t0->r_particle(i)->r_coordinate_true(j);
			if(version_flag == 2) dr = box_ti->r_particle(i)->r_coordinate_box(j) - box_t0->r_particle(i)->r_coordinate_box(j);
			if(version_flag == 3) dr = box_ti->r_particle(i)->r_coordinate_affine(j) - box_t0->r_particle(i)->r_coordinate_affine(j);

			if(com_flag) dr = dr - (box_ti->r_centre_of_mass(j) - box_t0->r_centre_of_mass(j));

			r2_acc += dr*dr;

		}

		q_dot_r = q*sqrt(r2_acc);
		sin_q_dot_r = sin(q_dot_r);
		if(q_dot_r == 0.0) {q_dot_r = 1.0; sin_q_dot_r = 1.0;}
		f_qt_acc += sin_q_dot_r/q_dot_r;
	}

	time = box_ti->r_timestamp() - box_t0->r_timestamp();

	printf("%f %f %f %f\n", box_t0->r_timestamp(), box_ti->r_timestamp(), f_qt_acc/box_ti->r_number_particles(), q);

}

void sys_info__scattering_functions::calc_overlap(double a, int com_flag) {

	double overlap_acc = 0.0, dr, r2_acc, r, time;

	if(com_flag) {
		box_ti->calc_centre_of_mass(1);
		box_t0->calc_centre_of_mass(1);
	}



	for(int i = 0; i < box_ti->r_number_particles(); i++) {

		r2_acc = 0.0;

		for(int j = 0; j < box_ti->r_dim(); j++) {
			
			dr = box_ti->r_particle(i)->r_coordinate_true(j) - box_t0->r_particle(i)->r_coordinate_true(j);

			if(com_flag) dr = dr - (box_ti->r_centre_of_mass(j) - box_t0->r_centre_of_mass(j));

			r2_acc += dr*dr;
		}

		r = sqrt(r2_acc);

		if(r <= a) overlap_acc += 1.0;

	}

	time = box_ti->r_timestamp() - box_t0->r_timestamp();

	printf("%f %f %f\n", time, overlap_acc/box_ti->r_number_particles(), a);

}

void sys_info__scattering_functions::calc_overlap_per_particle(double a, int com_flag) {

	double dr, r2_acc, r;

	overlap_per_particle.clear();
	overlap_per_particle.resize(box_ti->r_number_particles(), 0.0);

	if(com_flag) {
		box_ti->calc_centre_of_mass(1);
		box_t0->calc_centre_of_mass(1);
	}

	for(int i = 0; i < box_ti->r_number_particles(); i++) {

		r2_acc = 0.0;

		for(int j = 0; j < box_ti->r_dim(); j++) {
			
			dr = box_ti->r_particle(i)->r_coordinate_true(j) - box_t0->r_particle(i)->r_coordinate_true(j);

			if(com_flag) dr = dr - (box_ti->r_centre_of_mass(j) - box_t0->r_centre_of_mass(j));

			r2_acc += dr*dr;
		}

		r = sqrt(r2_acc);

		if(r <= a) overlap_per_particle[i] += 1.0;

	}

}

void sys_info__scattering_functions::coarse_grain_overlap_per_particle(double a, int com_flag, double rcg) {

	double r2, rcg2, rcg2sigma2, sigma_i;

	rcg2 = rcg*rcg;

	calc_overlap_per_particle(a, com_flag);

	overlap_coarse_grained.resize(box_ti->r_number_particles());
	
	for(int i = 0; i < box_ti->r_number_particles(); i++) {
		overlap_coarse_grained[i].clear();
		overlap_coarse_grained[i].resize(2, 0.0);
	}

	for(int i = 0; i < box_ti->r_number_particles(); i++) { //NB doing full double loop so that coarse graining radius for particle i depends on sigma_i
		
		sigma_i = box_ti->r_particle(i)->r_size();
		rcg2sigma2 = rcg*sigma_i*sigma_i;

		for(int j = 0; j < box_ti->r_number_particles(); j++) {

			r2 = box_ti->metric(i, j);
				
			if(r2 < rcg2sigma2) {
				overlap_coarse_grained[i][0] += overlap_per_particle[j];
				overlap_coarse_grained[i][1] += 1.0;

			}
		}
	}

}

void sys_info__scattering_functions::dump_number_particles_with_overlap_less_than(double overlap_cut, double a_in, int com_flag_in, double rcg_in) {

	double particle_count, time;

	particle_count = 0.0;

	coarse_grain_overlap_per_particle(a_in, com_flag_in, rcg_in);

	for(int i = 0; i < box_ti->r_number_particles(); i++) {

		if(overlap_coarse_grained[i][0]/overlap_coarse_grained[i][1] < overlap_cut) particle_count += 1.0;

	}

	 time = box_ti->r_timestamp() - box_t0->r_timestamp();
	 
	 printf("%f %f %d %f %f\n", time, particle_count, box_ti->r_number_particles(), a_in, rcg_in);

}

double sys_info__scattering_functions::return_number_particles_with_overlap_less_than(double overlap_cut, double a_in, int com_flag_in, double rcg_in) {

	double particle_count;

	particle_count = 0.0;

	coarse_grain_overlap_per_particle(a_in, com_flag_in, rcg_in);

	for(int i = 0; i < box_ti->r_number_particles(); i++) {

		if(overlap_coarse_grained[i][0]/overlap_coarse_grained[i][1] < overlap_cut) particle_count += 1.0;

	}

	return particle_count;

}

void sys_info__scattering_functions::dump_config_with_coarse_grained_overlap(double a, int com_flag, double rcg) {

	coarse_grain_overlap_per_particle(a, com_flag, rcg);

	printf("%d\n", box_ti->r_number_particles());

	for(int j = 0; j < box_ti -> r_dim(); j++) printf("%f ", box_ti->r_box_size(j));

	printf("config_with_coarse_grained_per_particle_overlaps_t0_coords\n");

	for(int i = 0; i < box_ti->r_number_particles(); i++) {

		printf("%f ",box_t0->r_particle(i)->r_size());

		for(int j = 0; j < box_ti->r_dim(); j++) printf("%f ", box_t0->r_particle(i)->r_coordinate_box(j));

		printf("%f\n", overlap_coarse_grained[i][0]/overlap_coarse_grained[i][1]);

	}

}

void sys_info__scattering_functions::dump_config_with_coarse_grained_overlap_II(double a, int com_flag, double rcg) {

	double coord_shift_com, coord_shift_com_pb;

	coarse_grain_overlap_per_particle(a, com_flag, rcg);

	printf("%d\n", box_ti->r_number_particles());

	for(int j = 0; j < box_ti -> r_dim(); j++) printf("%f ", box_ti->r_box_size(j));

	printf("config_with_coarse_grained_per_particle_overlaps_ti_coords\n");

	box_ti->calc_centre_of_mass(1);
	box_t0->calc_centre_of_mass(1);

	for(int i = 0; i < box_ti->r_number_particles(); i++) {

		printf("%f ",box_ti->r_particle(i)->r_size());

		for(int j = 0; j < box_ti->r_dim(); j++) {

			coord_shift_com = box_ti->r_particle(i)->r_coordinate_box(j) - (box_ti->r_centre_of_mass(j) - box_t0->r_centre_of_mass(j));

			if(coord_shift_com < 0.0) coord_shift_com_pb = coord_shift_com + box_ti->r_box_size(j);
			else if(coord_shift_com >= box_ti->r_box_size(j)) coord_shift_com_pb = coord_shift_com - box_ti->r_box_size(j);
			else coord_shift_com_pb = coord_shift_com;

			printf("%f ", coord_shift_com_pb);
		}

		printf("%f\n", overlap_coarse_grained[i][0]/overlap_coarse_grained[i][1]);

	}

}
