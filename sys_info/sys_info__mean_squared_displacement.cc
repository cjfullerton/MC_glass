#include "sys_info__mean_squared_displacement.h"

void sys_info__mean_squared_displacement::box_t0_in(box *i_box_t0) { box_t0 = i_box_t0; }

void sys_info__mean_squared_displacement::box_ti_in(box *i_box_ti) { box_ti = i_box_ti; }

void sys_info__mean_squared_displacement::calc_msd(int com_flag, int version_flag) {

	double dr, time, dr2_acc = 0.0,  msd_acc = 0.0, r2_acc = 0.0, df;

	if(com_flag) {
		box_ti->calc_centre_of_mass(version_flag);
		box_t0->calc_centre_of_mass(version_flag);
	}

	for(int i = 0; i < box_ti->r_number_particles(); i++) {

		dr2_acc = 0.0;

		for(int j = 0; j < box_ti->r_dim(); j++) {

			if(version_flag == 1) dr = box_ti->r_particle(i)->r_coordinate_true(j) - box_t0->r_particle(i)->r_coordinate_true(j);
			if(version_flag == 2) dr = box_ti->r_particle(i)->r_coordinate_box(j) - box_t0->r_particle(i)->r_coordinate_box(j);
			if(version_flag == 3) dr = box_ti->r_particle(i)->r_coordinate_affine(j) - box_t0->r_particle(i)->r_coordinate_affine(j);

			if(com_flag) dr = dr - (box_ti->r_centre_of_mass(j) - box_t0->r_centre_of_mass(j));

			dr2_acc += dr*dr;

		}

		msd_acc += dr2_acc;

	}

	time = box_ti->r_timestamp() - box_t0->r_timestamp();

	printf("#MSD_%d %f %f\n", version_flag, time, msd_acc/box_ti->r_number_particles());

}

void sys_info__mean_squared_displacement::calc_van_hove_function(int number_bins, int com_flag, int version_flag) {

	double dr, dr2_acc;
	int bin_no;

	if(com_flag) {
		box_ti->calc_centre_of_mass(version_flag);
		box_t0->calc_centre_of_mass(version_flag);
	}


	dr2_max = 0.0; dr2_min = 1000000000.0;
	MSD_list.clear();

	for(int i = 0; i < box_ti->r_number_particles(); i++) {

		dr2_acc = 0.0;

		for(int j = 0; j < box_ti->r_dim(); j++) {

			if(version_flag == 1) dr = box_ti->r_particle(i)->r_coordinate_true(j) - box_t0->r_particle(i)->r_coordinate_true(j);
			if(version_flag == 2) dr = box_ti->r_particle(i)->r_coordinate_box(j) - box_t0->r_particle(i)->r_coordinate_box(j);
			if(version_flag == 3) dr = box_ti->r_particle(i)->r_coordinate_affine(j) - box_t0->r_particle(i)->r_coordinate_affine(j);

			if(com_flag) dr = dr - (box_ti->r_centre_of_mass(j) - box_t0->r_centre_of_mass(j));

			dr2_acc += dr*dr;

		}

		if(dr2_acc > dr2_max) dr2_max = dr2_acc;
		if(dr2_acc < dr2_min) dr2_min = dr2_acc;

		MSD_list.push_back(dr2_acc);

	}

	bin_size = (dr2_max - dr2_min)/number_bins;

	MSD_count.clear();
	MSD_count.resize(number_bins + 1, 0.0);

	for(int i = 0; i < MSD_list.size(); i++) {

		bin_no = floor((MSD_list[i]-dr2_min)/bin_size);
		MSD_count[bin_no] += 1.0;

	}
	

}

void sys_info__mean_squared_displacement::print_van_hove_function() {

	double time;

	time = box_ti->r_timestamp() - box_t0->r_timestamp();

	for(int i = 0; i < MSD_count.size(); i++) if(MSD_count[i] > 0) printf("%f %f %f %f %d\n", time, dr2_min + i*bin_size, MSD_count[i], bin_size, MSD_list.size());

	printf("\n\n");

}

void sys_info__mean_squared_displacement::calc_displacement_distribution(int number_bins, int com_flag) {

	double dr;
	int bin_no;

	if(com_flag) {
		box_ti->calc_centre_of_mass(1);
		box_t0->calc_centre_of_mass(1);
	}


	disp_max = -1000000000.0; disp_min = 1000000000.0;
	disp_list.clear();

	for(int i = 0; i < box_ti->r_number_particles(); i++) {

		for(int j = 0; j < box_ti->r_dim(); j++) {

			dr = box_ti->r_particle(i)->r_coordinate_true(j) - box_t0->r_particle(i)->r_coordinate_true(j);

			if(com_flag) dr = dr - (box_ti->r_centre_of_mass(j) - box_t0->r_centre_of_mass(j));
			
			if(dr > disp_max) disp_max = dr;
			if(dr < disp_min) disp_min = dr;

			disp_list.push_back(dr);

		}

	}

	bin_size = (disp_max - disp_min)/number_bins;

	disp_count.clear();
	disp_count.resize(number_bins + 1, 0.0);

	for(int i = 0; i < disp_list.size(); i++) {

		bin_no = floor((disp_list[i] - disp_min)/bin_size);
		disp_count[bin_no] += 1.0;

	}

}

void sys_info__mean_squared_displacement::print_displacement_distribution() {

	double time;

	time = box_ti->r_timestamp() - box_t0->r_timestamp();

	for(int i = 0; i < disp_count.size(); i++) if(disp_count[i] > 0) printf("%f %f %f %f %d\n", time, disp_min + i*bin_size, disp_count[i], bin_size, disp_list.size());

	printf("\n\n");

}
