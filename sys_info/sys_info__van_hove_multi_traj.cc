#include "sys_info__van_hove_multi_traj.h"

void sys_info__van_hove_multi_traj::box_t0_add(box *i_box_t0) { box_t0.push_back(i_box_t0); }

void sys_info__van_hove_multi_traj::box_ti_add(box *i_box_ti) { box_ti.push_back(i_box_ti); }

void sys_info__van_hove_multi_traj::calc_displacement_distribution_multi_traj(int number_bins, int com_flag) {

	double dr;
	int bin_no;

	if(com_flag) {
		for(int box_pair = 0; box_pair < box_ti.size(); box_pair++ ) {
			box_ti[box_pair]->calc_centre_of_mass(1);
			box_t0[box_pair]->calc_centre_of_mass(1);
		}
	}


	disp_max = -1000000000.0; disp_min = 1000000000.0;
	disp_list.clear();

	for(int box_pair = 0; box_pair < box_ti.size(); box_pair++) {
		for(int i = 0; i < box_ti[box_pair]->r_number_particles(); i++) {
	
			for(int j = 0; j < box_ti[box_pair]->r_dim(); j++) {
	
				dr = box_ti[box_pair]->r_particle(i)->r_coordinate_true(j) - box_t0[box_pair]->r_particle(i)->r_coordinate_true(j);
	
				if(com_flag) dr = dr - (box_ti[box_pair]->r_centre_of_mass(j) - box_t0[box_pair]->r_centre_of_mass(j));
				
				if(dr > disp_max) disp_max = dr;
				if(dr < disp_min) disp_min = dr;
	
				disp_list.push_back(dr);
	
			}
	
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

void sys_info__van_hove_multi_traj::print_displacement_distribution_multi_traj() {

	double time;

	time = box_ti[0]->r_timestamp() - box_t0[0]->r_timestamp();

	for(int i = 0; i < disp_count.size(); i++) if(disp_count[i] > 0) printf("%f %f %f %f %d\n", time, disp_min + i*bin_size, disp_count[i], bin_size, disp_list.size());

	printf("\n\n");

}
