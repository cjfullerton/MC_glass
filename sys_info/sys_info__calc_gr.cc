#include "sys_info__calc_gr.h"

void sys_info__calc_gr::box_in(box *i_box_1) { box_1 = i_box_1; }

void sys_info__calc_gr::calc_nr() {

	double data_min, data_max, r, r_r, sigma_ij;
	int number_bins, bin_no, bin_no_r;

	ave_diam = 0.0;
	ave_diam_3 = 0.0;

	data_min = 0.0; data_max = box_1->r_box_size(0)/2.0;

	number_bins = floor((data_max-data_min)/bin_size_aim);
	bin_size = (data_max-data_min)/number_bins;

	nr.clear(); nr_r.clear();
	nr.resize(number_bins,0.0); nr_r.resize(number_bins,0.0);


	for(int i = 0; i < box_1->r_number_particles(); i++) {

		ave_diam = ave_diam + box_1->r_particle(i)->r_size();
		ave_diam_3  = ave_diam_3 + pow(box_1->r_particle(i)->r_size(),3.0);

		for(int j = i+1; j < box_1->r_number_particles(); j++) {

				sigma_ij = box_1->r_particle_interaction->interaction_width(box_1->r_particle(i)->r_size(), box_1->r_particle(j)->r_size());
	
				r = sqrt(box_1->metric(i,j));
				r_r = sqrt(box_1->metric(i,j))/sigma_ij;
	
				if(r < data_max){
					bin_no = floor(r/bin_size);
					nr[bin_no] += 2.0;
				}
				if(r_r < data_max){
					bin_no_r = floor(r_r/bin_size);
					nr_r[bin_no_r] += 2.0;
				}
		}
	}

	ave_diam = ave_diam/box_1->r_number_particles();
	ave_diam_3 = ave_diam_3/box_1->r_number_particles();

}

void sys_info__calc_gr::dump_nr() {

	for(int i = 0; i < nr.size(); i++) printf("%f %f %f %f %f %f %f %f %f %f\n", box_1->r_timestamp(), i*bin_size, nr[i], nr_r[i], 4.0*PI*(i*bin_size)*(i*bin_size)*bin_size, box_1->r_number_density(), box_1->r_packing_fraction(), ave_diam, ave_diam_3, box_1->r_pressure());

}

