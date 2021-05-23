#include "sys_info__est_p_from_gr.h"

void sys_info__est_p_from_gr::box_in(box *i_box_1) { box_1 = i_box_1; }

void sys_info__est_p_from_gr::calc_gr(double i_bin_size_aim) {

	double r, r_r, sigma_ij;
	int number_bins, bin_no, bin_no_r;

	ave_diam = 0.0;
	ave_diam_3 = 0.0;

	data_min = 0.0; data_max = box_1->r_box_size(0)/2.0;

	number_bins = floor((data_max-data_min)/i_bin_size_aim);
	bin_size = (data_max-data_min)/number_bins;

	nr.clear(); nr_r.clear();
	nr.resize(number_bins,0.0); nr_r.resize(number_bins,0.0);


	for(int i = 0; i < box_1->r_number_particles(); i++) {
		ave_diam = ave_diam + box_1->r_particle(i)->r_size();
		ave_diam_3  = ave_diam_3 + pow(box_1->r_particle(i)->r_size(),3.0);
		for(int j = i+1; j < box_1->r_number_particles(); j++) {

				sigma_ij = (box_1->r_particle(i)->r_size() + box_1->r_particle(j)->r_size())/2.0;
	
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

	gr_r.clear(); gr_r.resize(nr_r.size(),0.0);

	for(int i = 0; i < nr_r.size(); i++) {
		gr_r[i] = nr_r[i]/(4.0*PI*(i*bin_size)*(i*bin_size)*bin_size)/box_1->r_number_particles()/(box_1->r_number_density()*ave_diam_3);
	//	printf("%f %f\n", data_min+i*bin_size, gr_r[i]);
	}

}

void sys_info__est_p_from_gr::fit_gr(double i_bin_size, double lower_limit, double upper_limit) {

	double x_acc = 0.0, x_2_acc = 0.0, y_acc = 0.0, xy_acc = 0.0, count_acc = 0.0;

	calc_gr(i_bin_size);

	for(int i = 0; i < gr_r.size(); i++) {

		if(gr_r[i] > 0) {

			if((data_min + i*bin_size) >= lower_limit && (data_min + i*bin_size) <= upper_limit) {

				//printf("%f %f\n", data_min + i*bin_size, gr_r[i]);

				x_acc += data_min + i*bin_size;
				x_2_acc += (data_min+i*bin_size)*(data_min+i*bin_size);
				y_acc += gr_r[i];
				xy_acc += (data_min + i*bin_size)*gr_r[i];
				count_acc += 1.0;
			}

		}

	}

	alpha = (count_acc*xy_acc - y_acc*x_acc)/(count_acc*x_2_acc - x_acc*x_acc);

	beta = (y_acc - alpha*x_acc)/count_acc;
	
	//printf("x:%f x2: %f y:%f xy:%f\n", x_acc, x_2_acc, y_acc, xy_acc);
	//printf("a: %f b: %f\n", alpha, beta);

}

double sys_info__est_p_from_gr::estimate_pressure(double i_bin_size, double lower_limit, double upper_limit) {

	double estimated_pressure, gr_contact;

	fit_gr(i_bin_size, lower_limit, upper_limit); 

	gr_contact = alpha + beta;

	estimated_pressure = (1.0 + 4.0*box_1->r_packing_fraction()*gr_contact)*box_1->r_number_density();

	return estimated_pressure;

}
