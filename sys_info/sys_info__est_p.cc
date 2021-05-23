#include "sys_info__est_p.h"

//NB When I am confused about this at a future date - rho in the paper where this method comes from (eppenga & frenkel - monte carlo study
//of the isotropic and nematic phases of infinitely thin hard platelets) is defined as rho = (N/V)sigma^3.

void sys_info__est_p::box_in(box *i_box_1) {box_1 = i_box_1;}

double sys_info__est_p::smallest_scaled_separation(int part_1) {

	double smallest_scaled_separation, r_2, s_part_1, s_i, scaled_separation;
	vector <int> neighbour_list;

	smallest_scaled_separation = 1.0;
	for(int j = 0; j < box_1->r_dim(); j++) smallest_scaled_separation += box_1->r_box_size(j)*box_1->r_box_size(j);

	//printf("%f ", smallest_scaled_separation);
	
	s_part_1 = box_1->r_particle(part_1)->r_size();

	neighbour_list = box_1->r_particle_closer_than_list(part_1, 10.0);

	for(int i = 0; i < neighbour_list.size(); i++) {

		r_2 = box_1->metric(part_1, neighbour_list[i]);
		s_i = box_1->r_particle(neighbour_list[i])->r_size();

		scaled_separation = sqrt(r_2)/((s_part_1+s_i)/2.0);

		//printf("%f ", scaled_separation);

		if(scaled_separation < smallest_scaled_separation) smallest_scaled_separation = scaled_separation;

	}

	//printf("%f\n", smallest_scaled_separation);

	return smallest_scaled_separation;
}

double sys_info__est_p::delta_rho_until_overlap(int part_1) {

	double rho_old, rho_new, delta_rho, scale_factor_to_overlap;

	rho_old = box_1->r_number_density()*box_1->r_ave_sigma_3();

	scale_factor_to_overlap = 1.0/smallest_scaled_separation(part_1);

	rho_new = box_1->r_number_particles()*box_1->r_ave_sigma_3();
	for(int j = 0; j < box_1->r_dim(); j++) rho_new = rho_new/(box_1->r_box_size(j) * scale_factor_to_overlap);

	delta_rho = rho_new - rho_old;

	return delta_rho;
}

void sys_info__est_p::calc_delta_rho_distro(double i_bin_size) {

	double delta_rho_calc;
	int number_bins, bin_no;

	//printf("STARTING RHO DISTRO\n");

	delta_rho_list.clear();

	delta_rho_min = 1.0; delta_rho_max = 0.0;

	for(int i = 0; i < box_1->r_number_particles(); i++) {

		//printf("%d ", i);

		delta_rho_calc = delta_rho_until_overlap(i);

		//printf("%f \n", delta_rho_calc);

		delta_rho_list.push_back(delta_rho_calc);
		
		if(delta_rho_calc > delta_rho_max) delta_rho_max = delta_rho_calc;
		if(delta_rho_calc < delta_rho_min) delta_rho_min = delta_rho_calc;

	}

	number_bins = floor((delta_rho_max-delta_rho_min)/i_bin_size);
	bin_size = (delta_rho_max-delta_rho_min)/number_bins;	

	delta_rho_count.clear();
	delta_rho_count.resize(number_bins + 1, 0.0);

	for(int i = 0; i < delta_rho_list.size(); i++) {

		bin_no = floor(delta_rho_list[i]/bin_size);
		delta_rho_count[bin_no] += 1.0;

	}

	for(int i = 0; i < delta_rho_count.size(); i++) delta_rho_count[i] = delta_rho_count[i]/bin_size/delta_rho_list.size();

	//printf("CALCULATED RHO DISTRO!\n");
}

void sys_info__est_p::print_delta_rho_distro() {

	for(int i = 0; i < delta_rho_count.size(); i++) if(delta_rho_count[i] > 0.0) printf("%f %f %f %f %d\n", box_1->r_timestamp(), delta_rho_min + i*bin_size, delta_rho_count[i], bin_size, delta_rho_list.size());
	printf("\n\n");
}

double sys_info__est_p::fit_rho_distro(double i_bin_size) {

	double x_acc = 0.0, x_2_acc = 0.0, y_acc = 0.0, xy_acc = 0.0, count_acc = 0.0, alpha, beta;

	if(!rho_distro_exists) calc_delta_rho_distro(i_bin_size);

	for(int i = 0; i < delta_rho_count.size(); i++) {

		if(delta_rho_count[i] > 0) {

			if(i*bin_size < 0.02) {

				x_acc += delta_rho_min + i*bin_size;
				x_2_acc += (delta_rho_min+i*bin_size)*(delta_rho_min+i*bin_size);
				y_acc += log(delta_rho_count[i]);
				xy_acc += (delta_rho_min + i*bin_size)*log(delta_rho_count[i]);
				count_acc += 1.0;
			}

		}

	}

	alpha = (count_acc*xy_acc - y_acc*x_acc)/(count_acc*x_2_acc - x_acc*x_acc);

	beta = (y_acc - x_acc)/count_acc;

	return alpha;

}

double sys_info__est_p::fit_rho_distro_II(double i_bin_size, double threshold) {

	double x_acc = 0.0, x_2_acc = 0.0, y_acc = 0.0, xy_acc = 0.0, alpha_min, alpha_max, alpha_middle, function_value;
	
	if(!rho_distro_exists) calc_delta_rho_distro(i_bin_size);

	for(int i = 0; i < delta_rho_count.size(); i++) {

		if(delta_rho_count[i] > 0) {

			x_acc += i*bin_size;
			x_2_acc += i*bin_size*i*bin_size;
			y_acc += log(delta_rho_count[i]);
			xy_acc += i*bin_size*log(delta_rho_count[i]);

		}

	}

	alpha_max = 100.0, alpha_min = 0.0; alpha_middle = (alpha_min + alpha_max)/2.0;

	function_value = function_to_zero(x_acc, x_2_acc, y_acc, xy_acc, alpha_middle);

	//printf("%f %f %f %f %f func: %f\n", alpha_middle, x_acc, x_2_acc, y_acc, xy_acc, function_value);

	while(abs(function_value) > threshold) {

		//printf("%f %f %f %f %f func: %f\n", alpha_middle, x_acc, x_2_acc, y_acc, xy_acc, function_value);

		if(function_value < 0.0) alpha_min = alpha_middle;
		else if (function_value > 0.0) alpha_max = alpha_middle;

		alpha_middle = (alpha_min + alpha_max)/2.0;

		function_value = function_to_zero(x_acc, x_2_acc, y_acc, xy_acc, alpha_middle);

		//printf("func: %f\n", function_value);

	}

	return alpha_middle;

}

double sys_info__est_p::function_to_zero(double i_x_acc, double i_x_2_acc, double i_y_acc, double i_xy_acc, double alpha_test) {

	double return_value;

	return_value = -i_y_acc + log(alpha_test)*(double(delta_rho_count.size())-alpha_test*i_x_acc) + alpha_test*alpha_test*i_x_2_acc + alpha_test*(i_xy_acc - i_x_acc);

	return return_value;

}

double sys_info__est_p::r_est_p(double i_bin_size) {

	double alpha_est, pressure_est, rho, kT, ave_sigma_3;
	
	alpha_est = -fit_rho_distro(i_bin_size); //NB careful with the sign of of alpha! v I function returns negative of what is in frenkel paper...

	rho = box_1->r_number_density();

	ave_sigma_3 = box_1->r_ave_sigma_3();

	kT = box_1->r_temperature();

	pressure_est = rho*kT*(1.0 + rho*alpha_est*ave_sigma_3/2.0);

	return pressure_est;

}
