#include "sys_info__calc_z.h"

void sys_info__calc_z::box_in(box *i_box_1) { box_1 = i_box_1; }

void sys_info__calc_z::calc_nr() {

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

}

void sys_info__calc_z::dump_nr() {

	for(int i = 0; i < nr.size(); i++) printf("#gr %f %f %f %f %f %f %f %f %f %f\n", box_1->r_timestamp(), i*bin_size, nr[i], nr_r[i], 4.0*PI*(i*bin_size)*(i*bin_size)*bin_size, box_1->r_number_density(), box_1->r_packing_fraction(), ave_diam, ave_diam_3, box_1->r_pressure());

}

void sys_info__calc_z::calc_ssf_I(double q) {

	double int_acc = 0.0, r, gr_i;

	for(int i = 1; i < nr_r.size(); i++) {

 		gr_i = nr_r[i]/(4.0*PI*(i*bin_size)*(i*bin_size)*bin_size)/box_1->r_number_particles()/(box_1->r_number_density()*ave_diam_3);

		r = i*bin_size;

		int_acc += r*sin(q*r)*(gr_i-1.0)*bin_size;

	}

	printf("%f %f %f\n", box_1->r_timestamp(), q, 1.0 + 4.0*PI*box_1->r_number_density()/q*int_acc);

}

void sys_info__calc_z::calc_ssf_II(double q) {

	double int_acc = 0.0, r, gr_i;

	for(int i = 1; i < nr_r.size(); i++) {

 		gr_i = nr[i]/(4.0*PI*(i*bin_size)*(i*bin_size)*bin_size)/box_1->r_number_particles();

		r = i*bin_size;

		int_acc += r*sin(q*r)*(gr_i-1.0)*bin_size;

	}

	printf("%f %f %f\n", box_1->r_timestamp(), q, 1.0 + 4.0*PI*box_1->r_number_density()/q*int_acc);

}

void sys_info__calc_z::dump_z() {

	double phi, rho, sigma, gr_contact, Z;

	rho = box_1->r_number_density();

	phi = box_1->r_packing_fraction();

	gr_contact = nr[floor(1.0/bin_size)]/box_1->r_number_particles()/(4.0*PI*1.0*1.0*bin_size)/(rho*ave_diam_3);

	Z = 1.0 + 4.0*phi*gr_contact;

	printf("%f %f\n", phi, Z);

}

double sys_info__calc_z::r_z() {

	double phi, rho, sigma, gr_contact, Z;

	rho = box_1->r_number_density();

	phi = box_1->r_packing_fraction();

	gr_contact = nr[floor(1.0/bin_size)]/box_1->r_number_particles()/(4.0*PI*1.0*1.0*bin_size)/(rho*ave_diam_3);

	Z = 1.0 + 4.0*phi*gr_contact;

	return Z;

}

void sys_info__calc_z::dump_integral_gr_under_first_peak() {

	double gr_integral_acc = 0.0, r, rho, approx_z_integral_acc = 0.0;
	int bin_iterator;

	calc_nr();

	bin_iterator = 1;

	rho = box_1->r_number_density();

	while(r < 1.5) {
		
		r = bin_iterator*bin_size;
		gr_integral_acc += nr_r[bin_iterator]/box_1->r_number_particles()/(4.0*PI*r*r*bin_size)/(rho*ave_diam_3)*bin_size;
		approx_z_integral_acc += (1.0 + 4.0*box_1->r_packing_fraction()*nr_r[bin_iterator]/box_1->r_number_particles()/(4.0*PI*r*r*bin_size)/(rho*ave_diam_3))*bin_size;
		bin_iterator++;

	}

	printf("#gr_int %f %f %f %f\n", box_1->r_timestamp(), gr_integral_acc, approx_z_integral_acc, box_1->r_packing_fraction());

}
