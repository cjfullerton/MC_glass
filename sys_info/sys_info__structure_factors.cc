#include "sys_info__structure_factors.h"

void sys_info__structure_factors::box_t0_in(box *i_box_t0) { box_t0 = i_box_t0;}

void sys_info__structure_factors::box_ti_in(box *i_box_ti) { box_ti = i_box_ti;}

void sys_info__structure_factors::calc_static_structure_factor(double bin_size_aim) {

	vector <double> q;
	vector <vector <double> > s_q;
	vector <vector <double> > binned_s_q;
	vector <double> data_point;
	vector <int> max_q_vectors;
	double cos_q_r_acc = 0.0, sin_q_r_acc = 0.0, q_size, q_min, q_max, bin_size;
	int no_bins, bin_no, c_0_I, c_1_I, c_2_I;

	s_q.clear(); binned_s_q.clear();

	q.resize(box_ti->r_dim(), 0.0);
	max_q_vectors.resize(box_ti->r_dim(), 10);
	data_point.resize(2);

	for(int j = 0; j < box_ti->r_dim(); j++) q[j] = 2.0*PI/box_ti->r_box_size(j);

	for(int c_0 = 0; c_0 < max_q_vectors[0]; c_0++) { 
		for(int c_1 = 0; c_1 < max_q_vectors[1]; c_1++) { 
			for(int c_2 = 0; c_2 < max_q_vectors[2]; c_2++) { 
				
				sin_q_r_acc = 0.0;
				cos_q_r_acc = 0.0;

				for(int i = 0; i < box_ti -> r_number_particles(); i++) {

					cos_q_r_acc += cos(c_0*q[0]*box_ti->r_particle(i)->r_coordinate_box(0) + c_1*q[1]*box_ti->r_particle(i)->r_coordinate_box(1) + c_2*q[2]*box_ti->r_particle(i)->r_coordinate_box(2));
					sin_q_r_acc += sin(c_0*q[0]*box_ti->r_particle(i)->r_coordinate_box(0) + c_1*q[1]*box_ti->r_particle(i)->r_coordinate_box(1) + c_2*q[2]*box_ti->r_particle(i)->r_coordinate_box(2));

				}

				q_size = sqrt(c_0*q[0]*c_0*q[0] + c_1*q[1]*c_1*q[1] + c_2*q[2]*c_2*q[2]);
				data_point[0] = q_size;
				data_point[1] = (cos_q_r_acc*cos_q_r_acc + sin_q_r_acc*sin_q_r_acc)/box_ti->r_number_particles();

				s_q.push_back(data_point);

			}
		}
	}

	c_0_I = max_q_vectors[0]-1;
	c_1_I = max_q_vectors[1]-1;
	c_2_I = max_q_vectors[2]-1;

	q_min = 0.0;

	q_max = sqrt(c_0_I*q[0]*c_0_I*q[0] + c_1_I*q[1]*c_1_I*q[1] + c_2_I*q[2]*c_2_I*q[2]);
	no_bins = floor((q_max-q_min)/bin_size_aim);
	bin_size = (q_max-q_min)/no_bins;

	binned_s_q.resize(no_bins+1);
	for(int i = 0; i < binned_s_q.size(); i++) binned_s_q[i].resize(2, 0.0);

	for(int i = 1; i < s_q.size(); i++) {

		bin_no = floor(s_q[i][0]/bin_size);
		binned_s_q[bin_no][0] += s_q[i][1];
		binned_s_q[bin_no][1] += 1.0;
	}

	for(int i = 0; i < binned_s_q.size(); i++) {

		if(binned_s_q[i][1] > 0.0) printf("%f %f %f\n", box_ti->r_timestamp(), bin_size*i, binned_s_q[i][0]/binned_s_q[i][1]);
		
	}

}
