#ifndef SYSINFO_ESTP_H
#define SYSINFO_ESTP_H

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include "sys_info.h"

#define PI 3.141592654

using namespace std;

class sys_info__est_p : public sys_info {\
	
	private:
		box* box_1;
		vector <double> delta_rho_list;
		vector <double> delta_rho_count;

		double bin_size;
		double rho_delta_min;
		double rho_delta_min;
		int rho_distro_exists;

	public:
		sys_info__est_p() {
			rho_distro_exists = 0;
		}
		~sys_info__est_p() {;}

		void box_in(box *i_box_1);
		double smallest_scaled_separation(int part_1);
		double delta_rho_until_overlap(int part_1);
		void calc_delta_rho_distro(double i_bin_size);
		void print_delta_rho_distro();
		double fit_rho_distro(double i_bin_size);
		double fit_rho_distro_II(double i_bin_size, double threshold);
		double function_to_zero(double i_x_acc, double i_x_2_acc, double i_y_acc, double i_xy_acc, double alpha_test);
		double r_est_p(double i_bin_size);

};

#endif
