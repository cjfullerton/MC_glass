#ifndef SYSINFO_CALCZ_H
#define SYSINFO_CALCZ_H

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include "sys_info.h"

#define PI 3.141592654

using namespace std;

class sys_info__calc_z : public sys_info {

	private:
		box *box_1;
		double count;
		double bin_size_aim;
		double bin_size;
		double ave_diam;
		double ave_diam_3;
		int no_of_bins;

		vector <double> nr;
		vector <double> nr_r;

	public:
		sys_info__calc_z(double i_bin_size_aim) {

			bin_size_aim = i_bin_size_aim;

		}
		~sys_info__calc_z(){;}

		void box_in(box *i_box_1);
		void setup_arrays();
		void calc_nr();
		void calc_ssf_I(double q);
		void calc_ssf_II(double q);
		void dump_nr();
		void dump_z();
		double r_z();
		void dump_integral_gr_under_first_peak();
};

#endif
