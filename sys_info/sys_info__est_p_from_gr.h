#ifndef SYSINFO_ESTPFROMGR_H
#define SYSINFO_ESTPFROMGR_H

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include "sys_info.h"

#define PI 3.141592654

using namespace std;

class sys_info__est_p_from_gr : public sys_info {

	private:
		box *box_1;
		double count;
		double bin_size;
		double data_min;
		double data_max;
		double ave_diam;
		double ave_diam_3;
		double alpha;
		double beta;
		int no_of_bins;

		vector <double> nr;
		vector <double> nr_r;
		vector <double> gr_r;

	public:
		sys_info__est_p_from_gr() {;}
		~sys_info__est_p_from_gr(){;}

		void box_in(box *i_box_1);
		void setup_arrays();
		void calc_gr(double i_bin_size_aim);
		void fit_gr(double i_bin_size_aim, double lower_limit, double upper_limit);
		double estimate_pressure(double i_bin_size_aim, double lower_limit, double upper_limit);
};

#endif
