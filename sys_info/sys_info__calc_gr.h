#ifndef SYSINFO_CALCGR_H
#define SYSINFO_CALCGR_H

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include "sys_info.h"

#define PI 3.141592654

using namespace std;

class sys_info__calc_gr : public sys_info {

	private:
		box *box_1;
		double bin_size_aim;
		double bin_size;
		double ave_diam;
		double ave_diam_3;
		int no_of_bins;

		vector <double> nr;
		vector <double> nr_r;

	public:
		sys_info__calc_gr(double i_bin_size_aim) {

			bin_size_aim = i_bin_size_aim;

		}
		~sys_info__calc_gr(){;}

		void box_in(box *i_box_1);
		void calc_nr();
		void dump_nr();
};

#endif
