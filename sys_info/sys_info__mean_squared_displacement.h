#ifndef SYSINFO_MSD_H
#define SYSINFO_MSD_H

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include "sys_info.h"

using namespace std;

class sys_info__mean_squared_displacement : public sys_info {

	private:
		box *box_t0;
		box *box_ti;

		double bin_size;
		double dr2_min;
		double disp_max;
		double dr2_max;
		double disp_min;
		vector <double> MSD_count;
		vector <double> MSD_list;
		vector <double> disp_count;
		vector <double> disp_list;

	public:
		sys_info__mean_squared_displacement() {;}
		~sys_info__mean_squared_displacement() {;}

		void box_t0_in(box *i_box_t0);
		void box_ti_in(box *i_box_ti);

		void calc_msd(int com_flag, int version_flag);
		void calc_van_hove_function(int number_bins, int com_flag, int version_flag);
		void print_van_hove_function();
		void calc_displacement_distribution(int number_bins, int com_flag);
		void print_displacement_distribution();

};

#endif
