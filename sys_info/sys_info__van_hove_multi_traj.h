#ifndef SYSINFO_MSD_H
#define SYSINFO_MSD_H

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include "sys_info.h"

using namespace std;

class sys_info__van_hove_multi_traj : public sys_info {

	private:
		vector <box*> box_t0;
		vector <box*> box_ti;

		double bin_size;
		double disp_max;
		double disp_min;
		vector <double> disp_count;
		vector <double> disp_list;

	public:
		sys_info__van_hove_multi_traj() {;}
		~sys_info__van_hove_multi_traj() {;}

		void box_t0_add(box *i_box_t0);
		void box_ti_add(box *i_box_ti);

		void calc_displacement_distribution_multi_traj(int number_bins, int com_flag);
		void print_displacement_distribution_multi_traj();

};

#endif
