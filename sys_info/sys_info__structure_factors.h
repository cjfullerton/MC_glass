#ifndef SYSINFO_STRUCTFACT_H
#define SYSINFO_STRUCTFACT_H

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include "sys_info.h"

#define PI 3.141592654

using namespace std;

class sys_info__structure_factors : public sys_info {

	private:
		box *box_t0;
		box *box_ti;

	public:
		sys_info__structure_factors() {;}
		~sys_info__structure_factors() {;}

		void box_t0_in(box *i_box_t0);
		void box_ti_in(box *i_box_ti);

		void calc_static_structure_factor(double bin_size_aim);

};
#endif
