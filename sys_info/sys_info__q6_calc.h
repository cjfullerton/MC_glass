#ifndef SYSINFO_Q6CALC_H
#define SYSINFO_Q6CALC_H

#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>
#include "sys_info.h"
#include "gsl/gsl_sf_legendre.h"

using namespace std;

class sys_info__q6_calc : public sys_info {
	
	private:
		box *box_ti;
		vector < vector <int> > bond_list;
		vector < vector <int> > bond_array;

		int bond_list_exists;

	public:
		sys_info__q6_calc() {;}
		~sys_info__q6_calc() {;}

		void box_ti_in(box *i_box_ti);
		void q6_cor(double bin, double maxr, int q_r_flag, const char *tag);
		double dxfn(int dim, int part_1, int part_2);
		void create_bond_list();
};

#endif
