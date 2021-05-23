#ifndef SYSINFO_SCATTFUNC_H
#define SYSINFO_SCATTFUNC_H

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include "sys_info.h"

#define PI 3.141592654

using namespace std;

class sys_info__scattering_functions : public sys_info {

	private:
		box *box_t0;
		box *box_ti;

		vector < double > overlap_per_particle;
		vector < vector < double > > overlap_coarse_grained;

		

	public:
		sys_info__scattering_functions() {;}
		~sys_info__scattering_functions() {;}

		void box_t0_in(box *i_box_t0);
		void box_ti_in(box *i_box_ti);

		void calc_fkt(double q, int com_flag, int version_flag);
		void calc_fkt_II(double q, int com_flag, int version_flag);
		void calc_overlap(double a, int com_flag);
		void calc_overlap_per_particle(double a, int com_flag);
		void coarse_grain_overlap_per_particle(double a, int com_flag, double rcg);
		void dump_number_particles_with_overlap_less_than(double overlap_cut, double a_in, int com_flag_in, double rcg_in);
		double return_number_particles_with_overlap_less_than(double overlap_cut, double a_in, int com_flag_in, double rcg_in);
		void dump_config_with_coarse_grained_overlap(double a, int com_flag, double rcg);
		void dump_config_with_coarse_grained_overlap_II(double a, int com_flag, double rcg);

};
#endif
