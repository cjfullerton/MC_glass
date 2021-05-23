#ifndef INTEGRATOR_MC_NVT_HS_H
#define INTEGRATOR_MC_NVT_HS_H

#include "integrator.h"
#include "box.h"
#include "JRand_dyn.h"

#include <cmath>
#include <cstdio>

using namespace std;

class integrator_MC_NVT_HS : public integrator {

	private:
		box *box_1;
		JRand_dyn *rand_1;
		double MC_dr;
		double move_count;
		double MC_sweep_count; 
		double accept_count;
	
	public:
		integrator_MC_NVT_HS(box *i_box_1, double i_MC_dr, JRand_dyn *i_rand_1) {

			MC_dr = i_MC_dr;

			box_1 = i_box_1;
			rand_1 = i_rand_1;

			move_count = 0.0; MC_sweep_count = 0.0; accept_count = 0.0;

		}

		~integrator_MC_NVT_HS() {;}

		double run_dynamics(double moves);

};

#endif
