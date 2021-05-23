#ifndef INTEGRATOR_MC_NPT_H
#define INTEGRATOR_MC_NPT_H

#include "integrator.h"
#include "box__neighbour_list__cell_list.h"
#include "JRand_dyn.h"

#include <cmath>
#include <cstdio>

using namespace std;

class integrator_MC_NPT : public integrator {

	protected:
		box__neighbour_list__cell_list *box_1;
		JRand_dyn *rand_1;
		double MC_dr;
		double MC_dV;
		double prob_v_move;
		double attempted_moves;
		double attempted_particle_moves;
		double attempted_volume_moves;
		double successful_moves;
		double successful_particle_moves;
		double successful_volume_moves;
		vector<double> displacement;
		vector<double> rev_displacement;
	
	public:
		integrator_MC_NPT(box__neighbour_list__cell_list *i_box_1, double i_MC_dr, double i_MC_dV, JRand_dyn *i_rand_1) {

			MC_dr = i_MC_dr;
			MC_dV = i_MC_dV;

			box_1 = i_box_1;
			rand_1 = i_rand_1;

			displacement.resize(box_1->r_dim(),0.0);
			rev_displacement.resize(box_1->r_dim(),0.0);

			prob_v_move = 1.0/box_1->r_number_particles();

			successful_moves = 0.0; successful_particle_moves = 0.0; successful_volume_moves = 0.0;
			attempted_moves = 0.0; attempted_particle_moves = 0.0; attempted_volume_moves = 0.0;

		}

		~integrator_MC_NPT() {;}

		double run_dynamics(double moves);
		double run_dynamics_particle_volume(double moves);
		double run_dynamics_particle(double moves);
		double run_dynamics_volume(double moves);
		void reset_moves();
		void dump_moves();
		double particle_move();
		double volume_move();
};

#endif
