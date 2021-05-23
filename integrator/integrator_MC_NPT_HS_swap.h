#ifndef INTEGRATOR_MC_NPT_HS_SWAP_H
#define INTEGRATOR_MC_NPT_HS_SWAP_H

#include "integrator.h"
#include "integrator_MC_NPT_HS.h"
#include "box__neighbour_list__cell_list__HS.h"
#include "JRand_dyn.h"

#include <cmath>
#include <cstdio>

//SWAP MOVE TO GO WITH NPT HS MC INTEGRATOR

using namespace std;

class integrator_MC_NPT_HS_swap : private integrator_MC_NPT_HS {

	private:
		double prob_swap_move;
		double attempted_swap_moves;
		double successful_swap_moves;
		double delta_size_cut_2;

	public:
		integrator_MC_NPT_HS_swap(box__neighbour_list__cell_list__HS *i_box_1, double i_MC_dr, double i_MC_dV, double i_delta_size_cut, double i_prob_swap_move, JRand_dyn *i_rand_1) : integrator_MC_NPT_HS(i_box_1, i_MC_dr, i_MC_dV, i_rand_1) {

			prob_swap_move = i_prob_swap_move;
			delta_size_cut_2 = i_delta_size_cut*i_delta_size_cut;

			attempted_swap_moves = 0.0;
			successful_swap_moves = 0.0;
		}

		~integrator_MC_NPT_HS_swap() {;}

		double run_dynamics(double moves);
		double run_dynamics_particle_volume_swap(double moves);
		double run_dynamics_particle_swap(double moves);
		double swap_move();
		void reset_moves();
		void dump_moves();
};

#endif
