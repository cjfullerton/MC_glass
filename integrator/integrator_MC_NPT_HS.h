#ifndef INTEGRATOR_MC_NPT_HS_H
#define INTEGRATOR_MC_NPT_HS_H

#include "integrator.h"
#include "box__neighbour_list__cell_list__HS.h"
#include "JRand_dyn.h"

#include <cmath>
#include <cstdio>

//SPECIAL FASTER MONTE CARLO INTEGRATOR DESIGNED FOR HARD SPHERE SYSTEM WITH HARD INTERACTIONS ONLY

using namespace std;

class integrator_MC_NPT_HS : public integrator {

	protected:
		box__neighbour_list__cell_list__HS *box_1;
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
		vector<int> prev_smallest_separated_pair;
	

	public:
		integrator_MC_NPT_HS(box__neighbour_list__cell_list__HS *i_box_1, double i_MC_dr, double i_MC_dV, JRand_dyn *i_rand_1) {

			MC_dr = i_MC_dr;
			MC_dV = i_MC_dV;

			box_1 = i_box_1;
			rand_1 = i_rand_1;

			displacement.resize(box_1->r_dim(),0.0);
			rev_displacement.resize(box_1->r_dim(),0.0);

			prev_smallest_separated_pair.resize(2,0);

			successful_moves = 0.0; successful_particle_moves = 0.0; successful_volume_moves = 0.0;
			attempted_moves = 0.0; attempted_particle_moves = 0.0; attempted_volume_moves = 0.0;

			prob_v_move = 1.0/box_1->r_number_particles();

		}

		~integrator_MC_NPT_HS() {;}

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
