#include "integrator_MC_NVT_HS.h"


double integrator_MC_NVT_HS::run_dynamics(double moves) {
    
	int k;
	vector<double> displacement;
	vector<double> rev_displacement;
	double dr, mc_move_energy, mc_move_energy_new, dE, Prob;
	displacement.resize(box_1->r_dim(),0.0);
	rev_displacement.resize(box_1->r_dim(),0.0);
	move_count = 0.0;
	accept_count = 0.0;
    
	for(int l=0; l<(moves); ++l){
        
		k=rand_1->prob_gen()*(box_1->r_number_particles());

		if(k==box_1->r_number_particles()) k = k - 1;
		
		mc_move_energy = box_1->particle_energy(k);

		for(int j=0; j<box_1->r_dim(); j++) {
				displacement[j]=MC_dr*0.5*(2.0*rand_1->prob_gen()-1.0);
				rev_displacement[j]= -displacement[j];
		}

		box_1->displace_particle(k, displacement);
        	
		mc_move_energy_new = box_1->particle_energy(k);
		dE=mc_move_energy_new-mc_move_energy;

		if(dE <= 0.0) accept_count = accept_count + 1.0;
		else box_1->displace_particle(k, rev_displacement);

		move_count = move_count + 1.0;

	}

	MC_sweep_count = MC_sweep_count + move_count/box_1->r_number_particles();

	return accept_count;

}
