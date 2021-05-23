#include "integrator_MC_NPT.h"


double integrator_MC_NPT::run_dynamics(double moves) { return run_dynamics_particle_volume(moves); }

double integrator_MC_NPT::run_dynamics_particle_volume(double moves) {

	double move_acc = 0.0;

	while(move_acc < moves) {
		if(rand_1->prob_gen() <= prob_v_move) {
			successful_volume_moves += volume_move();
			attempted_volume_moves += 1.0;
		}
		else {
			successful_particle_moves += particle_move();
			attempted_particle_moves += 1.0;
		}

		move_acc += 1.0;
	}

	box_1->c_timestamp(box_1->r_timestamp() + move_acc/double(box_1->r_number_particles()));

	successful_moves = successful_particle_moves + successful_volume_moves;
	attempted_moves = attempted_particle_moves + attempted_volume_moves;

	return successful_moves;
}

double integrator_MC_NPT::run_dynamics_particle(double moves) {

	double move_acc = 0.0;

	while(move_acc < moves) {
		successful_moves += particle_move();
		move_acc += 1.0;
	}

	box_1->c_timestamp(box_1->r_timestamp() + move_acc/double(box_1->r_number_particles()));

	return successful_moves;

}

double integrator_MC_NPT::run_dynamics_volume(double moves) {

	double move_acc = 0.0;

	while(move_acc < moves) {
		successful_moves += volume_move();
		move_acc += 1.0;
	}

	box_1->c_timestamp(box_1->r_timestamp() + move_acc/double(box_1->r_number_particles()));
	
	return successful_moves;

}

void integrator_MC_NPT::dump_moves() { printf("#MC_MOVES: %f %f %f %f %f %f\n", attempted_moves, successful_moves, attempted_particle_moves, successful_particle_moves, attempted_volume_moves, successful_volume_moves); }

void integrator_MC_NPT::reset_moves() {

	successful_moves = 0.0; attempted_moves = 0.0;
	successful_particle_moves = 0.0; attempted_particle_moves = 0.0;
	successful_volume_moves = 0.0; attempted_volume_moves = 0.0;

}

double integrator_MC_NPT::particle_move() {

	int k, accept_move = 0;
	double dr, mc_move_energy, mc_move_energy_new, dE, Prob;

	k=rand_1->prob_gen()*(box_1->r_number_particles());

	if(k==box_1->r_number_particles()) k = k - 1;
		
	mc_move_energy = box_1->particle_energy(k);

	for(int j=0; j<box_1->r_dim(); j++) {\
		displacement[j]=MC_dr*0.5*(2.0*rand_1->prob_gen()-1.0);
		rev_displacement[j]= -displacement[j];
	}

	box_1->displace_particle(k, displacement);
        	
	mc_move_energy_new = box_1->particle_energy(k);
	dE=mc_move_energy_new-mc_move_energy;

	if(mc_move_energy != mc_move_energy) { box_1->displace_particle(k, rev_displacement); } 
	else if(dE <= 0.0) { accept_move = 1;}
	else{
		Prob = exp(-dE/box_1->r_temperature());

		if(rand_1->prob_gen()<=Prob){ accept_move=1; }
		else { box_1->displace_particle(k, rev_displacement); }
	}

	return double(accept_move);

}

double integrator_MC_NPT::volume_move() {

	double mc_move_energy, mc_move_energy_new;
	double dV, V, scale_factor, dH, Prob;
	int accept_move = 0;
		
	mc_move_energy = box_1->total_energy_calc();

	dV = MC_dV*0.5*(2.0*rand_1->prob_gen()-1.0);
	V = box_1->r_volume();

	scale_factor = pow((V+dV)/V,1.0/3.0);
		
	box_1->rescale_all(scale_factor);

	mc_move_energy_new = box_1->total_energy_calc();

	dH = (mc_move_energy_new - mc_move_energy) 
		+ (box_1->r_pressure())*dV 
		- (box_1->r_temperature())*(box_1->r_number_particles())*log((V+dV)/V);

	//printf("dH: %f\n", dH);

	//accept/reject move
		
	if(dH !=dH) { box_1->rescale_all(1.0/scale_factor); }
	else if(dH <= 0.0) {accept_move = 1;}
	else {

		Prob = exp(-dH/box_1->r_temperature());
			
		if(rand_1->prob_gen()<=Prob){ accept_move = 1; }
		else { box_1->rescale_all(1.0/scale_factor); }

	}

	return double(accept_move);

}
