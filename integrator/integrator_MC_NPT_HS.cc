#include "integrator_MC_NPT_HS.h"

double integrator_MC_NPT_HS::run_dynamics(double moves) { return run_dynamics_particle_volume(moves); }

double integrator_MC_NPT_HS::run_dynamics_particle_volume(double moves) {

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

double integrator_MC_NPT_HS::run_dynamics_particle(double moves) {

	double move_acc = 0.0;

	while(move_acc < moves) {
	//	printf("!\n");
		successful_moves += particle_move();
		move_acc += 1.0;
	}

	box_1->c_timestamp(box_1->r_timestamp() + move_acc/double(box_1->r_number_particles()));

	return successful_moves;

}

double integrator_MC_NPT_HS::run_dynamics_volume(double moves) {

	double move_acc = 0.0;

	while(move_acc < moves) {
		successful_moves += volume_move();
		move_acc += 1.0;
	}

	box_1->c_timestamp(box_1->r_timestamp() + move_acc/double(box_1->r_number_particles()));
	
	return successful_moves;

}

void integrator_MC_NPT_HS::dump_moves() { printf("#MC_MOVES: %f %f %f %f %f %f\n", attempted_moves, successful_moves, attempted_particle_moves, successful_particle_moves, attempted_volume_moves, successful_volume_moves); }

void integrator_MC_NPT_HS::reset_moves() {

	successful_moves = 0.0; attempted_moves = 0.0;
	successful_particle_moves = 0.0; attempted_particle_moves = 0.0;
	successful_volume_moves = 0.0; attempted_volume_moves = 0.0;

}

double integrator_MC_NPT_HS::particle_move() {

	int k, accept_move = 0;
	double mc_move_energy, prev_smallest_scaled_separation;

	k=rand_1->prob_gen()*(box_1->r_number_particles());

	if(k==box_1->r_number_particles()) k = k - 1;

	prev_smallest_separated_pair[0] = box_1->r_smallest_separated_pair()[0];
	prev_smallest_separated_pair[1] = box_1->r_smallest_separated_pair()[1];

	for(int j=0; j<box_1->r_dim(); j++) {
		displacement[j]=MC_dr*0.5*(2.0*rand_1->prob_gen()-1.0);
		rev_displacement[j]= -displacement[j];
	}

	prev_smallest_scaled_separation = box_1->r_smallest_scaled_separation(); 

	box_1->displace_particle(k, displacement);
        	
	mc_move_energy = box_1->particle_energy(k);

	if(mc_move_energy != mc_move_energy) {
		box_1->displace_particle(k, rev_displacement);
		box_1->c_smallest_scaled_separation(prev_smallest_scaled_separation);
		box_1->c_smallest_separated_pair(0,prev_smallest_separated_pair[0]);
		box_1->c_smallest_separated_pair(1,prev_smallest_separated_pair[1]);
	}
	else {
		accept_move = 1;
		if(k == prev_smallest_separated_pair[0] || k == prev_smallest_separated_pair[1]) box_1->find_smallest_separation();

	}

	return double(accept_move);

}

double integrator_MC_NPT_HS::volume_move() {

	double mc_move_energy, mc_move_energy_new;
	double dV, V, scale_factor, dH, Prob, rand;
	int accept_move = 0; //ASSUME REJECTION

	dV = MC_dV*0.5*(2.0*rand_1->prob_gen()-1.0);
	V = box_1->r_volume();
		
	scale_factor = pow((V+dV)/V,1.0/3.0);

	dH = 0.0;

	if(box_1->r_smallest_scaled_separation()*scale_factor >= 1.0) {

		mc_move_energy_new = 0.0; mc_move_energy = 0.0;

		dH = (mc_move_energy_new - mc_move_energy) 
			+ (box_1->r_pressure())*dV 
			- (box_1->r_temperature())*(box_1->r_number_particles())*log((V+dV)/V);
		
		if(dH <= 0.0) accept_move = 1;
		else {
			Prob = exp(-dH/box_1->r_temperature());
			rand = rand_1->prob_gen();
			if(rand <=  Prob) accept_move = 1; 
		}
		
	}


	if(accept_move) {
		box_1->rescale_all(scale_factor);
	}

	return double(accept_move);

}
