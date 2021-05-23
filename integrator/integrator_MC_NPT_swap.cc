#include "integrator_MC_NPT_swap.h"

double integrator_MC_NPT_swap::run_dynamics(double moves) { return run_dynamics_particle_volume_swap(moves); }

double integrator_MC_NPT_swap::run_dynamics_particle_volume_swap(double moves) {

	double rand_number_for_move_type, move_acc = 0.0;

	while(move_acc < moves) {

		rand_number_for_move_type = rand_1->prob_gen();

		if(rand_number_for_move_type <= prob_v_move) {
			successful_volume_moves += volume_move();
			attempted_volume_moves += 1.0;
		//	printf("VOLUME: %f/%f\n", successful_volume_moves, attempted_volume_moves);
		}
		else if(rand_number_for_move_type <= prob_v_move + prob_swap_move) {
			successful_swap_moves += swap_move();
			attempted_swap_moves += 1.0;
		//	printf("SWAP: %f/%f\n", successful_swap_moves, attempted_swap_moves);
		}
		else {
			successful_particle_moves += particle_move();
			attempted_particle_moves += 1.0;
		//	printf("PARTICLE: %f/%f\n", successful_particle_moves, attempted_particle_moves);
		}

		move_acc += 1.0;
	}

	box_1->c_timestamp(box_1->r_timestamp() + move_acc/double(box_1->r_number_particles()));

	successful_moves = successful_particle_moves + successful_volume_moves + successful_swap_moves;
	attempted_moves = attempted_particle_moves + attempted_volume_moves + attempted_swap_moves;;

	return successful_moves;
}

double integrator_MC_NPT_swap::run_dynamics_particle_swap(double moves) {

	double move_acc = 0.0;

	while(move_acc < moves) {
		if(rand_1->prob_gen() <= prob_swap_move) {
			successful_swap_moves += swap_move();
			attempted_swap_moves += 1.0;
		}
		else {
			successful_particle_moves += particle_move();
			attempted_particle_moves += 1.0;
		}

		move_acc += 1.0;
	}

	box_1->c_timestamp(box_1->r_timestamp() + move_acc/double(box_1->r_number_particles()));

	successful_moves = successful_particle_moves + successful_swap_moves;
	attempted_moves = attempted_particle_moves + attempted_swap_moves;

	return successful_moves;
}

double integrator_MC_NPT_swap::swap_move() {

	int part_1, part_2;
	double mc_move_energy_new, mc_move_energy_old, dE, Prob, size_part_1_old, size_part_2_old, delta_size, delta_size_2, accept_move = 0.0;

	part_1=rand_1->prob_gen()*(box_1->r_number_particles());
	if(part_1==box_1->r_number_particles()) part_1 = part_1 - 1;

	part_2=rand_1->prob_gen()*(box_1->r_number_particles());
	if(part_2==box_1->r_number_particles()) part_2 = part_2 - 1;
	
	mc_move_energy_old = box_1->particle_energy(part_1) + box_1->particle_energy(part_2);
	
	size_part_1_old = box_1->r_particle(part_1)->r_size();
	size_part_2_old = box_1->r_particle(part_2)->r_size();

	delta_size = size_part_1_old - size_part_2_old;
	delta_size_2 = delta_size*delta_size;

	if(delta_size_2 >= delta_size_cut_2) return accept_move;

	box_1->r_particle(part_1)->c_size(size_part_2_old);
	box_1->r_particle(part_2)->c_size(size_part_1_old);

	mc_move_energy_new = box_1->particle_energy(part_1) + box_1->particle_energy(part_2);
	dE=mc_move_energy_new-mc_move_energy_old;

	if(mc_move_energy_new != mc_move_energy_new) {
		box_1->r_particle(part_1)->c_size(size_part_1_old);
		box_1->r_particle(part_2)->c_size(size_part_2_old);
	}
	else if(dE <= 0.0) { 
		accept_move = 1;
		box_1->check_particle_nlist_and_update(part_1);
		box_1->check_particle_nlist_and_update(part_2);
	}
	else {
		Prob = exp(-dE/box_1->r_temperature());

		if(rand_1->prob_gen()<=Prob){ 
			accept_move = 1;
			box_1->check_particle_nlist_and_update(part_1);
			box_1->check_particle_nlist_and_update(part_2);
		}

		else { 
			box_1->r_particle(part_1)->c_size(size_part_1_old);
			box_1->r_particle(part_2)->c_size(size_part_2_old);
		}
	}

	return accept_move;

}

void integrator_MC_NPT_swap::reset_moves() {

	successful_moves = 0.0; attempted_moves = 0.0;
	successful_particle_moves = 0.0; attempted_particle_moves = 0.0;
	successful_volume_moves = 0.0; attempted_volume_moves = 0.0;
	successful_swap_moves = 0.0; attempted_swap_moves = 0.0;

}

void integrator_MC_NPT_swap::dump_moves() { printf("%f %f %f %f %f %f %f %f\n", attempted_moves, successful_moves, attempted_particle_moves, successful_particle_moves, attempted_volume_moves, successful_volume_moves, attempted_swap_moves, attempted_particle_moves); }

void integrator_MC_NPT_swap::print_moves(char out_file_name[], double i_time) {

	FILE *acceptance_out;

	acceptance_out = fopen(out_file_name, "a");

	fprintf(acceptance_out, "%f %f %f %f %f %f %f %f %f\n", i_time, attempted_moves, successful_moves, attempted_particle_moves, successful_particle_moves, attempted_volume_moves, successful_volume_moves, attempted_swap_moves, successful_swap_moves);

	fclose(acceptance_out);
}
