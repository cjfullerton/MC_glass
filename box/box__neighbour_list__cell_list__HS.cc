#include "box__neighbour_list__cell_list__HS.h"

double box__neighbour_list__cell_list__HS::r_smallest_scaled_separation() { return smallest_scaled_separation; }

void box__neighbour_list__cell_list__HS::c_smallest_scaled_separation(double new_smallest_scaled_separation) { smallest_scaled_separation = new_smallest_scaled_separation; }
//THIS IS AN ODD LOOKING FUNCTION, BUT IT SHOULD BE USED IF A PARTICLE DISPLACEMENT IS REJECTED AND UNDONE (IN CASE THE DISPLACEMENT REDUCED THE SMALLEST SCALED SEPARATION.

vector <int> box__neighbour_list__cell_list__HS::r_smallest_separated_pair() { return smallest_separated_pair; }

void box__neighbour_list__cell_list__HS::c_smallest_separated_pair(int no_in_pair, int particle) { smallest_separated_pair[no_in_pair] = particle; }

double box__neighbour_list__cell_list__HS::calc_scaled_separation(int i_part_1, int i_part_2) {

	double r_ij_2, scaled_separation, int_width;

	r_ij_2 = metric(i_part_1, i_part_2);
	int_width = part_interaction_1->interaction_width(particle_list[i_part_1]->r_size(), particle_list[i_part_2]->r_size());

	scaled_separation = sqrt(r_ij_2)/int_width;

	if(scaled_separation < smallest_scaled_separation) {
		smallest_scaled_separation = scaled_separation;
		smallest_separated_pair[0] = i_part_1; smallest_separated_pair[1] = i_part_2;
	}

	return scaled_separation;
}

double box__neighbour_list__cell_list__HS::calc_scaled_separation(int i_part_1, int i_part_2, double i_rij_2) {

	double r_ij_2, int_width, scaled_separation;

	r_ij_2 = i_rij_2;
	int_width = part_interaction_1 -> interaction_width(particle_list[i_part_1]->r_size(), particle_list[i_part_2]->r_size());

	scaled_separation = sqrt(r_ij_2)/int_width;

	if(scaled_separation < smallest_scaled_separation) {
		smallest_scaled_separation = scaled_separation;
		smallest_separated_pair[0] = i_part_1; smallest_separated_pair[1] = i_part_2;
	}

	return scaled_separation;
}

void box__neighbour_list__cell_list__HS::find_smallest_separation() {

	int cell, neigh, neigh_cell, part_2;

	smallest_scaled_separation = 1.0;
	for(int j = 0; j<dim; j++) smallest_scaled_separation = smallest_scaled_separation + box_size[j]*box_size[j];
	smallest_scaled_separation = sqrt(smallest_scaled_separation);

	smallest_separated_pair.resize(2,0);

	if(nlist_flag && nlist_done) { //Do neighbour list procedure
		for(int part_1 = 0; part_1 < number_particles; part_1++) {
			for(int j = 0; j < nlist[part_1].size(); j++) {
				part_2 = nlist[part_1][j];
				if(part_1 != part_2) calc_scaled_separation(part_1, part_2);
			}
		}
	}
	else if(cells_flag && cells_done) { //Do calc with cells

		for(int part_1; part_1 < number_particles; part_1++) {

			cell = part_cell[part_1][0];

			for(int k = 0; k < dim*dim*dim; k++) {
				neigh_cell = cell_neighbours[cell][k];
				for(int j = 0; j<cell_contents[neigh_cell].size(); j++) {
					part_2 = cell_contents[neigh_cell][j];
					if(part_2 != part_1) calc_scaled_separation(part_1, part_2);
					}
				
				}
			}
		}
	else { //Do normal procedure...
		for(int part_1 = 0; part_1 < number_particles; part_1++) {
			for(int part_2=part_1 + 1; part_2<number_particles; part_2++){
				if(part_1 != part_2) calc_scaled_separation(part_1, part_2);
			}
		}
	}
}

double box__neighbour_list__cell_list__HS::interaction_energy(int part_1, int part_2) {

	double energy, r_ij_2;

	r_ij_2 = metric(part_1, part_2);

	energy = part_interaction_1->interaction_energy(r_ij_2, particle_list[part_1]->r_size(), particle_list[part_2]->r_size());

	calc_scaled_separation(part_1, part_2, r_ij_2);

	return energy;

}

double box__neighbour_list__cell_list__HS::particle_energy(int part_1){
	
	double total_energy=0.0;
	int cell, neigh_cell, part_2;

	if(nlist_flag && nlist_done) { //Do neighbour list procedure
		for(int j = 0; j < nlist[part_1].size(); j++) {
			total_energy = total_energy + interaction_energy(part_1,nlist[part_1][j]);
			if(total_energy != total_energy) break;
		}
	}
	else if(cells_flag && cells_done) { //Do calc with cells

		cell = part_cell[part_1][0];

		for(int k = 0; k < dim*dim*dim; k++) {
			neigh_cell = cell_neighbours[cell][k];
			for(int j = 0; j<cell_contents[neigh_cell].size(); j++) {
				part_2 = cell_contents[neigh_cell][j];
				if(part_2 != part_1) { 
					total_energy = total_energy + interaction_energy(part_1, part_2);
				}
				
			}
			
		}
	}
	else { //Do normal procedure...
		for(int part_2=0; part_2<number_particles; part_2++){
			if( part_2 != part_1) { 
				total_energy = total_energy + interaction_energy(part_1, part_2);
				if(interaction_energy(part_1,part_2) != interaction_energy(part_1,part_2)) printf("OH DEAR %d %d %f %f %f\n", part_1, part_2, metric(part_1, part_2), particle_list[part_1]->r_size(), particle_list[part_2]->r_size());
			}
		}
	}

	total_energy = total_energy + boundary_1->boundary_int_energy(particle_list[part_1]->r_coordinate_box_vec(), particle_list[part_1]->r_size());
	total_energy = total_energy + ext_field_1->field_energy(particle_list[part_1]->r_coordinate_box_vec(), particle_list[part_1]->r_size());
    
	return total_energy;

}

void box__neighbour_list__cell_list__HS::rescale_all(double scale_factor) {
	
	double old_density, dr;
	int rebuild_nlist_flag;

	old_density = r_number_density();
	rebuild_nlist_flag = 0;
	
	for(int i = 0; i < number_particles; i++) {

		if(nlist_flag) dr = 0.0;
		for(int j = 0; j < dim; j++) {

			if(nlist_flag) {
				disp_since_rebuild[i][j] = disp_since_rebuild[i][j]/scale_factor;
				dr = dr + disp_since_rebuild[i][j]*disp_since_rebuild[i][j];
			}

			particle_list[i]->c_coordinate_box(j, scale_factor*(particle_list[i]->r_coordinate_box(j)));
			particle_list[i]->c_coordinate_affine(j, scale_factor*(particle_list[i]->r_coordinate_affine(j)));
		}
		if(nlist_flag) {
			nlist_rebuild_limit = (OUTER_RANGE - particle_list[i]->r_size())/2.0;
			nlist_rebuild_limit = nlist_rebuild_limit*nlist_rebuild_limit;
			if(dr > nlist_rebuild_limit) rebuild_nlist_flag = 1;
		}
	}

	for(int j = 0; j < dim; j++) box_size[j] = box_size[j]*scale_factor;

	smallest_scaled_separation = scale_factor*smallest_scaled_separation;

	packing_fraction = r_packing_fraction();

	if(cells_flag) {
		setup_cell_arrays();
		setup_cell_neighbours();
		assign_all_to_cells();
	}

	if(nlist_flag && rebuild_nlist_flag) { setup_nlist(); }

}

void box__neighbour_list__cell_list__HS::c_packing_fraction(double packing_fraction_new) {

	double scale_factor;

	scale_factor = pow(r_packing_fraction()/packing_fraction_new,1.0/3.0);
	
	rescale_all(scale_factor);

}


void box__neighbour_list__cell_list__HS::generate_initial(){

	int i = 0;
	double part_en, rand_coord;


	for(int j=0; j<dim; j++) {
		rand_coord = rand_1->prob_gen()*box_size[j];
		particle_list[i]->c_coordinate_box(j, rand_coord);
		particle_list[i]->c_coordinate_true(j, rand_coord);
		particle_list[i]->c_coordinate_affine(j, rand_coord);
	}
	part_interaction_1->assign_type(particle_list[i]);

	part_en = particle_energy(i);
	while(part_en != part_en) {
		for(int j=0; j<dim; j++) {
			rand_coord = rand_1->prob_gen()*box_size[j];
			particle_list[i]->c_coordinate_box(j, rand_coord);
			particle_list[i]->c_coordinate_true(j, rand_coord);
			particle_list[i]->c_coordinate_affine(j, rand_coord);
		}
		part_en = particle_energy(i);
	}

	i++;
	

	while(i<number_particles) {
	
		for(int j=0; j<dim; j++) {
			rand_coord = rand_1->prob_gen()*box_size[j];
			particle_list[i]->c_coordinate_box(j, rand_coord);
			particle_list[i]->c_coordinate_true(j, rand_coord);
			particle_list[i]->c_coordinate_affine(j, rand_coord);
		}
		part_interaction_1->assign_type(particle_list[i]);

		part_en = particle_energy(i);
		while(part_en != part_en) {
			for(int j=0; j<dim; j++) {
				rand_coord = rand_1->prob_gen()*box_size[j];
				particle_list[i]->c_coordinate_box(j, rand_coord);
				particle_list[i]->c_coordinate_true(j, rand_coord);
				particle_list[i]->c_coordinate_affine(j, rand_coord);
			}
			part_en = particle_energy(i);
		}

		i++;
	
	}

	if(cells_flag) assign_all_to_cells();
	if(nlist_flag) {
		setup_nlist_arrays();
		setup_nlist();
	}

	find_smallest_separation();
}

void box__neighbour_list__cell_list__HS::generate_initial(FILE *in_file) {

	read_configuration(in_file);

	if(cells_flag) assign_all_to_cells();
	if(nlist_flag) {
		setup_nlist_arrays();
		setup_nlist();
	}

	find_smallest_separation();

}
