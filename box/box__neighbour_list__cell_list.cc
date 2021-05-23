#include "box__neighbour_list__cell_list.h"

int box__neighbour_list__cell_list::r_dim() { return dim; }

double box__neighbour_list__cell_list::r_box_size(int i_dim) { return box_size[i_dim]; }

int box__neighbour_list__cell_list::r_pb_flag(int i_dim) { return pb_flag[i_dim]; } 

double box__neighbour_list__cell_list::r_temperature() { return temperature; }

void box__neighbour_list__cell_list::c_temperature(double new_temperature) { temperature = new_temperature; }

double box__neighbour_list__cell_list::r_pressure() { return pressure; }

void box__neighbour_list__cell_list::c_pressure(double new_pressure) { pressure = new_pressure; }

double box__neighbour_list__cell_list::r_volume() {

	double volume = 1.0;

	for(int j = 0; j < dim; j++) volume = volume*box_size[j];

	return volume;

}

double box__neighbour_list__cell_list::r_number_density() {

	return number_particles/r_volume();

}

double box__neighbour_list__cell_list::r_packing_fraction() { return PI/6.0*ave_sigma3*r_number_density(); }

double box__neighbour_list__cell_list::r_ave_sigma_3() { return ave_sigma3; }

void box__neighbour_list__cell_list::c_box_size(vector <double> new_box_size) { box_size = new_box_size; }

void box__neighbour_list__cell_list::c_packing_fraction(double packing_fraction_new) {

	double scale_factor;

	scale_factor = pow(r_packing_fraction()/packing_fraction_new,1.0/3.0);
	
	rescale_all(scale_factor);

	reset_displacements();

}

double box__neighbour_list__cell_list::r_timestamp() { return timestamp; }

void box__neighbour_list__cell_list::c_timestamp(double new_timestamp) { timestamp = new_timestamp; }

int box__neighbour_list__cell_list::r_number_particles() { return number_particles; }

particle* box__neighbour_list__cell_list::r_particle(int part_i) { return particle_list[part_i]; }

particle_interaction* box__neighbour_list__cell_list::r_particle_interaction() { return part_interaction_1; }

vector <particle*> box__neighbour_list__cell_list::r_p_vector() { return particle_list; }

void box__neighbour_list__cell_list::copy_coordinates(vector <particle*> coords_in) {

	if(cells_flag == 1) {
		setup_cell_arrays();
		setup_cell_neighbours();
	}


	for(int i = 0; i<number_particles; i++) {
		particle_list[i]->c_type(coords_in[i]->r_type());
		particle_list[i]->c_size(coords_in[i]->r_size());
		for(int j = 0; j<dim; j++) {
			particle_list[i] -> c_coordinate_box(j, coords_in[i]->r_coordinate_box(j));
			particle_list[i] -> c_coordinate_true(j, coords_in[i]->r_coordinate_true(j));
			particle_list[i] -> c_coordinate_affine(j, coords_in[i]->r_coordinate_affine(j));
		}
	}
		
	if(cells_flag) assign_all_to_cells();
	if(nlist_flag) {
		setup_nlist_arrays();
		setup_nlist();
	}
}

void box__neighbour_list__cell_list::generate_initial(){

	int i = 0;
	double part_en, rand_coord_j;


	for(int j=0; j<dim; j++) {
		rand_coord_j = rand_1->prob_gen()*box_size[j];
		particle_list[i]->c_coordinate_box(j, rand_coord_j);
		particle_list[i]->c_coordinate_true(j, rand_coord_j);
		particle_list[i]->c_coordinate_affine(j, rand_coord_j);
	}

	part_interaction_1->assign_type(particle_list[i]);

	part_en = particle_energy(i);
	while(part_en != part_en) {
		for(int j=0; j<dim; j++) {
			rand_coord_j = rand_1->prob_gen()*box_size[j];
			particle_list[i]->c_coordinate_box(j, rand_coord_j);
			particle_list[i]->c_coordinate_true(j, rand_coord_j);
			particle_list[i]->c_coordinate_affine(j, rand_coord_j);
		}
		part_en = particle_energy(i);
	}

	i++;
	

	while(i<number_particles) {
	
		for(int j=0; j<dim; j++) {
			rand_coord_j = rand_1->prob_gen()*box_size[j];
			particle_list[i]->c_coordinate_box(j, rand_coord_j);
			particle_list[i]->c_coordinate_true(j, rand_coord_j);
			particle_list[i]->c_coordinate_affine(j, rand_coord_j);
		}

		part_interaction_1->assign_type(particle_list[i]);
	
		part_en = particle_energy(i);
		while(part_en != part_en) {
			for(int j=0; j<dim; j++) {
				rand_coord_j = rand_1->prob_gen()*box_size[j];
				particle_list[i]->c_coordinate_box(j, rand_coord_j);
				particle_list[i]->c_coordinate_true(j, rand_coord_j);
				particle_list[i]->c_coordinate_affine(j, rand_coord_j);
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
}

void box__neighbour_list__cell_list::generate_initial(FILE *in_file) {

	read_configuration(in_file);

	if(cells_flag) {
		setup_cell_arrays();
		setup_cell_neighbours();
		assign_all_to_cells();
	}
	if(nlist_flag) {
		setup_nlist_arrays();
		setup_nlist();
	}

}

double box__neighbour_list__cell_list::metric_true(int part_1, int part_2){

	double r, delta, coord_true_1, coord_true_2;

	r = 0.0;
	
	for(int j = 0; j<dim; j++) 
	{
		coord_true_1 = particle_list[part_1]->r_coordinate_true(j);
		coord_true_2 = particle_list[part_2]->r_coordinate_true(j);

		delta = coord_true_1 - coord_true_2;
		r = r + delta*delta;
	}

	return r;
}

double box__neighbour_list__cell_list::metric(int part_1, int part_2){
	
	double r, delta;
    
	r=0.0;
    
	for(int j = 0; j < dim; j++)
	{
		delta = particle_list[part_1]->r_coordinate_box(j)-particle_list[part_2]->r_coordinate_box(j);
		if(delta > box_size[j]*0.5 && pb_flag[j] == 1) delta = delta - box_size[j];
		if(delta <= -box_size[j]*0.5 && pb_flag[j] == 1) delta = delta + box_size[j];
		r = r + delta*delta;
	}
    
	return r;
}

void box__neighbour_list__cell_list::displace_particle(int part_1, vector <double> displacement) {

	for(int j = 0; j < displacement.size(); j++) {

		particle_list[part_1]->c_coordinate_box(j,particle_list[part_1]->r_coordinate_box(j)+displacement[j]);
		particle_list[part_1]->c_coordinate_true(j,particle_list[part_1]->r_coordinate_true(j)+displacement[j]);
		particle_list[part_1]->c_coordinate_affine(j,particle_list[part_1]->r_coordinate_affine(j)+displacement[j]);

		if(particle_list[part_1]->r_coordinate_box(j) < 0.0 && pb_flag[j] == 1){
			particle_list[part_1]->c_coordinate_box(j,particle_list[part_1]->r_coordinate_box(j) + box_size[j]);
		}
            
		if(particle_list[part_1]->r_coordinate_box(j) >= box_size[j] && pb_flag[j] == 1){
			particle_list[part_1]->c_coordinate_box(j, particle_list[part_1]->r_coordinate_box(j) - box_size[j]);
		}

		if(nlist_flag) disp_since_rebuild[part_1][j] = disp_since_rebuild[part_1][j] + displacement[j];
	}


	if(cells_flag) check_and_update_cell(part_1);

	if(nlist_flag) check_particle_nlist_and_update(part_1);

}

void box__neighbour_list__cell_list::check_particle_nlist_and_update(int part_1) {

	double dr2;

	if(nlist_flag) {

		dr2 = 0.0;
		for(int j = 0; j<dim; j++) dr2 = dr2 + disp_since_rebuild[part_1][j]*disp_since_rebuild[part_1][j];

		nlist_rebuild_limit = SKIN/2.0;
		nlist_rebuild_limit = nlist_rebuild_limit*nlist_rebuild_limit;

		if(dr2 > nlist_rebuild_limit) setup_nlist();
	}
}

void box__neighbour_list__cell_list::rescale_all(double scale_factor) {
	
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
			nlist_rebuild_limit = (SKIN + part_interaction_1->interaction_range(particle_list[i]->r_size())/2.0 - part_interaction_1->interaction_range_max()/2.0)/2.0;
			nlist_rebuild_limit = nlist_rebuild_limit*nlist_rebuild_limit;
			if(dr > nlist_rebuild_limit) rebuild_nlist_flag = 1;
		}
	}

	for(int j = 0; j < dim; j++) box_size[j] = box_size[j]*scale_factor;

	packing_fraction = r_packing_fraction();

	if(cells_flag) {
		setup_cell_arrays();
		setup_cell_neighbours();
		assign_all_to_cells();
	}

	if(nlist_flag && rebuild_nlist_flag) { setup_nlist(); }

}

double box__neighbour_list__cell_list::total_energy_calc(){

	double energy_acc = 0.0;

	for(int part_1 = 0; part_1<number_particles; ++part_1) energy_acc = energy_acc + particle_energy(part_1);

	return energy_acc;
}

double box__neighbour_list__cell_list::particle_energy(int part_1){
	
	double total_energy=0.0;
	int cell, neigh_cell, part_2, boundary_on, field_on;

	boundary_on = 0; field_on = 0;
	

	if(nlist_flag && nlist_done) { //Do neighbour list procedure
		for(int j = 0; j < nlist[part_1].size(); j++) {
			total_energy = total_energy + interaction_energy(part_1,nlist[part_1][j]);
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
				if(interaction_energy(part_1,part_2) != interaction_energy(part_1,part_2)) printf("%d %d %f %f\n", part_1, part_2, interaction_energy(part_1, part_2), total_energy);
			}
		}
	}

	if(boundary_on) total_energy = total_energy + boundary_1->boundary_int_energy(particle_list[part_1]->r_coordinate_box_vec(), particle_list[part_1]->r_size());
	if(field_on) total_energy = total_energy + ext_field_1->field_energy(particle_list[part_1]->r_coordinate_box_vec(), particle_list[part_1]->r_size());
    
	return total_energy;

}

double box__neighbour_list__cell_list::interaction_energy(int part_1, int part_2) {

	double energy, r_ij_2, s_i, s_j;

	r_ij_2 = metric(part_1, part_2);
	s_i = particle_list[part_1]->r_size();
	s_j = particle_list[part_2]->r_size();

	energy = part_interaction_1->interaction_energy(r_ij_2, s_i, s_j);

	return energy;

}

int box__neighbour_list__cell_list::bond(int part_1, int part_2) {

	double bond;

	bond = part_interaction_1->bonded(metric(part_1, part_2), particle_list[part_1]->r_size(), particle_list[part_2]->r_size());

	return bond;

}

vector <int> box__neighbour_list__cell_list::r_particle_bond_list(int part_1) {
	
	int cell, neigh_cell, part_2;
	vector <int> bond_list_i;

	if(nlist_flag && nlist_done) { //Do neighbour list procedure
		for(int j = 0; j < nlist[part_1].size(); j++) {
			if(bond(part_1, part_2)) bond_list_i.push_back(nlist[part_1][j]);
		}
	}
	else if(cells_flag && cells_done) { //Do calc with cells

		cell = part_cell[part_1][0];

		for(int k = 0; k < dim*dim*dim; k++) {
			neigh_cell = cell_neighbours[cell][k];
			for(int j = 0; j<cell_contents[neigh_cell].size(); j++) {
				part_2 = cell_contents[neigh_cell][j];
				if(part_2 != part_1) { 
					if(bond(part_1, part_2)) bond_list_i.push_back(part_2);
				}
			}
		}
	}
	else { //Do normal procedure...
		for(int part_2=0; part_2<number_particles; part_2++){
			if( part_2 != part_1) {
				if(bond(part_1, part_2)) bond_list_i.push_back(part_2);
			}
		}
	}
    
	return bond_list_i;

}

vector <int> box__neighbour_list__cell_list::r_particle_closer_than_list(int part_1, double range) {
	
	int cell, neigh_cell, part_2;
	vector <int> bond_list_i;
	double range_2;

	range_2 = range*range;

	if(nlist_flag && nlist_done) { //Do neighbour list procedure
		for(int j = 0; j < nlist[part_1].size(); j++) {
			if(metric(part_1, nlist[part_1][j]) < range_2) bond_list_i.push_back(nlist[part_1][j]);
		}
	}
	else if(cells_flag && cells_done) { //Do calc with cells

		cell = part_cell[part_1][0];

		for(int k = 0; k < dim*dim*dim; k++) {
			neigh_cell = cell_neighbours[cell][k];
			for(int j = 0; j<cell_contents[neigh_cell].size(); j++) {
				part_2 = cell_contents[neigh_cell][j];
				if(part_2 != part_1) { 
					if(metric(part_1, part_2) < range_2) bond_list_i.push_back(part_2);
				}
			}
		}
	}
	else { //Do normal procedure...
		for(int part_2=0; part_2<number_particles; part_2++){
			if( part_2 != part_1) {
				if(metric(part_1, part_2) < range_2) bond_list_i.push_back(part_2);
			}
		}
	}
    
	return bond_list_i;

}

void box__neighbour_list__cell_list::setup_cell_arrays() {

	number_cells.resize(dim);
	cell_size.resize(dim);

	cell_size_aim = 1.1 + MC_dr;
	if(nlist_flag) cell_size_aim = SKIN + part_interaction_1->interaction_range_max();

	total_cells = 1;
	
	for(int j =0; j<dim; j++) {
		number_cells[j] = box_size[j]/cell_size_aim;
		cell_size[j] = box_size[j]/number_cells[j];
		total_cells = total_cells*number_cells[j];
	}

	cell_neighbours.resize(total_cells);
	cell_contents.resize(total_cells);
	for(int i = 0; i<total_cells; i++) {
		cell_neighbours[i].resize(dim*dim*dim,0);
	}
	part_cell.resize(number_particles);
	for(int i = 0; i<number_particles; i++) part_cell[i].resize(2,0);
}

void box__neighbour_list__cell_list::setup_cell_neighbours() {

	int cell_no, neigh_no;
	vector<int> mod_0;
	vector<int> mod_1;
	vector<int> mod_2;

	mod_0.resize(dim,0);
	mod_1.resize(dim,0);
	mod_2.resize(dim,0);

	for(int c_0 = 0; c_0<number_cells[0]; c_0++) {
		for(int c_1 = 0; c_1<number_cells[1]; c_1++) {
			for(int c_2 = 0; c_2<number_cells[2]; c_2++) {

				mod_0[0] = -1; mod_0[1] = 0; mod_0[2] = 1; 
				mod_1[0] = -1; mod_1[1] = 0; mod_1[2] = 1; 
				mod_2[0] = -1; mod_2[1] = 0; mod_2[2] = 1; 

				if(c_0 == 0) mod_0[0] = (number_cells[0] - 1.0);
				if(c_0 == number_cells[0]-1) mod_0[2] = -(number_cells[0] - 1.0);
				if(c_1 == 0) mod_1[0] = (number_cells[1] - 1.0);
				if(c_1 == number_cells[1]-1) mod_1[2] = -(number_cells[1] - 1.0);
				if(c_2 == 0) mod_2[0] = (number_cells[2] - 1.0);
				if(c_2 == number_cells[2] - 1) mod_2[2] = -(number_cells[2] - 1.0);

				cell_no = c_2 + number_cells[2]*(c_1 + number_cells[1]*c_0);
				
				for(int l = 0; l<3; l++) {
					for(int m = 0; m<3; m++) {
						for(int n = 0; n<3; n++) {

							neigh_no = l + dim*(m + dim*n);
							cell_neighbours[cell_no][neigh_no] = cell_no + mod_2[l] + mod_1[m]*number_cells[2] + mod_0[n]*number_cells[2]*number_cells[1];
						}
					}
				}
			}
		}
	}

}

int box__neighbour_list__cell_list::cell_number(int part){

	int cell_no;
	vector<int> cell_coords;

	cell_coords.resize(dim);

	for(int j = 0; j<dim; j++) cell_coords[j] = floor(particle_list[part]->r_coordinate_box(j)/cell_size[j]);

	cell_no = cell_coords[2] + number_cells[2]*(cell_coords[1] + number_cells[1]*cell_coords[0]);

	return cell_no;
}

void box__neighbour_list__cell_list::assign_to_cell(int part){

	int cell_no, cell_occup;
	cell_no = cell_number(part);

	cell_contents[cell_no].push_back(part);

	part_cell[part][0] = cell_no;
	part_cell[part][1] = (cell_contents[cell_no].size()-1);
}

void box__neighbour_list__cell_list::remove_from_cell(int part) {

	int cell_no = part_cell[part][0];
	int pos = part_cell[part][1];
	int occup = (cell_contents[cell_no].size()-1);


	if(occup == pos) {
		cell_contents[cell_no].pop_back();
	}
	else {
		part_cell[cell_contents[cell_no][occup]][1] = pos;
		cell_contents[cell_no][pos] = cell_contents[cell_no][occup];
		cell_contents[cell_no].pop_back();
	}
}

void box__neighbour_list__cell_list::check_and_update_cell(int part) {

	int cell_no_stored = part_cell[part][0];
	int cell_no_current = cell_number(part);

	if(cell_no_stored != cell_no_current) {

		remove_from_cell(part);
		assign_to_cell(part);

	}
}

void box__neighbour_list__cell_list::assign_all_to_cells() {

	for(int i = 0; i<number_particles; i++) {
		assign_to_cell(i);
	}

	cells_done = 1;
}

void box__neighbour_list__cell_list::setup_nlist_arrays() {
	
	nlist_cut_2 = (SKIN + part_interaction_1->interaction_range_max());
	nlist_cut_2 = nlist_cut_2*nlist_cut_2;

	nlist.resize(number_particles);

	disp_since_rebuild.resize(number_particles);
	for(int i = 0; i<number_particles; i++) disp_since_rebuild[i].assign(dim,0.0);

}

void box__neighbour_list__cell_list::setup_nlist() {

	int cell_part_2, cell, neigh_cell, nlist_cut_local_2_i, nlist_cut_local_2_j;
	double r_2, nlist_cut_local;

	for(int i = 0; i<number_particles; i++){
		disp_since_rebuild[i].clear();
		disp_since_rebuild[i].resize(dim, 0.0);
		nlist[i].resize(0);
	}

	if(cells_flag && cells_done) {
		for(int part_1 = 0; part_1<number_particles; part_1++) {

			cell = part_cell[part_1][0];

			for(int k = 0; k < dim*dim*dim; k++) {
				neigh_cell = cell_neighbours[cell][k];
				for(int j = 0; j<cell_contents[neigh_cell].size(); j++) {
					cell_part_2 = cell_contents[neigh_cell][j];
					if(cell_part_2 != part_1) { 
						r_2 = metric(part_1, cell_part_2);
						nlist_cut_local = SKIN + part_interaction_1->interaction_width(particle_list[part_1]->r_size(), particle_list[cell_part_2]->r_size());
						if(r_2 <= nlist_cut_local*nlist_cut_local) {
							nlist[part_1].push_back(cell_part_2);
						}
					}
				}

			}
		}
	}
	else {
		for(int part_1 = 0; part_1 < number_particles; part_1++) {

			for(int part_2 = part_1 + 1; part_2 < number_particles; part_2++) {

				nlist_cut_local = SKIN + part_interaction_1->interaction_width(particle_list[part_1]->r_size(), particle_list[part_2]->r_size());

				r_2 = metric(part_1, part_2);

				if(r_2 < nlist_cut_local*nlist_cut_local) {
					nlist[part_1].push_back(part_2);
					nlist[part_2].push_back(part_1);
				}
			}
		}
	}
	nlist_done = 1;
}

void box__neighbour_list__cell_list::dump_configuration(char out_file_name[], int seed, char traj_start[], double i_dr, double i_dv, char type_stamp[]) {

	FILE *traj_out;

	traj_out = fopen(out_file_name,"a");

	fprintf(traj_out, "%d\n", number_particles);

	fprintf(traj_out, "#t %f ", r_timestamp());
	fprintf(traj_out, "#BOX ");
	for(int j = 0; j < dim; j++) fprintf(traj_out, "%f ", box_size[j]);
	fprintf(traj_out, "#T %f ", r_temperature());
	fprintf(traj_out, "#P %f ", r_pressure());
	fprintf(traj_out, "#V %f ", r_volume());
	fprintf(traj_out, "#RHO %f ", r_number_density());
	fprintf(traj_out, "#PHI %f ", r_packing_fraction());
	fprintf(traj_out, "#RNG_SEED %d ", seed);
	fprintf(traj_out, "#TRAJ_START %s ", traj_start);
	fprintf(traj_out, "#dr %f ", i_dr);
	fprintf(traj_out, "#dv %f ", i_dv);
	fprintf(traj_out, "#type_stamp %s ", type_stamp);
	fprintf(traj_out, "\n");

	for(int i = 0; i<number_particles; ++i) {

		fprintf(traj_out, "%.12lf ", particle_list[i]->r_size());

		for(int j = 0; j<dim; j++) fprintf(traj_out, "%.12lf ", particle_list[i]->r_coordinate_box(j));

		for(int j = 0; j<dim; j++) fprintf(traj_out, "%.12lf ", particle_list[i]->r_coordinate_true(j));

		for(int j = 0; j<dim; j++) fprintf(traj_out, "%.12lf ", particle_list[i]->r_coordinate_affine(j));

		fprintf(traj_out, "\n");

	}

	fclose(traj_out);
}

void box__neighbour_list__cell_list::dump_restart(char out_file_name[], int jobid, char traj_start[], double i_dr, double i_dv, char type_stamp[]) {
 
	FILE *traj_out;
	
	traj_out = fopen(out_file_name,"w");

	fprintf(traj_out, "%d\n", number_particles);
	
	fprintf(traj_out, "#t %f ", r_timestamp());
	fprintf(traj_out, "#BOX ");
	for(int j = 0; j < dim; j++) fprintf(traj_out, "%f ", box_size[j]);
	fprintf(traj_out, "#T %f ", r_temperature());
	fprintf(traj_out, "#P %f ", r_pressure());
	fprintf(traj_out, "#JOBID %d ", jobid);
	fprintf(traj_out, "#TRAJ_START %s ", traj_start);
	fprintf(traj_out, "#dr %f ", i_dr);
	fprintf(traj_out, "#dv %f ", i_dv);
	fprintf(traj_out, "#type_stamp %s ", type_stamp);
	fprintf(traj_out, "\n");
 
	for(int i = 0; i<number_particles; ++i) {
 
		fprintf(traj_out, "%.12lf ", particle_list[i]->r_size());
 
		for(int j = 0; j<dim; j++) fprintf(traj_out, "%.12lf ", particle_list[i]->r_coordinate_box(j));
 
		for(int j = 0; j<dim; j++) fprintf(traj_out, "%.12lf ", particle_list[i]->r_coordinate_true(j));
 
		for(int j = 0; j<dim; j++) fprintf(traj_out, "%.12lf ", particle_list[i]->r_coordinate_affine(j));
 
		fprintf(traj_out, "\n");
 
	}
 
	fclose(traj_out);
}


void box__neighbour_list__cell_list::ignore_configuration(FILE *in_file) {

	char line[1024];

	fgets(line,1000,in_file);
	fgets(line,1000,in_file);

	for(int i = 0; i<number_particles; i++) fgets(line,1000,in_file);
}

void box__neighbour_list__cell_list::read_configuration(FILE *in_file) {

	char dump[1024];
	double dmp, coord_box_in, coord_true_in, coord_affine_in, p_size, sigma3_acc = 0.0, largest_sigma;

	largest_sigma = 0.0;

	fgets(dump,1000,in_file);

	fgets(dump,1000,in_file);

	sscanf(dump, "%*s %lf %*s %lf %lf %lf %*s %lf %*s %lf %*s %*lf %*s %*lf %*s %*lf %*s %*lf", &timestamp, &box_size[0], &box_size[1], &box_size[2], &temperature, &pressure);

	for(int i=0; i<number_particles; i++) {

		fscanf(in_file,"%lf ",&p_size);
		particle_list[i]->c_type(1);
		particle_list[i]->c_size(p_size);
		if(p_size > largest_sigma) largest_sigma = p_size;
		sigma3_acc = sigma3_acc + p_size*p_size*p_size;

		for(int j = 0; j < dim; j++) {
			fscanf(in_file, "%lf ", &coord_box_in);
			if(coord_box_in >= box_size[j] && pb_flag[j]) coord_box_in = coord_box_in - box_size[j];
			if(coord_box_in < 0.0 && pb_flag[j]) coord_box_in = coord_box_in + box_size[j];
			particle_list[i]->c_coordinate_box(j, coord_box_in);
		}

		for(int j = 0; j < dim; j++) {
			fscanf(in_file, "%lf ", &coord_true_in);
			particle_list[i]->c_coordinate_true(j, coord_true_in);
		}

		for(int j = 0; j < dim; j++) {
			fscanf(in_file, "%lf ", &coord_affine_in);
			particle_list[i]->c_coordinate_affine(j, coord_affine_in);
		}
	
		fscanf(in_file, "\n");
	}

	ave_sigma3 = sigma3_acc/number_particles;
	packing_fraction = PI/6.0*r_number_density()*ave_sigma3;
	part_interaction_1->c_interaction_range(largest_sigma);
}

void box__neighbour_list__cell_list::calc_centre_of_mass(int version_flag) {

	centre_of_mass.resize(dim, 0.0);

	for(int i = 0; i < number_particles; i++) {
		for(int j = 0; j < dim; j++) {
			if(version_flag == 1) centre_of_mass[j] += particle_list[i] -> r_coordinate_true(j);
			if(version_flag == 2) centre_of_mass[j] += particle_list[i] -> r_coordinate_box(j);
			if(version_flag == 3) centre_of_mass[j] += particle_list[i] -> r_coordinate_affine(j);
		}
	}

	for(int j = 0; j < dim; j++) centre_of_mass[j] = centre_of_mass[j]/number_particles;

}

void box__neighbour_list__cell_list::reset_displacements() {

	for(int i = 0; i < number_particles; i++) {
		for(int j = 0; j < dim; j++) {
			particle_list[i]->c_coordinate_true(j, particle_list[i]->r_coordinate_box(j));
			particle_list[i]->c_coordinate_affine(j, particle_list[i]->r_coordinate_box(j));
		}
	}

}

double box__neighbour_list__cell_list::r_centre_of_mass(int i_dim) { return centre_of_mass[i_dim]; }
