#include "MC_NPT_swap_poly_sticky_spheres_test.h"

int main(int argc, char *argv[]){

	vector <int> pb_flag;
	vector <double> box_size;
	vector <double> time_series;
	vector <int> step;
	vector <vector <double> > time_out;
	int seed, next_step, jobid;
	double L, phi_i, pressure_f, time_step, log_step, time_end, count, waiting_time, simulation_time, time_current, dv, dr, lambda, u;
	FILE *config_in;
	FILE *sys_info_out;
	char file_name[1024] = {0};
	char sys_info_file_name[1024] = {0};
	vector <std::string> file_names;

	stringstream convert_pressure_f(argv[2]); convert_pressure_f >> pressure_f;
	stringstream convert_u(argv[3]); convert_u >> u;
	stringstream convert_lambda(argv[4]); convert_lambda >> lambda;
	stringstream convert_dr(argv[5]); convert_dr >> dr;
	stringstream convert_dv(argv[6]); convert_dv >> dv;

	stringstream convert_jobid(argv[8]); convert_jobid >> jobid;

	JRand_dyn *rand_1 = new JRand_dyn();
	seed = (unsigned int)(jobid+time(NULL))%1000;
	rand_1->srandom(seed);

	count = 0;
	waiting_time = 0.0;
	while(waiting_time <= WAITING_TIME_MAX) {

		time_series.clear();
		time_series.push_back(waiting_time);
		time_out.push_back(time_series);
		count += 1;
		waiting_time = pow(10.0, count);

	}
	
	step.resize(time_out.size(), 0);

	count = -5;
	simulation_time = 0.0;
	while(simulation_time <= SIMULATION_TIME_MAX) {

		log_step = exp(count);
		for(int i = 0; i < time_out.size(); i++) time_out[i].push_back(time_out[i][time_out[i].size()-1] + log_step);
		simulation_time += log_step;
		count += 0.66;

	}


	for(int i = 0; i < time_out.size(); i++) {

		stringstream convert_w_time; convert_w_time << time_out[i][0]; const std::string str_w_time = convert_w_time.str();

		strcpy(file_name,getenv("PWD"));
		strcat(file_name,"/trajectories");
		mkdir(file_name, S_IRWXU);
		strcat(file_name,"/");
		strcat(file_name,argv[0]);
		strcat(file_name,"__");
		strcat(file_name, argv[8]);
		strcat(file_name, "__");
		strcat(file_name, str_w_time.c_str());
		strcat(file_name, ".xyz");

		file_names.push_back(file_name);

	}

	strcpy(sys_info_file_name,getenv("PWD"));
	strcat(sys_info_file_name,"/sys_info_data");
	mkdir(sys_info_file_name, S_IRWXU);
	strcat(sys_info_file_name, "/");
	strcat(sys_info_file_name, argv[0]);
	strcat(sys_info_file_name, "__sys_info__");
	strcat(sys_info_file_name, argv[8]);
	strcat(sys_info_file_name, ".txt");

	pb_flag.resize(DIM); box_size.resize(DIM);
	pb_flag[0] = PB_X; pb_flag[1] = PB_Y; pb_flag[2] = PB_Z;
	box_size[0] = BOX_SIZE_X; box_size[1] = BOX_SIZE_Y; box_size[2] = BOX_SIZE_Z;

	part_poly_sphere *particle_1 = new part_poly_sphere(DIM);
	particle_interaction_sticky_spheres *int_poly_sticky_1 = new particle_interaction_sticky_spheres(lambda, u, 0.0, rand_1);
	no_field *field_1 = new no_field();
	no_boundary *boundary_1 = new no_boundary();

	config_in = fopen(argv[1], "r");

	box__neighbour_list__cell_list *box_ti = new box__neighbour_list__cell_list(NUMBER_PARTICLES, TEMPERATURE, DR, box_size, pb_flag, particle_1, int_poly_sticky_1, rand_1, field_1, boundary_1);
	box_ti -> generate_initial(config_in);

	phi_i = box_ti->r_packing_fraction();

	fclose(config_in);

	integrator_MC_NPT_swap *integrator_1 = new integrator_MC_NPT_swap(box_ti, DR, dv, rand_1);

	time_current = 0.0;

	sys_info_out = fopen(sys_info_file_name, "a");
	fprintf(sys_info_out, "#RNG_Seed %d #phi_i %f #pressure_f %f #u %f #lambda %f #dr %f #Init_config %s #type_stamp %s\n", seed, phi_i, pressure_f, u, lambda, dr, argv[1], argv[7]);
	fprintf(sys_info_out, "#KEY #time #packing_fraction #pressure #rho #total_energy\n\n");
	fclose(sys_info_out);

	box_ti->c_timestamp(0.0);
	box_ti -> c_pressure(pressure_f);

	while(time_current < SIMULATION_TIME_MAX) {

		time_step = SIMULATION_TIME_MAX - time_current;

		for(int i = 0; i < time_out.size(); i++) {
			if(time_out[i][step[i]] - time_current <= time_step) {
				time_step = time_out[i][step[i]] - time_current;
				next_step = i;
			}
		}

		step[next_step] += 1;

		integrator_1->run_dynamics_particle_volume_swap(time_step*double(NUMBER_PARTICLES));

		time_current = time_current + time_step;

		box_ti -> dump_configuration(strdup(file_names[next_step].c_str()), seed, argv[1], DR, dv, argv[7]);

		sys_info_out = fopen(sys_info_file_name, "a");
		fprintf(sys_info_out, "%f %f %f %f %f\n", box_ti->r_timestamp(), box_ti->r_packing_fraction(), box_ti->r_pressure(), box_ti->r_number_density(), box_ti->total_energy_calc());
		fclose(sys_info_out);



	}

	return 0;
}
