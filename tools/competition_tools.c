#include "competition_tools.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "orbit_calculator/itin_tool.h"
#include "orbit_calculator/transfer_calc.h"
#include "file_io.h"

enum FILE_TYPE {COMP_FILE_PLANET, COMP_FILE_COMET, COMP_FILE_ASTEROID};

void parse_celestial_body_line(const char *line, enum FILE_TYPE type, Body *new_body) {
	// Ignore comment or empty lines
	if(line[0] == '#' || strlen(line) < 3)
		return;
	
	if(type == COMP_FILE_PLANET) {
		// Use sscanf to parse CSV values
		sscanf(line,
			   "%d,%63[^,],%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
			   &new_body->id,
			   new_body->name,
			   &new_body->mu,
			   &new_body->radius,
			   &new_body->orbit.a,
			   &new_body->orbit.e,
			   &new_body->orbit.i,
			   &new_body->orbit.raan,
			   &new_body->orbit.arg_peri,
			   &new_body->orbit.ta,
			   &new_body->scale_height
		);
		new_body->color[0] = 1;
		new_body->color[1] = 0;
		new_body->color[2] = 1;
	} else {
		// Use sscanf to parse CSV values
		sscanf(line,
			   "%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
			   &new_body->id,
			   &new_body->orbit.a,
			   &new_body->orbit.e,
			   &new_body->orbit.i,
			   &new_body->orbit.raan,
			   &new_body->orbit.arg_peri,
			   &new_body->orbit.ta,
			   &new_body->scale_height
		);
		if(type == COMP_FILE_COMET) {
			sprintf(new_body->name, "C%d", new_body->id);
			new_body->color[0] = 0;
			new_body->color[1] = 0.3;
			new_body->color[2] = 1;
		} else {
			sprintf(new_body->name, "A%d", new_body->id);
			new_body->color[0] = 0.5;
			new_body->color[1] = 0.5;
			new_body->color[2] = 0.5;
		}
	}
	
	new_body->mu *= 1e9;
	new_body->radius *= 1e3;
	new_body->orbit.a *= 1e3;
	new_body->orbit.i = deg2rad(new_body->orbit.i);
	new_body->orbit.raan = deg2rad(new_body->orbit.raan);
	new_body->orbit.arg_peri = deg2rad(new_body->orbit.arg_peri);
	new_body->orbit.ta = calc_true_anomaly_from_mean_anomaly(new_body->orbit, deg2rad(new_body->orbit.ta));
}

void load_competition_file(CelestSystem *system, char *filename, enum FILE_TYPE type) {
	FILE *file = fopen(filename, "r");
	if(!file) {
		perror("Failed to open file");
		free(system);
		return;
	}
	char line[256];  // Buffer for each line
	// skip first line
	fgets(line, sizeof(line), file);
	
	while(fgets(line, sizeof(line), file)) {
		Body *body = new_body();
		parse_celestial_body_line(line, type, body);
		body->orbit.cb = system->cb;
		system->bodies[system->num_bodies] = body;
		system->num_bodies++;
	}
	fclose(file);
}

CelestSystem * load_competition_system(char *directory) {
	CelestSystem *system = new_system();
	system->num_bodies = 0;
	system->prop_method = ORB_ELEMENTS;
	sprintf(system->name, "Altaira System");
	system->ut0 = 0;
	
	system->cb = new_body();
	system->cb->id = 0;
	sprintf(system->cb->name, "Altaira");
	system->cb->color[0] = 1;
	system->cb->color[1] = 1;
	system->cb->color[2] = 0.3;
	system->cb->mu = 139348062043.343e9;
	system->cb->system = system;
	system->bodies = malloc(3000*sizeof(Body*));
	
	char filename[256];
	sprintf(filename, "%sgtoc13_planets.csv", directory);
	load_competition_file(system, filename, COMP_FILE_PLANET);
//	sprintf(filename, "%sgtoc13_comets.csv", directory);
//	load_competition_file(system, filename, COMP_FILE_COMET);
//	sprintf(filename, "%sgtoc13_asteroids.csv", directory);
//	load_competition_file(system, filename, COMP_FILE_ASTEROID);
	
	system->home_body = system->bodies[0];
	
	return system;
}

struct ItinStep * attach_initial_competition_state(struct ItinStep *step) {
	struct ItinStep *ptr = get_first(step);
	struct ItinStep *new_step = malloc(sizeof(struct ItinStep));
	new_step->body = NULL;
	new_step->next = malloc(sizeof(struct ItinStep *));
	new_step->next[0] = ptr;
	new_step->prev = NULL;
	new_step->had_low_perihelion = false;
	new_step->num_next_nodes = 1;
	ptr->prev = new_step;
	
	
	double phi, kappa;
	Vector3 v_inf_dep = subtract_vec3(ptr->next[0]->v_dep, ptr->v_body);
	Vector3 v_inf_arr;
	double x_error, y_error, z_error;
	Orbit orbit_arr;
	Orbit orbit_dep;
	OSV osv_dep;
	DataArray2 *x_error_data = data_array2_create();
	DataArray2 *y_error_data = data_array2_create();
	DataArray2 *z_error_data = data_array2_create();
	
	data_array2_insert_new(z_error_data, deg2rad(30), 1e20);
	data_array2_insert_new(z_error_data, deg2rad(-30), -1e20);
	do {
		phi = root_finder_monot_func_next_x(z_error_data);
		data_array2_clear(y_error_data);
		data_array2_insert_new(y_error_data, deg2rad(30), 1e20);
		data_array2_insert_new(y_error_data, deg2rad(-30), -1e20);
		do {
			kappa = root_finder_monot_func_next_x(y_error_data);
			printf("Phi: %f    Kappa: %f\n", rad2deg(phi), rad2deg(kappa));
			v_inf_arr = vec3(cos(phi),0,sin(phi));
			v_inf_arr = rotate_vector_around_axis(v_inf_arr, vec3(0,0,1), kappa);
			
			v_inf_arr = scale_vec3(v_inf_arr, mag_vec3(v_inf_dep));
			ptr->v_arr = add_vec3(v_inf_arr, ptr->v_body);
			
			orbit_arr = constr_orbit_from_osv(ptr->r, ptr->v_arr, ptr->body->orbit.cb);
			orbit_dep = orbit_arr;
			data_array2_clear(x_error_data);
			data_array2_insert_new(x_error_data, orbit_dep.ta, osv_from_orbit(orbit_dep).r.x + 200*AU);
			Orbit orbit_neg_inf = propagate_orbit_time(orbit_arr, -1e10);
			data_array2_insert_new(x_error_data, orbit_neg_inf.ta, osv_from_orbit(orbit_neg_inf).r.x + 200*AU);
			
			do {
				orbit_dep.ta = root_finder_monot_func_next_x(x_error_data);
				osv_dep = osv_from_orbit(orbit_dep);
				x_error = osv_dep.r.x - (-200*AU);
				data_array2_insert_new(x_error_data, orbit_dep.ta, osv_from_orbit(orbit_dep).r.x + 200*AU);
			} while(fabs(x_error) > 1);
			
			y_error = osv_dep.v.y;
			data_array2_insert_new(y_error_data, kappa, y_error);
		} while(fabs(y_error) > 1e-4);
		
		z_error = osv_dep.v.z;
		data_array2_insert_new(z_error_data, phi, z_error);
	} while(fabs(z_error) > 1e-4);
	
	data_array2_free(x_error_data);
	data_array2_free(y_error_data);
	data_array2_free(z_error_data);
	
	new_step->r = osv_from_orbit(orbit_dep).r;
	
	double dt = fabs(calc_orbit_time_since_periapsis(orbit_arr)-calc_orbit_time_since_periapsis(orbit_dep));
	new_step->date = ptr->date - dt/86400.0;
	
	ptr->v_dep = osv_from_orbit(orbit_dep).v;
	
	printf("dt: %f\n", ptr->date - new_step->date);
	print_vec3(v_inf_dep);
	print_vec3(v_inf_arr);
	print_vec3(new_step->r);
	print_vec3(subtract_vec3(new_step->r,vec3(-200*AU, 0, 0)));
	print_vec3(scale_vec3(new_step->r, 1/AU));
	print_vec3(ptr->v_dep);
	
	print_vec3(subtract_vec3(propagate_osv_time(osv_dep, ptr->body->orbit.cb, dt).r, ptr->r));
	print_vec3(subtract_vec3(propagate_osv_time(osv_from_orbit(orbit_arr), ptr->body->orbit.cb, -dt).r, osv_dep.r));
//	print_vec3(subtract_vec3(propagate_osv_time(osv_from_orbit(orbit_arr), ptr->body->orbit.cb, -dt).v, osv_dep.v));
//	print_vec3(propagate_osv_time(osv_from_orbit(orbit_arr), ptr->body->orbit.cb, -dt).v);
	print_vec3(subtract_vec3(osv_from_orbit(propagate_orbit_time(propagate_orbit_time(orbit_arr, -dt), dt)).r, ptr->r));
	print_vec3(subtract_vec3(osv_from_orbit(propagate_orbit_time(propagate_orbit_time(orbit_arr, -dt), dt)).r, ptr->r));
	OSV temp_osv = osv_from_orbit(orbit_arr);
	temp_osv.v = scale_vec3(temp_osv.v, -1);
	print_vec3(subtract_vec3(propagate_osv_time(temp_osv, ptr->body->orbit.cb, dt).r, osv_dep.r));

	OSV osv_arr = osv_from_orbit(orbit_arr);
	Lambert3 solution = calc_lambert3(ptr->r, osv_arr.r, dt, ptr->body->orbit.cb);
	
	OSV osv_arr2 = (OSV) {solution.r0, scale_vec3(solution.v0,-1)};
	OSV osv_dep2 = (OSV) {solution.r1, scale_vec3(solution.v1, -1)};
	
	print_vec3(subtract_vec3(osv_dep2.v, osv_dep.v));
	print_vec3(subtract_vec3(osv_arr2.v, ptr->v_arr));
	
//
//	print_orbit_info(orbit_dep);
//	print_orbit_info(orbit_arr);
	
	return new_step;
}

typedef struct Competition_Transfer {
	Body *body;
	double epoch;
	Vector3 r;
	double c3;
	double rp;
	struct Competition_Transfer *prev, *next;
} Competition_Transfer;

Competition_Transfer * get_first_comp_trans(Competition_Transfer *transfer) {
	if(transfer == NULL) return NULL;
	while(transfer->prev != NULL) transfer = transfer->prev;
	return transfer;
}

Competition_Transfer * get_last_comp_trans(Competition_Transfer *transfer) {
	if(transfer == NULL) return NULL;
	while(transfer->next != NULL) transfer = transfer->next;
	return transfer;
}

int get_num_comp_trans(Competition_Transfer *transfer) {
	if(transfer == NULL) return 0;
	transfer = get_first_comp_trans(transfer);
	int num_transfer = 1;
	while(transfer->next != NULL) {
		num_transfer++;
		transfer = transfer->next;
	}
	return num_transfer;
}

void free_competition_transfer_list(Competition_Transfer **transfers, int num_bodies) {
	for(int i = 0; i < num_bodies; i++) {
		Competition_Transfer *ptr = get_last_comp_trans(transfers[i]);
		if(ptr == NULL) continue;
		while(ptr->prev != NULL) {
			ptr = ptr->prev;
			free(ptr->next);
		}
		free(ptr);
	}
	free(transfers);
}

Vector3 calc_heliocentric_periapsis(Vector3 r_dep, Vector3 v_dep, Vector3 r_arr, Vector3 v_arr, CelestSystem *system) {
	Orbit orbit0 = constr_orbit_from_osv(r_dep, v_dep, system->cb);
	Orbit orbit1 = constr_orbit_from_osv(r_arr, v_arr, system->cb);
	
	Vector3 r;
	
	if(orbit1.ta > orbit0.ta && orbit1.ta < 2*M_PI-orbit0.ta) {
		r = r_dep;
	} else if(orbit1.ta > orbit0.ta) {
		r = r_arr;
	} else {
		orbit0.ta = 0;
		OSV osv = osv_from_orbit(orbit0);
		r = osv.r;
	}
	return r;
}

Competition_Transfer ** build_competition_transfer_from_itin(struct ItinStep *arr_step, CelestSystem *system) {
	struct ItinStep *itin_ptr = arr_step;
	
	Competition_Transfer **transfers = malloc(3000 *  sizeof(Competition_Transfer*));
	for(int i = 0; i < 3000; i++) transfers[i] = NULL;
	
	Competition_Transfer *transfer = malloc(sizeof(Competition_Transfer));
	transfer->prev = NULL;
	transfer->epoch = itin_ptr->date;
	transfer->body = itin_ptr->body;
	transfer->r = itin_ptr->r;
	double v_inf = mag_vec3(subtract_vec3(itin_ptr->v_arr, itin_ptr->v_body));
	transfer->c3 = v_inf*v_inf;
	transfer->rp = 0;
	transfer->next = transfers[transfer->body->id];
	transfers[transfer->body->id] = transfer;
	if(transfer->next != NULL) transfer->next->prev = transfer;
	
	while(itin_ptr != NULL) {
		transfer = malloc(sizeof(Competition_Transfer));
		transfer->prev = NULL;
		
		transfer->r = calc_heliocentric_periapsis(
				itin_ptr->prev->r, itin_ptr->v_dep,
				itin_ptr->r, itin_ptr->v_arr,
				system);
		transfer->rp = mag_vec3(transfer->r);
		
		transfer->epoch = (itin_ptr->date+itin_ptr->prev->date);
		transfer->body = itin_ptr->body->orbit.cb;
		Orbit orbit = constr_orbit_from_osv(itin_ptr->prev->r, itin_ptr->v_dep, itin_ptr->body->orbit.cb);
		transfer->c3 = -transfer->body->mu/orbit.a;
		transfer->next = transfers[transfer->body->id];
		transfers[transfer->body->id] = transfer;
		if(transfer->next != NULL) transfer->next->prev = transfer;
		
		if(itin_ptr->prev->prev == NULL) break;
		
		struct ItinStep *next_itin = itin_ptr;
		itin_ptr = itin_ptr->prev;
		
		transfer = malloc(sizeof(Competition_Transfer));
		transfer->prev = NULL;
		
		transfer->epoch = itin_ptr->date;
		transfer->body = itin_ptr->body;
		transfer->r = itin_ptr->r;
		v_inf = mag_vec3(subtract_vec3(itin_ptr->v_arr, itin_ptr->v_body));
		transfer->c3 = v_inf*v_inf;
		transfer->rp = get_flyby_periapsis(itin_ptr->v_arr, next_itin->v_dep, itin_ptr->v_body, itin_ptr->body);
		transfer->next = transfers[transfer->body->id];
		transfers[transfer->body->id] = transfer;
		if(transfer->next != NULL) transfer->next->prev = transfer;
	}
	
	struct ItinStep *next_itin = itin_ptr;
	itin_ptr = itin_ptr->prev;
	
	transfer = malloc(sizeof(Competition_Transfer));
	transfer->prev = NULL;
	
	transfer->epoch = itin_ptr->date;
	transfer->body = itin_ptr->body;
	transfer->r = itin_ptr->r;
	v_inf = mag_vec3(subtract_vec3(next_itin->v_dep, itin_ptr->v_body));
	transfer->c3 = v_inf*v_inf;
	transfer->rp = 0;
	transfer->next = transfers[transfer->body->id];
	transfers[transfer->body->id] = transfer;
	if(transfer->next != NULL) transfer->next->prev = transfer;
	
	return transfers;
}


bool is_grand_tour() {
	return false;
}


double calc_seasonal_penalty(Competition_Transfer *transfer) {
	double s = 0;
	Competition_Transfer *prev_ptr = transfer->prev;
	while(prev_ptr != NULL) {
		Vector3 r_i = norm_vec3(transfer->r);
		Vector3 r_j = norm_vec3(prev_ptr->r);
		double acosd = rad2deg(acos(dot_vec3(r_i, r_j)));
		s += exp(-(acosd*acosd)/50);
		prev_ptr = prev_ptr->prev;
	}
	return 0.1 + (0.9/(1+10*s));
}

double calc_flyby_velocity_penalty(Competition_Transfer *transfer) {
	double v_inf = sqrt(transfer->c3)/1e3;
	return 0.2 + exp(-v_inf/13) / (1+exp(-5*(v_inf-1.5)));;
}

double get_competition_flyby_score(Competition_Transfer *transfer) {
	if(transfer == 0) return 0;
	double w = transfer->body->scale_height;
	double s = calc_seasonal_penalty(transfer);
	double f = calc_flyby_velocity_penalty(transfer);
	return w*s*f;
}


double get_competition_body_score(Competition_Transfer *transfer) {
	if(transfer == 0) return 0;
	double score = 0;
	while(transfer != NULL) {
		score += get_competition_flyby_score(transfer);
		transfer = transfer->next;
	}
	return score;
}

double get_itin_competition_score(struct ItinStep *arr_step, CelestSystem *system) {
	double b = is_grand_tour() ? 1.2 : 1;
	
	Competition_Transfer **transfers = build_competition_transfer_from_itin(arr_step, system);
	
	double score = 0;
	Competition_Transfer *ptr;
	
	for(int i = 1; i < 3000; i++) {
		score += get_competition_body_score(transfers[i]);
	}

	ptr = transfers[0];

	bool was_below_0_05 = false;
	while(ptr != NULL) {
		if(ptr->rp/AU < 0.05) {
			if(was_below_0_05 || ptr->rp/AU < 0.01) { free_competition_transfer_list(transfers, 3000); return 0; }
			else was_below_0_05 = true;
		}
		ptr = ptr->next;
	}
	
	for(int i = 1; i < 3000; i++) {
		ptr = transfers[i];
		while(ptr != NULL) {
			if(ptr->rp > 0) {
				if(ptr->rp/ptr->body->radius - 1 < 0.1 || ptr->rp/ptr->body->radius - 1 > 100) {
					free_competition_transfer_list(transfers, 3000);
					return 0;
				}
			}
			ptr = ptr->next;
		}
	}
	
	free_competition_transfer_list(transfers, 3000);
	
	return score;
}

void print_itin_competition_score(struct ItinStep *arr_step, CelestSystem *system) {
	double b = is_grand_tour() ? 1.2 : 1;
	
	Competition_Transfer **transfers = build_competition_transfer_from_itin(arr_step, system);
	
	double score = 0;
	Competition_Transfer *ptr;
	
	for(int i = 1; i < 3000; i++) {
		score += get_competition_body_score(transfers[i]);
	}

	printf("\n---\n");
	for(int i = 0; i < 3000; i++) {
		int num_trans = get_num_comp_trans(transfers[i]);
		if(num_trans > 0) printf("%d  %d\n", i, num_trans);
	}

	ptr = transfers[0];

	while(ptr != NULL) {
		printf("%f  %f  %f\n", ptr->rp/AU, ptr->c3/1e9, ptr->epoch);
		ptr = ptr->next;
	}
	printf("--\n");
	for(int i = 1; i < 3000; i++) {
		ptr = transfers[i];
		while(ptr != NULL) {
			printf("%s  %f  %f  %f  | %f  s:%f  f:%f  w:%f |\n",
				   ptr->body->name, mag_vec3(ptr->r)/AU, sqrt(ptr->c3), ptr->rp/ptr->body->radius-1,
				   get_competition_flyby_score(ptr), calc_seasonal_penalty(ptr), calc_flyby_velocity_penalty(ptr),
				   ptr->body->scale_height);
			ptr = ptr->next;
		}
		if(transfers[i] != NULL) printf("%s:  %f\n", transfers[i]->body->name, get_competition_body_score(transfers[i]));
	}
	
	
	ptr = transfers[0];
	
	bool was_below_0_05 = false;
	while(ptr != NULL) {
		if(ptr->rp/AU < 0.05) {
			if(was_below_0_05 || ptr->rp/AU < 0.01) { free_competition_transfer_list(transfers, 3000); return; }
			else was_below_0_05 = true;
		}
		ptr = ptr->next;
	}
	
	for(int i = 1; i < 3000; i++) {
		ptr = transfers[i];
		while(ptr != NULL) {
			if(ptr->rp > 0) {
				if(ptr->rp/ptr->body->radius-1 < 0.1 || ptr->rp/ptr->body->radius-1 > 100) { free_competition_transfer_list(transfers, 3000); return; }
			}
			ptr = ptr->next;
		}
	}
	
	
	printf("--\n");
	printf("SCORE: %f\n", score);
	printf("--\n");
	
	free_competition_transfer_list(transfers, 3000);
}

void run_competition_calc(char *load_filename, char *store_filename, CelestSystem *system) {
	struct Itin_Calc_Data calc_data;
	struct Itin_Calc_Results ic_results;
	
	FILE *file = fopen(load_filename, "r");
	if(!file) {
		perror("Failed to open file");
		return;
	}
	char line[256];  // Buffer for each line
	// skip first line
	fgets(line, sizeof(line), file);
	fgets(line, sizeof(line), file);
	
	enum CompetitionCalcType {CCT_ItinFromT0, CCT_ItinFromFb, CCT_Sequence};
	
	int dep_body_id, arr_body_id, calc_type;
	
	sscanf(line,
		   "%lf,%lf,%lf,%lf,%lf,%lf,%d",
		   &calc_data.jd_min_dep,
		   &calc_data.jd_max_dep,
		   &calc_data.step_dep_date,
		   &calc_data.jd_max_arr,
		   &calc_data.max_duration,
		   &calc_data.dv_filter.max_totdv,
		   &calc_type
	);
	
	Body **fly_by_bodies;
	int counter = 0;
	
	switch(calc_type) {
		case CCT_ItinFromT0:
			fgets(line, sizeof(line), file);
			fgets(line, sizeof(line), file);
			sscanf(line, "%d,%d", &dep_body_id, &arr_body_id);
			
			fly_by_bodies = malloc(system->num_bodies * sizeof(Body*));
			fgets(line, sizeof(line), file);
			while(fgets(line, sizeof(line), file)) {
				int body_id;
				sscanf(line, "%d", &body_id);
				fly_by_bodies[counter++] = get_body_by_id(body_id, system);
			}
			
			
			struct ItinSequenceInfoToTarget seq_info_tt = {
					.system = system,
					.dep_body = get_body_by_id(dep_body_id, system),
					.arr_body = get_body_by_id(arr_body_id, system),
					.num_flyby_bodies = counter,
			};
			
			seq_info_tt.flyby_bodies = fly_by_bodies;
			calc_data.seq_info.to_target = seq_info_tt;
			break;
		case CCT_ItinFromFb:
			break;
		case CCT_Sequence:
			fly_by_bodies = malloc(100 * sizeof(Body*));
			fgets(line, sizeof(line), file);
			while(fgets(line, sizeof(line), file)) {
				int body_id;
				sscanf(line, "%d", &body_id);
				fly_by_bodies[counter++] = get_body_by_id(body_id, system);
			}
			
			
			struct ItinSequenceInfoSpecItin seq_info_spec;
			seq_info_spec.type = ITIN_SEQ_INFO_SPEC_SEQ;
			seq_info_spec.num_steps = counter;
			seq_info_spec.bodies = fly_by_bodies;
			seq_info_spec.system = system;
			calc_data.seq_info.spec_seq = seq_info_spec;
			
	}
	
	fclose(file);
	
	
	calc_data.dv_filter.max_depdv = calc_data.dv_filter.max_totdv;
	calc_data.dv_filter.max_satdv = calc_data.dv_filter.max_totdv;
	calc_data.dv_filter.last_transfer_type = TF_FLYBY;
	
	calc_data.dv_filter.dep_periapsis = 1e9;
	calc_data.dv_filter.arr_periapsis = 1e9;
	
	calc_data.num_deps_per_date = 500;
	calc_data.max_num_waiting_orbits = 0;
	
	ic_results = search_for_itineraries(calc_data);
	
	
	if(ic_results.departures == NULL || ic_results.num_deps == 0) return;
	store_itineraries_in_bfile(ic_results.departures, ic_results.num_nodes, ic_results.num_deps, ic_results.num_itins,
							   calc_data, system, store_filename, get_current_bin_file_type());
	for(int i = 0; i < ic_results.num_deps; i++) free_itinerary(ic_results.departures[i]);
	free(ic_results.departures);
	free(fly_by_bodies);
	if(ic_results.num_deps == 0) printf("No itineraries found!");
}

void store_competition_flyby_arc_arrival(FILE *file, struct ItinStep *step) {
	Vector3 r = scale_vec3(step->r, 1e-3);
	Vector3 v = scale_vec3(step->v_arr, 1e-3);
	Vector3 v_inf = subtract_vec3(v, scale_vec3(step->v_body, 1e-3));
	fprintf(file,
			"%d, "		// body_id: unique identifier for the body
			"%d, "		// flag: status or type flag
			"%lf, "		// epoch: seconds since reference epoch
			"%.9lf, "	// pos_x: position in X-axis (km)
			"%.9lf, "	// pos_y: position in Y-axis (km)
			"%.9lf, "	// pos_z: position in Z-axis (km)
			"%.12lf, "	// vel_x: velocity along X-axis (km/s)
			"%.12lf, "	// vel_y: velocity along Y-axis (km/s)
			"%.12lf, "	// vel_z: velocity along Z-axis (km/s)
			"%.9lf, "	// control_x: control input along X-axis
			"%.9lf, "	// control_y: control input along Y-axis
			"%.9lf\n",	// control_z: control input along Z-axis
			step->body->id, 0, step->date*86400, r.x, r.y, r.z, v.x, v.y, v.z, v_inf.x, v_inf.y, v_inf.z
	);
}

void store_competition_flyby_arc_departure(FILE *file, struct ItinStep *step0) {
	if(step0 == NULL || step0->next == NULL) return;
	struct ItinStep *step1 = step0->next[0];
	Vector3 r = scale_vec3(step0->r, 1e-3);
	Vector3 v = scale_vec3(step1->v_dep, 1e-3);
	Vector3 v_inf = subtract_vec3(v, scale_vec3(step0->v_body, 1e-3));
	fprintf(file,
			"%d, "		// body_id: unique identifier for the body
			"%d, "		// flag: status or type flag
			"%lf, "		// epoch: seconds since reference epoch
			"%.9lf, "	// pos_x: position in X-axis (km)
			"%.9lf, "	// pos_y: position in Y-axis (km)
			"%.9lf, "	// pos_z: position in Z-axis (km)
			"%.12lf, "	// vel_x: velocity along X-axis (km/s)
			"%.12lf, "	// vel_y: velocity along Y-axis (km/s)
			"%.12lf, "	// vel_z: velocity along Z-axis (km/s)
			"%.9lf, "	// control_x: control input along X-axis
			"%.9lf, "	// control_y: control input along Y-axis
			"%.9lf\n",	// control_z: control input along Z-axis
			step0->body->id, 0, step0->date*86400, r.x, r.y, r.z, v.x, v.y, v.z, v_inf.x, v_inf.y, v_inf.z
	);
}

void store_competition_conic_arc(FILE *file, struct ItinStep *step0) {
	if(step0 == NULL || step0->next == NULL) return;
	struct ItinStep *step1 = step0->next[0];
	Vector3 r = scale_vec3(step0->r, 1e-3);
	Vector3 v = scale_vec3(step1->v_dep, 1e-3);
	fprintf(file,
			"%d, "		// body_id: unique identifier for the body
			"%d, "		// flag: status or type flag
			"%lf, "		// epoch: seconds since reference epoch
			"%.9lf, "	// pos_x: position in X-axis (km)
			"%.9lf, "	// pos_y: position in Y-axis (km)
			"%.9lf, "	// pos_z: position in Z-axis (km)
			"%.12lf, "	// vel_x: velocity along X-axis (km/s)
			"%.12lf, "	// vel_y: velocity along Y-axis (km/s)
			"%.12lf, "	// vel_z: velocity along Z-axis (km/s)
			"%lf, "		// control_x: control input along X-axis
			"%lf, "		// control_y: control input along Y-axis
			"%lf\n",	// control_z: control input along Z-axis
			0, 0, step0->date*86400, r.x, r.y, r.z, v.x, v.y, v.z, 0.0, 0.0, 0.0
	);
	r = scale_vec3(step1->r, 1e-3);
	v = scale_vec3(step1->v_arr, 1e-3);
	fprintf(file,
			"%d, "		// body_id: unique identifier for the body
			"%d, "		// flag: status or type flag
			"%lf, "		// epoch: time or timestamp (e.g., seconds since reference epoch)
			"%.9lf, "	// pos_x: position in X-axis (km)
			"%.9lf, "	// pos_y: position in Y-axis (km)
			"%.9lf, "	// pos_z: position in Z-axis (km)
			"%.12lf, "	// vel_x: velocity along X-axis (km/s)
			"%.12lf, "	// vel_y: velocity along Y-axis (km/s)
			"%.12lf, "	// vel_z: velocity along Z-axis (km/s)
			"%lf, "		// control_x: control input along X-axis
			"%lf, "		// control_y: control input along Y-axis
			"%lf\n",	// control_z: control input along Z-axis
			0, 0, step1->date*86400, r.x, r.y, r.z, v.x, v.y, v.z, 0.0, 0.0, 0.0
	);
}

void store_competition_solution(char *filepath, struct ItinStep *step) {
	step = get_first(step);
	
	// Check if the string ends with ".itin"
	if(strlen(filepath) >= 3 && strcmp(filepath + strlen(filepath) - 4, ".csv") != 0) {
		// If not, append ".itin" to the string
		strcat(filepath, ".csv");
	}
	
	FILE *file;
	file = fopen(filepath, "w");
	
	fprintf(file, "#body_id, flag, epoch, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z, control_x, control_y, control_z\n");
	
	step = attach_initial_competition_state(step);
	store_competition_conic_arc(file, step);
	step = step->next[0];
	
	while(step->next != NULL) {
		if(step->prev != NULL) store_competition_flyby_arc_arrival(file, step);
		store_competition_flyby_arc_departure(file, step);
		store_competition_conic_arc(file, step);
		step = step->next[0];
	}
	store_competition_flyby_arc_arrival(file, step);
	
	step = get_first(step)->next[0];
	free(step->prev->next);
	free(step->prev);
	step->prev = NULL;
	
	fclose(file);
}