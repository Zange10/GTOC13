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
	if (line[0] == '#' || strlen(line) < 3)
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
	if (!file) {
		perror("Failed to open file");
		free(system);
		return;
	}
	char line[256];  // Buffer for each line
	// skip first line
	fgets(line, sizeof(line), file);
	
	while (fgets(line, sizeof(line), file)) {
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

void free_competition_transfer_list(Competition_Transfer ** transfers, int num_bodies) {
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
				itin_ptr->r,itin_ptr->v_arr,
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

double get_competition_body_score(Competition_Transfer * transfer) {
	if(transfer == 0) return 0;
	double w = transfer->body->scale_height;
	double score = 0;
	while(transfer != NULL) {
		double v_inf = sqrt(transfer->c3)/1e3;
		double s_sum = 0;
		Competition_Transfer *prev_ptr = transfer->prev;
		while(prev_ptr != NULL) {
			Vector3 r_i = norm_vec3(transfer->r);
			Vector3 r_j = norm_vec3(prev_ptr->r);
			double acosd = rad2deg(acos(dot_vec3(r_i, r_j)));
			s_sum += exp(-(acosd*acosd)/50);
			prev_ptr = prev_ptr->prev;
		}
		double s = 0.1 + (0.9/(1+10*s_sum));
		double f = 0.2 + exp(-v_inf/13) / (1+exp(-5*(v_inf-1.5)));
		score += s*f;
		transfer = transfer->next;
	}
	score *= w;
	return score;
}

double get_itin_competition_score(struct ItinStep *arr_step, CelestSystem *system) {
	double b = is_grand_tour() ? 1.2 : 1;
	
	Competition_Transfer ** transfers = build_competition_transfer_from_itin(arr_step, system);
	
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
				if(ptr->rp/ptr->body->radius-1 < 0.1 || ptr->rp/ptr->body->radius-1 > 100) { free_competition_transfer_list(transfers, 3000); return 0; }
			}
			ptr = ptr->next;
		}
	}
	
	free_competition_transfer_list(transfers, 3000);
	
	return score;
}

void print_itin_competition_score(struct ItinStep *arr_step, CelestSystem *system) {
	double b = is_grand_tour() ? 1.2 : 1;
	
	Competition_Transfer ** transfers = build_competition_transfer_from_itin(arr_step, system);
	
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
			printf("%s (%f)  %f  %f  %f\n", ptr->body->name, ptr->body->scale_height, mag_vec3(ptr->r)/AU, sqrt(ptr->c3), ptr->rp/ptr->body->radius-1);
			ptr = ptr->next;
		}
		if(transfers[i] != NULL) printf("%s (%f)  %f\n", transfers[i]->body->name, transfers[i]->body->scale_height, get_competition_body_score(transfers[i]));
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
	struct Itin_Calc_Data ic_calc_data;
	struct Itin_Calc_Results ic_results;
	
	FILE *file = fopen(load_filename, "r");
	if (!file) {
		perror("Failed to open file");
		return;
	}
	char line[256];  // Buffer for each line
	// skip first line
	fgets(line, sizeof(line), file);
	fgets(line, sizeof(line), file);
	
	int dep_body_id, arr_body_id;
	
	sscanf(line,
		   "%lf,%lf,%lf,%lf,%lf,%d,%d",
		   &ic_calc_data.jd_min_dep,
		   &ic_calc_data.jd_max_dep,
		   &ic_calc_data.jd_max_arr,
		   &ic_calc_data.max_duration,
		   &ic_calc_data.dv_filter.max_totdv,
		   &dep_body_id,
		   &arr_body_id
	);
	
	
	Body **fly_by_bodies = malloc(system->num_bodies * sizeof(Body*));
	fgets(line, sizeof(line), file);
	int idx = 0;
	while(fgets(line, sizeof(line), file)) {
		int body_id;
		sscanf(line, "%d", &body_id);
		fly_by_bodies[idx++] = get_body_by_id(body_id, system);
	}
	
	fclose(file);
	
	
	ic_calc_data.dv_filter.max_depdv = ic_calc_data.dv_filter.max_totdv;
	ic_calc_data.dv_filter.max_satdv = ic_calc_data.dv_filter.max_totdv;
	ic_calc_data.dv_filter.last_transfer_type = TF_FLYBY;
	
	ic_calc_data.dv_filter.dep_periapsis = 1e9;
	ic_calc_data.dv_filter.arr_periapsis = 1e9;
	
	ic_calc_data.num_deps_per_date = 500;
	ic_calc_data.step_dep_date = 1;
	ic_calc_data.max_num_waiting_orbits = 0;
	
	struct ItinSequenceInfoToTarget seq_info = {
			.system = system,
			.dep_body = get_body_by_id(dep_body_id, system),
			.arr_body = get_body_by_id(arr_body_id, system),
			.num_flyby_bodies = idx,
	};
	
	seq_info.flyby_bodies = fly_by_bodies;
	ic_calc_data.seq_info.to_target = seq_info;
	
	ic_results = search_for_itineraries(ic_calc_data);
	
	
	if(ic_results.departures == NULL || ic_results.num_deps == 0) return;
	store_itineraries_in_bfile(ic_results.departures, ic_results.num_nodes, ic_results.num_deps, ic_results.num_itins, ic_calc_data, system, store_filename, get_current_bin_file_type());
	for(int i = 0; i < ic_results.num_deps; i++) free_itinerary(ic_results.departures[i]);
	free(ic_results.departures);
	free(fly_by_bodies);
	if(ic_results.num_deps == 0) printf("No itineraries found!");
}