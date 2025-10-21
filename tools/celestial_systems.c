#include <string.h>
#include <stdlib.h>
#include <gtk/gtk.h>
#include "celestial_systems.h"
#include "tools/file_io.h"

CelestSystem **available_systems;

int num_available_systems = 0;

int parse_planet_line(const char *line, Body *new_body) {
	// Ignore comment or empty lines
	if (line[0] == '#' || strlen(line) < 3)
		return 0;
	
	// Use sscanf to parse CSV values
	int n = sscanf(line,
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
	
	new_body->mu *= 1e9;
	new_body->radius *= 1e3;
	new_body->orbit.a *= 1e3;
	new_body->orbit.i = deg2rad(new_body->orbit.i);
	new_body->orbit.raan = deg2rad(new_body->orbit.raan);
	new_body->orbit.arg_peri = deg2rad(new_body->orbit.arg_peri);
	new_body->orbit.ta = calc_true_anomaly_from_mean_anomaly(new_body->orbit, deg2rad(new_body->orbit.ta));
	
	return (n == 11); // Success if all fields parsed
}

int parse_comet_line(const char *line, Body *new_body) {
	// Ignore comment or empty lines
	if (line[0] == '#' || strlen(line) < 3)
		return 0;
	
	// Use sscanf to parse CSV values
	int n = sscanf(line,
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
	
	sprintf(new_body->name, "C%d", new_body->id);
	new_body->color[0] = 0;
	new_body->color[1] = 0.3;
	new_body->color[2] = 1;
	
	new_body->mu *= 1e9;
	new_body->radius *= 1e3;
	new_body->orbit.a *= 1e3;
	new_body->orbit.i = deg2rad(new_body->orbit.i);
	new_body->orbit.raan = deg2rad(new_body->orbit.raan);
	new_body->orbit.arg_peri = deg2rad(new_body->orbit.arg_peri);
	new_body->orbit.ta = calc_true_anomaly_from_mean_anomaly(new_body->orbit, deg2rad(new_body->orbit.ta));
	
	return (n == 11); // Success if all fields parsed
}

int parse_asteroid_line(const char *line, Body *new_body) {
	// Ignore comment or empty lines
	if (line[0] == '#' || strlen(line) < 3)
		return 0;
	
	// Use sscanf to parse CSV values
	int n = sscanf(line,
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
	
	sprintf(new_body->name, "A%d", new_body->id);
	new_body->color[0] = 0.5;
	new_body->color[1] = 0.5;
	new_body->color[2] = 0.5;
	
	new_body->mu *= 1e9;
	new_body->radius *= 1e3;
	new_body->orbit.a *= 1e3;
	new_body->orbit.i = deg2rad(new_body->orbit.i);
	new_body->orbit.raan = deg2rad(new_body->orbit.raan);
	new_body->orbit.arg_peri = deg2rad(new_body->orbit.arg_peri);
	new_body->orbit.ta = calc_true_anomaly_from_mean_anomaly(new_body->orbit, deg2rad(new_body->orbit.ta));
	
	return (n == 11); // Success if all fields parsed
}

CelestSystem * test_load(char *filename) {
	CelestSystem *system = new_system();
	system->num_bodies = 0;
	system->prop_method = ORB_ELEMENTS;
	sprintf(system->name, "Altaira System");
	system->ut0 = 0;
	
	system->cb = new_body();
	sprintf(system->cb->name, "Altaira");
	system->cb->color[0] = 1;
	system->cb->color[1] = 1;
	system->cb->color[2] = 0.3;
	system->cb->mu = 139348062043.343e9;
	system->cb->system = system;
	system->bodies = malloc(10000*sizeof(Body*));
	
	FILE *file = fopen(filename, "r");
	if (!file) {
		perror("Failed to open file");
		free(system);
		return NULL;
	}
	char line[256];  // Buffer for each line
	// skip first line
	fgets(line, sizeof(line), file);
	
	while (fgets(line, sizeof(line), file)) {
		Body *body = new_body();
		int success = parse_planet_line(line, body);
		printf("%d\n", success);
		body->orbit.cb = system->cb;
		system->bodies[system->num_bodies] = body;
		system->num_bodies++;
	}
	system->home_body = system->bodies[0];
	fclose(file);
	
//	file = fopen("../Celestial_Systems/gtoc13_comets.csv", "r");
//	if (!file) {
//		perror("Failed to open file");
//		free(system);
//		return NULL;
//	}
//
//	// skip first line
//	fgets(line, sizeof(line), file);
//
//	while (fgets(line, sizeof(line), file)) {
//		Body *body = new_body();
//		int success = parse_comet_line(line, body);
//		printf("%d\n", success);
//		body->orbit.cb = system->cb;
//		system->bodies[system->num_bodies] = body;
//		system->num_bodies++;
//	}
//	system->home_body = system->bodies[0];
//	fclose(file);

//	file = fopen("../Celestial_Systems/gtoc13_asteroids.csv", "r");
//	if (!file) {
//		perror("Failed to open file");
//		free(system);
//		return NULL;
//	}
//
//	// skip first line
//	fgets(line, sizeof(line), file);
//
//	while (fgets(line, sizeof(line), file)) {
//		Body *body = new_body();
//		int success = parse_asteroid_line(line, body);
//		printf("%d\n", success);
//		body->orbit.cb = system->cb;
//		system->bodies[system->num_bodies] = body;
//		system->num_bodies++;
//	}
//	system->home_body = system->bodies[0];
//	fclose(file);
	
	return system;
}


void init_available_systems(const char *directory) {
//	available_systems = init_available_systems_from_path(directory, &num_available_systems);
	
	available_systems = malloc(sizeof(struct System*));
	available_systems[num_available_systems++] = test_load("../Celestial_Systems/gtoc13_planets.csv");
}

int get_num_available_systems() {return num_available_systems;}
CelestSystem ** get_available_systems() {return available_systems;}

int is_available_system(CelestSystem *system) {
	if(system == NULL) return 0;
	for(int i = 0; i < num_available_systems; i++) {
		if(system == get_available_systems()[i]) return 1;
	}
	return 0;
}

CelestSystem * get_subsystem_from_system_and_id_rec(CelestSystem *system, int *id) {
	if(*id == 0) return system;
	(*id)--;
	for(int i = 0; i < system->num_bodies; i++) {
		if(system->bodies[i]->system != NULL) {
			CelestSystem *subsystem = get_subsystem_from_system_and_id_rec(system->bodies[i]->system, id);
			if(subsystem != NULL) return subsystem;
		}
	}
	return NULL;
}

CelestSystem * get_subsystem_from_system_and_id(CelestSystem *system, int id) {
	return get_subsystem_from_system_and_id_rec(system, &id);
}

CelestSystem * get_system_by_name(char *name) {
	for(int i = 0; i < get_num_available_systems(); i++) {
		if(strcmp(get_available_systems()[i]->name, name) == 0) return get_available_systems()[i];
	}
	
	return NULL;
}

void free_all_celestial_systems() {
	free_celestial_systems(available_systems, num_available_systems);
}
