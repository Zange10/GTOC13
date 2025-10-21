#include "competition_tools.h"
#include <stdio.h>
#include <string.h>

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
	sprintf(system->cb->name, "Altaira");
	system->cb->color[0] = 1;
	system->cb->color[1] = 1;
	system->cb->color[2] = 0.3;
	system->cb->mu = 139348062043.343e9;
	system->cb->system = system;
	system->bodies = malloc(10000*sizeof(Body*));
	
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