#include <string.h>
#include <stdlib.h>
#include <gtk/gtk.h>
#include "celestial_systems.h"
#include "tools/file_io.h"
#include "competition_tools.h"

CelestSystem **available_systems;

int num_available_systems = 0;


void init_available_systems(char *directory) {
	available_systems = malloc(sizeof(struct System*));
	available_systems[num_available_systems++] = load_competition_system(directory);
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
