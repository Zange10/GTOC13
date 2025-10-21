#ifndef CELESTIAL_BODIES
#define CELESTIAL_BODIES

#include "orbitlib.h"

void init_available_systems(const char *directory);

int get_num_available_systems();

CelestSystem ** get_available_systems();

int is_available_system(CelestSystem *system);

CelestSystem * get_subsystem_from_system_and_id(CelestSystem *system, int id);

CelestSystem * get_system_by_name(char *name);

void free_all_celestial_systems();

#endif