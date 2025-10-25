#ifndef KMAT_COMPETITION_TOOLS_H
#define KMAT_COMPETITION_TOOLS_H

#include "orbitlib.h"
#include "orbit_calculator/itin_tool.h"

#define AU 149597870691.0

CelestSystem * load_competition_system(char *directory);

double get_itin_competition_score(struct ItinStep *arr_step, CelestSystem *system);

void print_itin_competition_score(struct ItinStep *arr_step, CelestSystem *system);

Vector3 calc_heliocentric_periapsis(Vector3 r_dep, Vector3 v_dep, Vector3 r_arr, Vector3 v_arr, CelestSystem *system);

void run_competition_calc(char *load_filename, char *store_filename, CelestSystem *system);

void store_competition_solution(char *filepath, struct ItinStep *step);

struct ItinStep * attach_initial_competition_state(struct ItinStep *step);

#endif //KMAT_COMPETITION_TOOLS_H
