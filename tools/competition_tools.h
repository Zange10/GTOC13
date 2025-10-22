#ifndef KMAT_COMPETITION_TOOLS_H
#define KMAT_COMPETITION_TOOLS_H

#include "orbitlib.h"
#include "orbit_calculator/itin_tool.h"

CelestSystem * load_competition_system(char *directory);

double get_itin_competition_score(struct ItinStep *arr_step, CelestSystem *system);

void print_itin_competition_score(struct ItinStep *arr_step, CelestSystem *system);

#endif //KMAT_COMPETITION_TOOLS_H
