#include "poly_sphere.h"
#include "part_sticky_spheres.h"
#include "no_field.h"
#include "no_boundary.h"
#include "JRand_dyn.h"
#include "box__neighbour_list__cell_list.h"
#include "integrator_MC_NPT_swap.h"

#include <iostream>
#include <sstream>
#include <cstring>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#define DIM 3

#define NUMBER_PARTICLES 1000
#define TEMPERATURE 1.0

#define BOX_SIZE_X 8.5
#define BOX_SIZE_Y 8.5
#define BOX_SIZE_Z 8.5

#define PB_X 1
#define PB_Y 1
#define PB_Z 1

#define DR 0.1
#define DV 0.1

#define WAITING_TIME_MAX 0.0
#define SIMULATION_TIME_MAX 100.0
