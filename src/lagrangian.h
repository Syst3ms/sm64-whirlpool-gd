#pragma once

#include "types.h"

#define LAGR_PARTIAL_HEADER(var)\
    double lagr_partial_##var(union point *pt, double penalty_fac)

LAGR_PARTIAL_HEADER(x);
LAGR_PARTIAL_HEADER(z);
LAGR_PARTIAL_HEADER(xp);
LAGR_PARTIAL_HEADER(zp);

double compute_lagrangian(union point *pt, double penalty_fac);
double time_integrand_alone(v2d pos, v2d vel);
double theta(v2d pos, v2d vel);
double real_speed_norm(v2d pos, v2d vel);