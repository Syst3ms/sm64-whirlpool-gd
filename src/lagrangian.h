#pragma once

#include "types.h"

#define LAGR_PARTIAL_HEADER(var)\
    double lagr_partial_##var(union point *pt, double penalty_fac, double shift)

LAGR_PARTIAL_HEADER(x);
LAGR_PARTIAL_HEADER(z);
LAGR_PARTIAL_HEADER(xp);
LAGR_PARTIAL_HEADER(zp);

double compute_lagrangian(union point *pt, double penalty_fac, double shift);
void compute_lagrangian_and_constraint(
    union point *pt, double penalty_fac, double shift,
    double *lagr_out, double *constraint_out
);
double objective_from_positions(v2d *inter_pos, v2d start, v2d end, double penalty_fac);
double time_integrand_alone(v2d pos, v2d vel);
double theta(v2d pos, v2d vel);
double real_speed_norm(v2d pos, v2d vel);