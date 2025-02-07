#pragma once

#include "types.h"

#define LAGR_PARTIAL_HEADER(var) double lagr_partial_##var(union point *pt)

LAGR_PARTIAL_HEADER(x);
LAGR_PARTIAL_HEADER(z);
LAGR_PARTIAL_HEADER(xp);
LAGR_PARTIAL_HEADER(zp);

double lagr_partial_xp_with_side_eff(union point *pt, double *objective);
double compute_lagrangian(union point *pt);
double time_integrand_alone(v2d pos, v2d vel);