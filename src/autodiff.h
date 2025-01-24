#pragma once

#include "types.h"

void compute_lagrangian_with_intermediates(struct pt_vel *pt);
double time_integrand_alone(v2d pos, v2d vel);
void theta_ad1(struct ad1d x, struct ad1d z, struct ad1d xp, struct ad1d zp, struct ad1d *out, double *maxrange_out);