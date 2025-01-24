#pragma once

#include "types.h"

double lagr_partial_x(struct pt_vel *pt);
double lagr_partial_z(struct pt_vel *pt);
double lagr_partial_xp(struct pt_vel *pt);
double lagr_partial_zp(struct pt_vel *pt);

double theta(v2d pos, v2d vel);
double theta_maxrange(v2d pos, v2d vel, double *theta);
double real_speed_norm(v2d pos, v2d vel);