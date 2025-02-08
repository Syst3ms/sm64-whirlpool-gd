#pragma once

#include "types.h"

double square(double x);
v2d dot(v2d *a, v2d *b, size_t len);
double theta(v2d pos, v2d vel);
double real_speed_norm(v2d pos, v2d vel);