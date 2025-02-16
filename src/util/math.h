#pragma once

#include "../types.h"

unsigned short radians_to_au(double rad);

double fast_hypot(double dx, double dz);
double fast_hypot_v(v2d v);

v2d v2d_of(double a, double b);

v2d ortho(v2d a);
v2d vabs(v2d x);
double square(double x);