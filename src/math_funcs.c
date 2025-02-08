#define _USE_MATH_DEFINES

#include <float.h>
#include <math.h>
#include <stddef.h>

#include "parameters.h"
#include "util.h"
#include "math_funcs.h"

#define D_EPS cbrt(DBL_EPSILON) 

v2d dot(v2d *a, v2d *b, size_t len) {
    v2d res = {};

    for (size_t i = 0; i < len; i++) {
        res += a[i] * b[i];
    }

    return res;
}

double square(double x) {
    return x * x;
}

// the angle that mario should have to follow the path
double theta(v2d pos, v2d vel) {
    double x = pos[0], z = pos[1], xp = vel[0], zp = vel[1];
    double norm = fast_hypot(x, z),
           vel_sq_norm = xp * xp + zp * zp;
    double yaw_offset = - M_PI * 250 / (norm + 1000);
    double sin_o = sin(yaw_offset),
          cos_o = cos(yaw_offset),
          fac = -20 * (1/norm - 1/2000.0f) / (S * vel_sq_norm);
    double whirl_x_no_fac = cos_o * x + sin_o * z,
           whirl_z_no_fac = -sin_o * x + cos_o * z;
    double det = fac * (whirl_x_no_fac * zp - whirl_z_no_fac * xp);
    return acos(det) - atan2(zp, xp);
}

double real_speed_norm(v2d pos, v2d vel) {
    double x = pos[0], z = pos[1], xp = vel[0], zp = vel[1];
    double norm = fast_hypot(x, z);
    double yaw_offset = - M_PI * 250 / (norm + 1000);
    double sin_o = sin(yaw_offset),
          cos_o = cos(yaw_offset),
          fac = -20 * (1/norm - 1/2000.0f);
    double whirl_x = fac * (cos_o * x + sin_o * z),
           whirl_z = fac * (-sin_o * x + cos_o * z);
    double a = (whirl_x * zp - whirl_z * xp) / S,
           b = xp * xp + zp * zp;
    double sin_theta = (xp * sqrt(b - a * a) - a * zp) / b,
           cos_theta = (zp * sqrt(b - a * a) + a * xp) / b;
    return fast_hypot(S * sin_theta + whirl_x, S * cos_theta + whirl_z);
}