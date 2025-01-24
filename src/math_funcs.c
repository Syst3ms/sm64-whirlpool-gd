#define _USE_MATH_DEFINES

#include <float.h>
#include <math.h>

#include "parameters.h"
#include "util.h"
#include "math_funcs.h"

#define D_EPS cbrt(DBL_EPSILON) 

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

// the angle that mario should have to follow the path
double theta(v2d pos, v2d vel) {
    double x = pos[0], z = pos[1], xp = vel[0], zp = vel[1];
    double norm = fast_hypot(x, z),
            vel_norm = fast_hypot(xp, zp);
    double yaw_offset = - M_PI * 250 / (norm + 1000);
    double sin_o = sin(yaw_offset),
            cos_o = cos(yaw_offset),
            fac = (-20 / norm + 0.01) / (S * vel_norm);
    double whirl_x_no_fac = cos_o * x + sin_o * z,
            whirl_z_no_fac = -sin_o * x + cos_o * z;
    double det = fac * (whirl_x_no_fac * zp - whirl_z_no_fac * xp);
    return acos(det) - atan2(zp, xp);
}

// the angle that mario should have to follow the path
double theta_maxrange(v2d pos, v2d vel, double *theta) {
    double x = pos[0], z = pos[1], xp = vel[0], zp = vel[1];
    double norm = fast_hypot(x, z),
            vel_norm = fast_hypot(xp, zp);
    double yaw_offset = - M_PI * 250 / (norm + 1000);
    double sin_o = sin(yaw_offset),
            cos_o = cos(yaw_offset),
            fac = (-20 / norm + 0.01);
    double whirl_x_no_fac = cos_o * x + sin_o * z,
            whirl_z_no_fac = -sin_o * x + cos_o * z;
    double det = fac * (whirl_x_no_fac * zp - whirl_z_no_fac * xp);
    double t = *theta = acos(det / (S * vel_norm)) - atan2(zp, xp);
    return MAX_YAW_SPEED / POINTS * fabs(xp / (S * sin(t) + whirl_x_no_fac * fac));
}

double compute_lagrangian(struct pt_vel *pt) {
    double x = pt->x, z = pt->z, xp = pt->xp, zp = pt->zp;
    double norm = fast_hypot(x, z);
    double yaw_offset = - M_PI * 250 / (norm + 1000);
    double sin_o = sin(yaw_offset),
            cos_o = cos(yaw_offset),
            fac = -20 / norm + 0.01;
    double whirl_x_nf = cos_o * x + sin_o * z,
            whirl_z_nf = -sin_o * x + cos_o * z;
    double a = fac * (whirl_x_nf * zp - whirl_z_nf * xp) / S,
            b = xp * xp + zp * zp;
    double sin_theta = (xp * sqrt(b - a * a) - a * zp) / b;
    return fabs(xp / (S * sin_theta + fac * whirl_x_nf));
}

#define LAGR_PARTIAL(var)\
    double lagr_partial_##var(struct pt_vel *pt) {\
        double orig = pt->var;\
        pt->var *= D_FAC_UP;\
        double lagr_up = compute_lagrangian(pt);\
        pt->var = orig;\
        pt->var *= D_FAC_DOWN;\
        double lagr_down = compute_lagrangian(pt);\
        pt->var = orig;\
        return (lagr_up - lagr_down) / (2 * D_EPS * orig);\
    }\

LAGR_PARTIAL(x)
LAGR_PARTIAL(z)
LAGR_PARTIAL(xp)
LAGR_PARTIAL(zp)