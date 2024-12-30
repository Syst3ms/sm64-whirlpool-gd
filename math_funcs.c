#define _USE_MATH_DEFINES

#include <math.h>
#include <float.h>

#include "util.h"
#include "math_funcs.h"

#define S 28.0
#define D_EPS cbrt(DBL_EPSILON) 
#define PENALTY_FACTOR 100000
#define MAX_YAW_SPEED 75 * M_PI / 128

// the angle that mario should have to follow the path
double theta(double x, double xp, double z, double zp) {
    double norm = hypot(x, z);
    double yaw_offset = - M_PI * 250 / (norm + 1000);
    double sin_o = sin(yaw_offset),
        cos_o = cos(yaw_offset),
        fac = -20 * (1/norm - 1/2000.0f);
    double current_x = fac * (cos_o * x + sin_o * z),
        current_z = fac * (-sin_o * x + cos_o * z);
    double a = S * zp,
        b = -S * xp,
        c = current_x * zp - current_z * xp;
    return 2*atan((b + sqrt(a*a + b*b - c*c)) / (a - c));
}

void compute_theta_and_current(
    double x, double xp, double z, double zp,
    double *theta, double *x_current
) {
    double norm = hypot(x, z);
    double yaw_offset = - M_PI * 250 / (norm + 1000);
    double sin_o = sin(yaw_offset),
        cos_o = cos(yaw_offset),
        fac = -20 * (1/norm - 1/2000.0f);
    double x_cur = fac * (cos_o * x + sin_o * z),
           z_cur = fac * (-sin_o * x + cos_o * z);
    double a = S * zp,
        b = -S * xp,
        c = x_cur * zp - z_cur * xp;

    *theta = 2*atan((b + sqrt(a*a + b*b - c*c)) / (a - c));
    *x_current = x_cur;
}

double theta_partial_x(
    double x, double xp, double z, double zp
) {
    return (
        theta(x + D_EPS, xp, z, zp)
        - theta(x - D_EPS, xp, z, zp)
    ) / (2 * D_EPS);
}

double theta_partial_xp(
    double x, double xp, double z, double zp
) {
    return (
        theta(x, xp + D_EPS, z, zp)
        - theta(x, xp - D_EPS, z, zp)
    ) / (2 * D_EPS);
}

double theta_partial_z(
    double x, double xp, double z, double zp
) {
    return (
        theta(x, xp, z + D_EPS, zp)
        - theta(x, xp, z - D_EPS, zp)
    ) / (2 * D_EPS);
}

double theta_partial_zp(
    double x, double xp, double z, double zp
) {
    return (
        theta(x, xp, z, zp + D_EPS)
        - theta(x, xp, z, zp - D_EPS)
    ) / (2 * D_EPS);
}

double theta_p(
    double x, double xp, double xpp,
    double z, double zp, double zpp) {
    return xp * theta_partial_x(x, xp, z, zp) +
           xpp * theta_partial_xp(x, xp, z, zp) +
           zp * theta_partial_z(x, xp, z, zp) +
           zpp * theta_partial_zp(x, xp, z, zp);
}

// a function f(x,x',z,z',u) such that the lagrangian is f(x,x',z,z',θ*'(x,x',z,z'))
double lagr_intermediate(double x, double xp, double z, double zp, double u) {
    double norm = hypot(x, z);
    double yaw_offset = - M_PI * 250 / (norm + 1000);
    double sin_o = sin(yaw_offset),
          cos_o = cos(yaw_offset),
          fac = -20 * (1/norm - 1/2000.0f);
    double current_x = fac * (cos_o * x + sin_o * z),
          current_z = fac * (-sin_o * x + cos_o * z);
    double a = S * zp,
          b = -S * xp,
          c = current_x * zp - current_z * xp;
    double inner = (b + sqrt(a*a + b*b - c*c)) / (a - c);
    double cos_theta = 2 / (1 + inner * inner) - 1;
    double time_integrand = fabs(xp / (S * cos_theta + current_x));
    double penalty_term = fmax(0, (u / time_integrand) / MAX_YAW_SPEED - 1);
    return time_integrand + PENALTY_FACTOR * penalty_term;
}

double lagrangian_with_precomputed(
    double x, double xp, double xpp, double z, double zp, double zpp,
    double x_current, double theta
) {
    double time_integrand = fabs(xp / (S * cos(theta) + x_current));
    double thetap = theta_p(x, xp, xpp, z, zp, zpp);
    double penalty_term = fmax(0, (thetap / time_integrand) / MAX_YAW_SPEED - 1);
    return time_integrand + PENALTY_FACTOR * penalty_term;
}

double thetap_partial_x(
    double x, double xp, double xpp, double z, double zp, double zpp
) {
    return (
        theta_p(x + D_EPS, xp, xpp, z, zp, zpp)
        - theta_p(x - D_EPS, xp, xpp, z, zp, zpp)
    ) / (2 * D_EPS);
}

double thetap_partial_xp(
    double x, double xp, double xpp, double z, double zp, double zpp
) {
    return (
        theta_p(x, xp + D_EPS, xpp, z, zp, zpp)
        - theta_p(x, xp - D_EPS, xpp, z, zp, zpp)
    ) / (2 * D_EPS);
}

double thetap_partial_z(
    double x, double xp, double xpp, double z, double zp, double zpp
) {
    return (
        theta_p(x, xp, xpp, z + D_EPS, zp, zpp)
        - theta_p(x, xp, xpp, z - D_EPS, zp, zpp)
    ) / (2 * D_EPS);
}

double thetap_partial_zp(
    double x, double xp, double xpp, double z, double zp, double zpp
) {
    return (
        theta_p(x, xp, xpp, z, zp + D_EPS, zpp)
        - theta_p(x, xp, xpp, z, zp - D_EPS, zpp)
    ) / (2 * D_EPS);
}

double lagr_inter_partial_theta(double x, double xp, double z, double zp, double u) {
    return (
        lagr_intermediate(x, xp, z, zp, u + D_EPS)
        - lagr_intermediate(x, xp, z, zp, u - D_EPS)
    ) / (2 * D_EPS);
}

double lagr_partial_x(
    double x, double xp, double xpp, double z, double zp, double zpp,
    double u, double partial_u
) {
    double part = (
        lagr_intermediate(x + D_EPS, xp, z, zp, u)
        - lagr_intermediate(x - D_EPS, xp, z, zp, u)
    ) / (2 * D_EPS);
    return part + thetap_partial_x(x, xp, xpp, z, zp, zpp) * partial_u;
}

double lagr_partial_xp(
    double x, double xp, double xpp, double z, double zp, double zpp,
    double u, double partial_u
) {
    double part = (
        lagr_intermediate(x, xp + D_EPS, z, zp, u)
        - lagr_intermediate(x, xp - D_EPS, z, zp, u)
    ) / (2 * D_EPS);
    return part + thetap_partial_xp(x, xp, xpp, z, zp, zpp) * partial_u;
}

double lagr_partial_z(
    double x, double xp, double xpp, double z, double zp, double zpp,
    double u, double partial_u
) {
    double part = (
        lagr_intermediate(x, xp, z + D_EPS, zp, u)
        - lagr_intermediate(x, xp, z - D_EPS, zp, u)
    ) / (2 * D_EPS);
    return part + thetap_partial_z(x, xp, xpp, z, zp, zpp) * partial_u;
}

double lagr_partial_zp(
    double x, double xp, double xpp, double z, double zp, double zpp,
    double u, double partial_u
) {
    double part = (
        lagr_intermediate(x, xp, z, zp + D_EPS, u)
        - lagr_intermediate(x, xp, z, zp - D_EPS, u)
    ) / (2 * D_EPS);
    return part + thetap_partial_zp(x, xp, xpp, z, zp, zpp) * partial_u;
}

double objective(struct data *d) {
    double total = 0.0f;

    for (int i = 0; i < POINTS-1; i++) {
        total += d->lagrangian[i];
    }

    return total / (POINTS-1);
}

void recompute_dependent(struct data *d) {
    struct path *p = &d->p;
    derivative(d->xp, p->x, POINTS);
    derivative(d->zp, p->z, POINTS);
    derivative(d->xpp, d->xp, POINTS-1);
    derivative(d->zpp, d->zp, POINTS-1);
    d->xpp[POINTS-2] = d->xpp[POINTS-3];
    d->zpp[POINTS-2] = d->zpp[POINTS-3];

    for (int i = 0; i < POINTS-1; i++) {
        double x_current;
        double x = p->x[i+1],
               xp = d->xp[i],
               xpp = d->xpp[i],
               z = p->z[i+1],
               zp = d->zp[i],
               zpp = d->zpp[i];
        compute_theta_and_current(
            x, xp, z, zp,
            &d->theta[i], &x_current
        );
        d->partial_theta[i] = lagr_inter_partial_theta(x, xp, z, zp, d->theta[i]);
        d->lagrangian[i] = lagrangian_with_precomputed(
            x, xp, xpp, z, zp, zpp,
            x_current, d->theta[i]
        );
    }
}