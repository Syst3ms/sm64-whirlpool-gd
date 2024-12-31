#define _USE_MATH_DEFINES

#include <math.h>
#include <float.h>

#include "util.h"
#include "math_funcs.h"

#define S 28.0
#define D_EPS cbrt(DBL_EPSILON) 
#define PENALTY_FACTOR 0.001
#define MAX_YAW_SPEED (75*M_PI/128)

// the angle that mario should have to follow the path
double theta(double x, double xp, double z, double zp) {
    double norm = hypot(x, z);
    double yaw_offset = - M_PI * 250 / (norm + 1000);
    double sin_o = sin(yaw_offset),
        cos_o = cos(yaw_offset),
        fac = -20 * (1/norm - 1/2000.0f);
    double whirl_x = fac * (cos_o * x + sin_o * z),
           whirl_z = fac * (-sin_o * x + cos_o * z);
    double speed_norm = hypot(xp, zp),
           det = whirl_x * zp - whirl_z * xp;
    return acos(det / (S * speed_norm)) - atan2(zp, xp);
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
    double whirl_x = fac * (cos_o * x + sin_o * z),
           whirl_z = fac * (-sin_o * x + cos_o * z);
    double speed_norm = hypot(xp, zp),
           det = whirl_x * zp - whirl_z * xp;

    *theta = acos(det / (S * speed_norm)) - atan2(zp, xp);
    *x_current = whirl_x;
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

// a function f(x,x',z,z',u) such that the lagrangian is f(x,x',z,z',θ'(x,x',z,z'))
double lagr_intermediate(double x, double xp, double z, double zp, double u) {
    double norm = hypot(x, z);
    double yaw_offset = - M_PI * 250 / (norm + 1000);
    double sin_o = sin(yaw_offset),
          cos_o = cos(yaw_offset),
          fac = -20 * (1/norm - 1/2000.0f);
    double whirl_x = fac * (cos_o * x + sin_o * z),
           whirl_z = fac * (-sin_o * x + cos_o * z);
    double a = (whirl_x * zp - whirl_z * xp) / S,
           b = xp * xp + zp * zp;
    double sin_theta = (xp * sqrt(b - a * a) - a * zp) / b;
    double time_integrand = fabs(xp / (S * sin_theta + whirl_x)) / (POINTS-1);
    double penalty_term = fmax(0, (fabs(u) / time_integrand) / MAX_YAW_SPEED - 1);
    return time_integrand + PENALTY_FACTOR * penalty_term;
}

double lagrangian_with_precomputed(
    double xp, double x_current, double theta, double theta_p
) {
    double time_integrand = fabs(xp / (S * sin(theta) + x_current)) / (POINTS-1);
    double penalty_term = fmax(0, (fabs(theta_p) / time_integrand) / MAX_YAW_SPEED - 1);
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
    double theta_p, double partial_u
) {
    double part = (
        lagr_intermediate(x + D_EPS, xp, z, zp, theta_p)
        - lagr_intermediate(x - D_EPS, xp, z, zp, theta_p)
    ) / (2 * D_EPS);
    return part + thetap_partial_x(x, xp, xpp, z, zp, zpp) * partial_u;
}

double lagr_partial_xp(
    double x, double xp, double xpp, double z, double zp, double zpp,
    double theta_p, double partial_u
) {
    double part = (
        lagr_intermediate(x, xp + D_EPS, z, zp, theta_p)
        - lagr_intermediate(x, xp - D_EPS, z, zp, theta_p)
    ) / (2 * D_EPS);
    return part + thetap_partial_xp(x, xp, xpp, z, zp, zpp) * partial_u;
}

double lagr_partial_z(
    double x, double xp, double xpp, double z, double zp, double zpp,
    double theta_p, double partial_u
) {
    double part = (
        lagr_intermediate(x, xp, z + D_EPS, zp, theta_p)
        - lagr_intermediate(x, xp, z - D_EPS, zp, theta_p)
    ) / (2 * D_EPS);
    return part + thetap_partial_z(x, xp, xpp, z, zp, zpp) * partial_u;
}

double lagr_partial_zp(
    double x, double xp, double xpp, double z, double zp, double zpp,
    double theta_p, double partial_u
) {
    double part = (
        lagr_intermediate(x, xp, z, zp + D_EPS, theta_p)
        - lagr_intermediate(x, xp, z, zp - D_EPS, theta_p)
    ) / (2 * D_EPS);
    return part + thetap_partial_zp(x, xp, xpp, z, zp, zpp) * partial_u;
}

double objective(struct data *d) {
    double total = 0.0f;

    for (int i = 0; i < POINTS-1; i++) {
        total += d->lagrangian[i];
    }

    return total;
}

void recompute_dependent(struct data *d) {
    struct path *p = &d->p;
    derivative(d->xp, p->x, POINTS);
    derivative(d->zp, p->z, POINTS);
    derivative(d->lambda_p, d->lambda, POINTS);
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
        d->theta_p[i] = theta_p(x, xp, xpp, z, zp, zpp);
        d->partial_theta[i] = lagr_inter_partial_theta(
            x, xp, z, zp, d->theta_p[i]
        );
        d->lagrangian[i] = lagrangian_with_precomputed(
            xp, x_current, d->theta[i], d->theta_p[i]
        );
    }
}