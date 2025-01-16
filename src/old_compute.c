#define _USE_MATH_DEFINES

#include <float.h>
#include <math.h>

#include "parameters.h"
#include "util.h"
#include "math_funcs.h"

#define D_EPS cbrt(DBL_EPSILON) 

// FOR DEBUG PURPOSES

// the angle that mario should have to follow the path
double theta_old(double x, double xp, double z, double zp) {
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

double theta_partial_x_old(
    double x, double xp, double z, double zp
) {
    return (
        theta_old(x + D_EPS, xp, z, zp)
        - theta_old(x - D_EPS, xp, z, zp)
    ) / (2 * D_EPS);
}

double theta_partial_xp_old(
    double x, double xp, double z, double zp
) {
    return (
        theta_old(x, xp + D_EPS, z, zp)
        - theta_old(x, xp - D_EPS, z, zp)
    ) / (2 * D_EPS);
}

double theta_partial_z_old(
    double x, double xp, double z, double zp
) {
    return (
        theta_old(x, xp, z + D_EPS, zp)
        - theta_old(x, xp, z - D_EPS, zp)
    ) / (2 * D_EPS);
}

double theta_partial_zp_old(
    double x, double xp, double z, double zp
) {
    return (
        theta_old(x, xp, z, zp + D_EPS)
        - theta_old(x, xp, z, zp - D_EPS)
    ) / (2 * D_EPS);
}

double theta_p_old(double x, double z, double xp, double zp, double xpp, double zpp) {
    return xp * theta_partial_x_old(x, xp, z, zp) +
           xpp * theta_partial_xp_old(x, xp, z, zp) +
           zp * theta_partial_z_old(x, xp, z, zp) +
           zpp * theta_partial_zp_old(x, xp, z, zp);
}

// return 0 if t < mx, (t-mx)^2 otherwise
double soft_excess_old(double t, double mx) {
    double exc = fmax(t - mx, 0);
    return exc * exc;
}

double real_speed_norm_old(double x, double xp, double z, double zp) {
    double norm = hypot(x, z);
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
    return hypot(S * sin_theta + whirl_x, S * cos_theta + whirl_z);
}

// a function f(x,x',z,z',u) such that the lagrangian is f(x,x',z,z',θ'(x,x',z,z'))
double lagr_intermediate_old(double x, double xp, double z, double zp, double u) {
    double norm = hypot(x, z);
    double yaw_offset = - M_PI * 250 / (norm + 1000);
    double sin_o = sin(yaw_offset),
          cos_o = cos(yaw_offset),
          fac = -20 * (1/norm - 1/2000.0);
    double whirl_x = fac * (cos_o * x + sin_o * z),
           whirl_z = fac * (-sin_o * x + cos_o * z);
    double a = (whirl_x * zp - whirl_z * xp) / S,
           b = xp * xp + zp * zp;
    double sin_theta = (xp * sqrt(b - a * a) - a * zp) / b;
    double time_integrand = fabs(xp / (S * sin_theta + whirl_x));
    double penalty_term = soft_excess_old(fabs(u) / time_integrand, MAX_YAW_SPEED);
    return time_integrand + PENALTY_FACTOR * penalty_term;
}

double lagrangian_with_precomputed(
    double xp, double x_current, double theta, double theta_p, double *time_int
) {
    double time_integrand = fabs(xp / (S * sin(theta) + x_current));
    *time_int = time_integrand;
    double penalty_term = soft_excess_old(fabs(theta_p) / time_integrand, MAX_YAW_SPEED);
    return time_integrand + PENALTY_FACTOR * penalty_term;
}

double thetap_partial_x_old(
    double x, double xp, double xpp, double z, double zp, double zpp
) {
    return (
        theta_p_old(x + D_EPS, xp, xpp, z, zp, zpp)
        - theta_p_old(x - D_EPS, xp, xpp, z, zp, zpp)
    ) / (2 * D_EPS);
}

double thetap_partial_xp_old(
    double x, double xp, double xpp, double z, double zp, double zpp
) {
    return (
        theta_p_old(x, xp + D_EPS, xpp, z, zp, zpp)
        - theta_p_old(x, xp - D_EPS, xpp, z, zp, zpp)
    ) / (2 * D_EPS);
}

double thetap_partial_z_old(
    double x, double xp, double xpp, double z, double zp, double zpp
) {
    return (
        theta_p_old(x, xp, xpp, z + D_EPS, zp, zpp)
        - theta_p_old(x, xp, xpp, z - D_EPS, zp, zpp)
    ) / (2 * D_EPS);
}

double thetap_partial_zp_old(
    double x, double xp, double xpp, double z, double zp, double zpp
) {
    return (
        theta_p_old(x, xp, xpp, z, zp + D_EPS, zpp)
        - theta_p_old(x, xp, xpp, z, zp - D_EPS, zpp)
    ) / (2 * D_EPS);
}

double lagr_inter_partial_theta_old(double x, double xp, double z, double zp, double u) {
    return (
        lagr_intermediate_old(x, xp, z, zp, u + D_EPS)
        - lagr_intermediate_old(x, xp, z, zp, u - D_EPS)
    ) / (2 * D_EPS);
}

double lagr_partial_x_old(
    double x, double xp, double xpp, double z, double zp, double zpp,
    double theta_p, double partial_u
) {
    double part = (
        lagr_intermediate_old(x + D_EPS, xp, z, zp, theta_p)
        - lagr_intermediate_old(x - D_EPS, xp, z, zp, theta_p)
    ) / (2 * D_EPS);
    return part + thetap_partial_x_old(x, xp, xpp, z, zp, zpp) * partial_u;
}

double lagr_partial_xp_old(
    double x, double xp, double xpp, double z, double zp, double zpp,
    double theta_p, double partial_u
) {
    double part = (
        lagr_intermediate_old(x, xp + D_EPS, z, zp, theta_p)
        - lagr_intermediate_old(x, xp - D_EPS, z, zp, theta_p)
    ) / (2 * D_EPS);
    return part + thetap_partial_xp_old(x, xp, xpp, z, zp, zpp) * partial_u;
}

double lagr_partial_z_old(
    double x, double xp, double xpp, double z, double zp, double zpp,
    double theta_p, double partial_u
) {
    double part = (
        lagr_intermediate_old(x, xp, z + D_EPS, zp, theta_p)
        - lagr_intermediate_old(x, xp, z - D_EPS, zp, theta_p)
    ) / (2 * D_EPS);
    return part + thetap_partial_z_old(x, xp, xpp, z, zp, zpp) * partial_u;
}

double lagr_partial_zp_old(
    double x, double xp, double xpp, double z, double zp, double zpp,
    double theta_p, double partial_u
) {
    double part = (
        lagr_intermediate_old(x, xp, z, zp + D_EPS, theta_p)
        - lagr_intermediate_old(x, xp, z, zp - D_EPS, theta_p)
    ) / (2 * D_EPS);
    return part + thetap_partial_zp_old(x, xp, xpp, z, zp, zpp) * partial_u;
}