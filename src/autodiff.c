#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

#include "autodiff.h"
#include "old_compute.h"
#include "parameters.h"
#include "types.h"
#include "util.h"

struct autodiff {
    double dx, dz, dxp, dzp;
    double v;
} __attribute__((aligned(32)));

union ad_u {
    struct autodiff ad;
    struct {
        v4d part;
        double v;
    };
} __attribute__((aligned(32)));

struct autodiff hypot_xz(double x, double z) {
    union ad_u ret_u;
    double norm = sqrt(x * x + z * z);
    double inv_norm = 1 / norm;
    v4d val = {x, z, 0.0, 0.0};
    ret_u.part = val * inv_norm;
    ret_u.ad.v = norm;
    return ret_u.ad;
}

struct autodiff hypot_xpzp(double xp, double zp) {
    union ad_u ret_u;
    double norm = sqrt(xp * xp + zp * zp);
    double inv_norm = 1 / norm;
    v4d val = {0.0, 0.0, xp, zp};
    ret_u.part = val * inv_norm;
    ret_u.ad.v = norm;
    return ret_u.ad;
}

struct autodiff add_ad(struct autodiff a, struct autodiff b) {
    union ad_u *a_u = (union ad_u *) &a,
               *b_u = (union ad_u *) &b;
    union ad_u ret_u;
    ret_u.part = a_u->part + b_u->part;
    ret_u.ad.v = a.v + b.v;
    return ret_u.ad;
}

struct autodiff sub_ad(struct autodiff a, struct autodiff b) {
    union ad_u *a_u = (union ad_u *) &a,
               *b_u = (union ad_u *) &b;
    union ad_u ret_u;
    ret_u.part = a_u->part - b_u->part;
    ret_u.ad.v = a.v - b.v;
    return ret_u.ad;
}

struct autodiff add_scalar_ad(struct autodiff a, double s) {
    struct autodiff ret = a;
    ret.v += s;
    return ret;
}

struct autodiff mul_ad(struct autodiff a, struct autodiff b) {
    union ad_u *a_u = (union ad_u *) &a,
               *b_u = (union ad_u *) &b;
    union ad_u ret_u;
    ret_u.part = a.v * b_u->part + b.v * a_u->part;
    ret_u.ad.v = a.v * b.v;
    return ret_u.ad;
}

struct autodiff mul_scalar_ad(struct autodiff a, double s) {
    union ad_u *a_u = (union ad_u *) &a;
    union ad_u ret_u;
    ret_u.part = a_u->part * s;
    ret_u.ad.v = a.v * s;
    return ret_u.ad;
}

struct autodiff scalar_div_ad(double s, struct autodiff a) {
    union ad_u *a_u = (union ad_u *) &a;
    union ad_u ret_u;
    double fac = - s / (a.v * a.v);
    ret_u.part = a_u->part * fac;
    ret_u.ad.v = s / a.v;
    return ret_u.ad;
}

struct autodiff div_ad(struct autodiff a, struct autodiff b) {
    union ad_u ret_u;
    double den = 1 / b.v;
    double fac = a.v * den * den;
    union ad_u *a_u = (union ad_u *) &a,
               *b_u = (union ad_u *) &b;

    ret_u.ad.v = a.v * den;
    ret_u.part = a_u->part * den - b_u->part * fac;
    
    return ret_u.ad;
}

void sincos_ad(struct autodiff a, struct autodiff *sin_res, struct autodiff *cos_res) {
    double sin_a = sin(a.v), cos_a = cos(a.v);
    union ad_u *a_u = (union ad_u *) &a,
               *sin_u = (union ad_u *) sin_res,
               *cos_u = (union ad_u *) cos_res;

    sin_u->part = cos_a * a_u->part;
    sin_res->v = sin_a;

    cos_u->part = -sin_a * a_u->part;
    cos_res->v = cos_a;
}

struct autodiff dot_xz_ad(struct autodiff a, struct autodiff b, double x, double z) {
    union ad_u *a_u = (union ad_u *) &a, 
               *b_u = (union ad_u *) &b;
    union ad_u ret_u;

    v4d extra = {a.v, b.v, 0.0, 0.0};

    ret_u.part = a_u->part * x + b_u->part * z + extra;
    ret_u.ad.v = a.v * x + b.v * z;

    return ret_u.ad;
}

struct autodiff dot_xpzp_ad(struct autodiff a, struct autodiff b, double xp, double zp) {
    union ad_u *a_u = (union ad_u *) &a,
               *b_u = (union ad_u *) &b;

    union ad_u ret_u;

    v4d extra = {0.0, 0.0, a.v, b.v};

    ret_u.part = a_u->part * xp + b_u->part * zp + extra;
    ret_u.ad.v = a.v * xp + b.v * zp;

    return ret_u.ad;
}

struct autodiff dot_nx_z_ad(struct autodiff a, struct autodiff b, double x, double z) {
    union ad_u *a_u = (union ad_u *) &a,
               *b_u = (union ad_u *) &b;

    union ad_u ret_u;

    v4d extra = {-a.v, b.v, 0.0, 0.0};

    ret_u.part = - a_u->part * x + b_u->part * z + extra;
    ret_u.ad.v = - a.v * x + b.v * z;

    return ret_u.ad;
}

struct autodiff dot_nxp_zp_ad(struct autodiff a, struct autodiff b, double xp, double zp) {
    union ad_u *a_u = (union ad_u *) &a,
               *b_u = (union ad_u *) &b;

    union ad_u ret_u;

    v4d extra = {0.0, 0.0, -a.v, b.v};

    ret_u.part = - a_u->part * xp + b_u->part * zp + extra;
    ret_u.ad.v = - a.v * xp + b.v * zp;

    return ret_u.ad;
}

struct autodiff acos_ad(struct autodiff a) {
    union ad_u *a_u = (union ad_u *) &a;
    union ad_u ret_u;
    double fac = - 1 / sqrt(1 - a.v * a.v);
    ret_u.part = a_u->part * fac;
    ret_u.ad.v = acos(a.v);
    return ret_u.ad;
}

struct autodiff atan2_zpxp_ad(double xp, double zp) {
    union ad_u ret_u;
    double fac = 1 / (zp * zp + xp * xp);
    v4d v = {0.0, 0.0, -zp, xp};
    ret_u.part = v * fac;
    ret_u.ad.v = atan2(zp, xp);
    return ret_u.ad;
}

double soft_excess(double t, double mx) {
    double exc = fmax(t - mx, 0);
    return exc * exc;
}

double compute_lagrangian_(struct pt_vel *pt) {
    pt = __builtin_assume_aligned(pt, 16);
    double x = pt->x, z = pt->z, xp = pt->xp, zp = pt->zp, xpp = pt->xpp, zpp = pt->zpp;
    struct autodiff norm = hypot_xz(x, z), speed_norm = hypot_xpzp(xp, zp);
    double norm_d = length(x, z), speed_norm_d = length(xp, zp);
    struct autodiff yaw_offset = scalar_div_ad(- M_PI * 250.0, add_scalar_ad(norm, 1000.0));
    double yaw_offset_d = - M_PI * 250.0 / (norm_d + 1000.0);
    struct autodiff sin_o, cos_o;
    sincos_ad(yaw_offset, &sin_o, &cos_o);
    double sin_o_d = sin(yaw_offset_d), cos_o_d = cos(yaw_offset_d);
    struct autodiff fac = add_scalar_ad(scalar_div_ad(-20.0, norm), 0.01);
    double fac_d = -20.0 / norm_d + 0.01;
    struct autodiff whirl_x_no_fac = dot_xz_ad(cos_o, sin_o, x, z),
                    whirl_z_no_fac = dot_nx_z_ad(sin_o, cos_o, x, z);
    double whirl_x_no_fac_d = cos_o_d * x + sin_o_d * z,
           whirl_z_no_fac_d = - sin_o_d * x + cos_o_d * z;
    struct autodiff det = mul_ad(
        dot_nxp_zp_ad(whirl_z_no_fac, whirl_x_no_fac, xp, zp),
        fac
    );
    double det_d = fac_d * (-whirl_z_no_fac_d * xp + whirl_x_no_fac_d * zp);
    struct autodiff theta = sub_ad(
        acos_ad(div_ad(det, mul_scalar_ad(speed_norm, S))),
        atan2_zpxp_ad(xp, zp)
    );
    double theta_d = acos(det_d / (S * speed_norm_d)) - atan2(zp, xp);
    double theta_p = xp * theta.dx + zp * theta.dz + xpp * theta.dxp + zpp * theta.dzp;
    double theta_p_d = theta_p_old(x, z, xp, zp, xpp, zpp);
    double time_integrand = fabs(xp / (S * sin(theta.v) + whirl_x_no_fac.v * fac.v));
    double time_integrand_d = fabs(xp / (S * sin(theta_d) + whirl_x_no_fac_d * fac_d));
    double penalty_term = soft_excess(fabs(theta_p) / time_integrand, MAX_YAW_SPEED);
    double penalty_term_d = soft_excess(fabs(theta_p_d) / time_integrand_d, MAX_YAW_SPEED);
    double lagrangian = time_integrand + PENALTY_FACTOR * penalty_term;
    double lagrangian_d = time_integrand_d + PENALTY_FACTOR * penalty_term_d;
    double lagrangian_ref = lagr_intermediate_old(x, xp, z, zp, theta_p_d);
    (void) lagrangian_d, (void) lagrangian_ref;

    return lagrangian;
}

void compute_lagrangian_with_intermediates(struct pt_vel *pt) {
    double x = pt->x, z = pt->z, xp = pt->xp, zp = pt->zp;
    double norm = length(x, z),
           speed_norm = length(xp, zp);
    double yaw_offset = - M_PI * 250.0 / (norm + 1000.0);
    double sin_o = sin(yaw_offset),
          cos_o = cos(yaw_offset),
          fac = -20 * (1/norm - 1/2000.0);
    double whirl_x = fac * (cos_o * x + sin_o * z),
           whirl_z = fac * (-sin_o * x + cos_o * z);
    double det = whirl_x * zp - whirl_z * xp;
    pt->theta = acos(det / (S * speed_norm)) - atan2(zp, xp);
    pt->time_integrand = fabs(xp / (S * sin(pt->theta) + whirl_x));
}

double compute_lagrangian__(struct pt_vel *pt) {
    pt = __builtin_assume_aligned(pt, 16);
    double x = pt->x, z = pt->z, xp = pt->xp, zp = pt->zp, xpp = pt->xpp, zpp = pt->zpp;
    struct autodiff norm = hypot_xz(x, z), speed_norm = hypot_xpzp(xp, zp);
    struct autodiff yaw_offset = scalar_div_ad(- M_PI * 250.0, add_scalar_ad(norm, 1000.0));
    struct autodiff sin_o, cos_o;
    sincos_ad(yaw_offset, &sin_o, &cos_o);
    struct autodiff fac = add_scalar_ad(scalar_div_ad(-20.0, norm), 0.01);
    struct autodiff whirl_x_no_fac = dot_xz_ad(cos_o, sin_o, x, z),
                    whirl_z_no_fac = dot_nx_z_ad(sin_o, cos_o, x, z);
    struct autodiff det = mul_ad(
        dot_nxp_zp_ad(whirl_z_no_fac, whirl_x_no_fac, xp, zp),
        fac
    );
    struct autodiff theta = sub_ad(
        acos_ad(div_ad(det, mul_scalar_ad(speed_norm, S))),
        atan2_zpxp_ad(xp, zp)
    );
    double theta_p = xp * theta.dx + zp * theta.dz + xpp * theta.dxp + zpp * theta.dzp;
    double time_integrand = fabs(xp / (S * sin(theta.v) + whirl_x_no_fac.v * fac.v));
    double penalty_term = soft_excess(fabs(theta_p) / time_integrand, MAX_YAW_SPEED);
    double lagrangian = time_integrand + PENALTY_FACTOR * penalty_term;

    return lagrangian;
}

double time_integrand_alone(v2d pos, v2d vel) {
    double x = pos[0], z = pos[1], xp = vel[0], zp = vel[1];
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
    return fabs(xp / (S * sin_theta + whirl_x));
}

struct ad1d hypot_ad1(struct ad1d a, struct ad1d b) {
    double hyp = sqrt(a.v * a.v + b.v * b.v);
    struct ad1d ret = {
        hyp,
        (a.v * a.d + b.v * b.d) / hyp
    };
    return ret;
}

struct ad1d add_scalar_ad1(struct ad1d a, double s) {
    struct ad1d ret = a;
    ret.v += s;
    return ret;
}

struct ad1d mul_scalar_ad1(struct ad1d a, double s) {
    struct ad1d ret = a;
    ret.v *= s;
    ret.d *= s;
    return ret;
}

struct ad1d scalar_div_ad1(double s, struct ad1d a) {
    struct ad1d ret;
    ret.v = s / a.v;
    ret.d = - a.d * s / (a.v * a.v);
    return ret;
}

void sincos_ad1(struct ad1d a, struct ad1d *sin_res, struct ad1d *cos_res) {
    double sin_a = sin(a.v), cos_a = cos(a.v);

    sin_res->d = cos_a * a.d;
    sin_res->v = sin_a;

    cos_res->d = -sin_a * a.d;
    cos_res->v = cos_a;
}

struct ad1d div_ad1(struct ad1d a, struct ad1d b) {
    double den = 1 / b.v;
    double fac = a.v * den * den;

    struct ad1d ret;
    ret.v = a.v * den;
    ret.d = a.d * den - b.d * fac;
    
    return ret;
}

struct ad1d add_ad1(struct ad1d a, struct ad1d b) {
    struct ad1d ret;
    ret.v = a.v + b.v;
    ret.d = a.d + b.d;
    return ret;
}

struct ad1d sub_ad1(struct ad1d a, struct ad1d b) {
    struct ad1d ret;
    ret.v = a.v - b.v;
    ret.d = a.d - b.d;
    return ret;
}

struct ad1d mul_ad1(struct ad1d a, struct ad1d b) {
    struct ad1d ret;
    ret.v = a.v * b.v;
    ret.d = a.v * b.d + a.d * b.v;
    return ret;
}

struct ad1d acos_ad1(struct ad1d a) {
    struct ad1d ret;
    double fac = - 1 / sqrt(1 - a.v * a.v);
    ret.v = acos(a.v);
    ret.d = a.v * fac;
    return ret;
}

struct ad1d atan2_ad1(struct ad1d a, struct ad1d b) {
    struct ad1d ret;
    ret.v = atan2(a.v, b.v);
    ret.d = (b.v * a.d - a.v * b.d) / (a.v * a.v + b.v * b.v);
    return ret;
}

void theta_ad1(struct ad1d x, struct ad1d z, struct ad1d xp, struct ad1d zp, struct ad1d *out, double *maxrange_out) {
    struct ad1d norm = hypot_ad1(x, z),
                vel_norm = hypot_ad1(xp, zp);
    struct ad1d fac = add_scalar_ad1(scalar_div_ad1(-20, norm), 0.01);
    struct ad1d yaw_offset = scalar_div_ad1(-M_PI*250, add_scalar_ad1(norm, 1000));
    struct ad1d sin_o, cos_o;
    sincos_ad1(yaw_offset, &sin_o, &cos_o);
    struct ad1d whirl_x_nf = add_ad1(mul_ad1(cos_o, x), mul_ad1(sin_o, z)),
                whirl_z_nf = sub_ad1(mul_ad1(cos_o, z), mul_ad1(sin_o, x));
    struct ad1d det = mul_ad1(fac, sub_ad1(mul_ad1(whirl_x_nf, zp), mul_ad1(whirl_z_nf, xp)));
    struct ad1d theta = sub_ad1(
        acos_ad1(
            div_ad1(det, mul_scalar_ad1(vel_norm, S))
        ),
        atan2_ad1(zp, xp)
    );
    *out = theta;
    *maxrange_out = OVERSMOOTH_FACTOR * MAX_YAW_SPEED / POINTS * fabs(xp.v / (S * sin(theta.v) + whirl_x_nf.v * fac.v));
}