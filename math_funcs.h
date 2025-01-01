#pragma once

#define LAGR_PARTIAL(v) double lagr_partial_##v(\
        double x, double xp, double xpp, double z, double zp, double zpp,\
        double theta, double partial_theta\
    )

LAGR_PARTIAL(x);
LAGR_PARTIAL(xp);
LAGR_PARTIAL(z);
LAGR_PARTIAL(zp);
double lagr_inter_partial_theta(double x, double xp, double z, double zp, double u);
double lagrangian_with_precomputed(
    double xp, double x_current, double theta, double theta_p, double *penalty
);

void compute_theta_and_current(double x, double xp, double z, double zp, double *theta, double *x_current);
double theta(double x, double xp, double z, double zp);
double theta_p(double x, double xp, double xpp, double z, double zp, double zpp);
double real_speed_norm(double x, double xp, double z, double zp);
double objective(struct data *d);

void recompute_dependent(struct data *d);