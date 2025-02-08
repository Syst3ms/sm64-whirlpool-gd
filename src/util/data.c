#include <immintrin.h>

#include "./data.h"

#include "./math.h"
#include "../lagrangian.h"

void update_and_apply_momentum(
    struct data *d,
    struct momentum *mom,
    v2d *delta,
    double eps
) {
    for (int i = 0; i < POINTS-2; i++) {
        v2d del = delta[i];
        struct mom_point mp = mom->points[i];
        v2d m = mom->points[i].xz = BETA_1 * mp.xz + (1 - BETA_1) * del;
        v2d u = mom->points[i].ut = _mm_max_pd(mp.ut * BETA_2, vabs(del));
        d->points[i+1].pos -= eps * m / u;
    }
}

double objective(struct data *d, struct penalty_data *pdata) {
    double sum = 0.0;
    for (size_t i = 1; i < POINTS; i++) {
        sum += compute_lagrangian(&d->points[i], pdata->rho / 2.0);
    }
    return sum / (POINTS-1);
}

void recompute_dependent(struct data *d) {
    for (size_t i = 1; i < POINTS; i++) {
        d->points[i].vel = (d->points[i].pos - d->points[i-1].pos) * POINTS;
    }
    d->points[0].vel = d->points[1].vel;

    for (size_t i = 1; i < POINTS-1; i++) {
        d->points[i].acc = (d->points[i+1].vel - d->points[i].vel) * POINTS;
    }
    
    // zero-pad, maybe same-pad is better?
    v2d zero = {0.0, 0.0};
    d->points[POINTS-1].acc = zero;
}

