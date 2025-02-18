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

void recompute_dependent(struct data *d) {
    v2d prev_pos = d->points[0].pos,
        cur_pos = d->points[1].pos,
        next_pos = d->points[2].pos;

    d->points[0].vel = (4.0 * cur_pos - 3.0 * prev_pos - next_pos) * POINTS / 2.0;
    d->points[0].acc = (prev_pos + next_pos - 2.0 * cur_pos) * POINTS * POINTS;

    for (size_t i = 1; i < POINTS-2; i++) {
        d->points[i].vel = (next_pos - prev_pos) * POINTS / 2.0;
        d->points[i].acc = (prev_pos + next_pos - 2.0 * cur_pos) * POINTS * POINTS;

        prev_pos = cur_pos;
        cur_pos = next_pos;
        next_pos = d->points[i+2].pos;
    }

    d->points[POINTS-2].vel = (next_pos - prev_pos) * POINTS / 2.0;
    d->points[POINTS-1].vel = (prev_pos + 3.0 * next_pos - 4.0 * cur_pos) * POINTS / 2.0;
    d->points[POINTS-1].acc = d->points[POINTS-2].acc = (prev_pos + next_pos - 2.0 * cur_pos) * POINTS * POINTS;
}
