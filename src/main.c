#define _USE_MATH_DEFINES

#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>

#include "autodiff.h"
#include "parameters.h"
#include "util.h"
#include "renormalization.h"

void init_path(
    struct data *d,
    v2d start, v2d end,
    int start_i, int end_i
) {
    v2d denom = {end_i - start_i, end_i - start_i};
    v2d step = (end - start) / denom;
    v2d vel = step * POINTS;
    v2d zero = {0.0, 0.0};

    v2d pos = start;
    for (int i = start_i; i <= end_i; i++) {
        struct vec_pt *p = &d->points[i].v;
        p->pos = pos;
        p->vel = vel;
        p->acc = zero;
        ((struct pt_vel *) p)->time_integrand = time_integrand_alone(pos, vel);
        pos += step;
    }
}

// construct linear path from start to end, leaving a margin of length `margin`
// on either end
void initialize_path(
    struct data *d,
    double start_x, 
    double start_z, 
    double end_x, 
    double end_z,
    double margin
) {
    double dx = end_x - start_x,
           dz = end_z - start_z;
    double fac = margin / fast_hypot(dx, dz);
    v2d start = {start_x + dx * fac, start_z + dz * fac};
    v2d end = {end_x - dx * fac, end_z - dz * fac};
    init_path(d, start, end, 0, POINTS-1);
}

void initialize_path_inter(
    struct data *d,
    double start_x, double start_z,
    double inter_x, double inter_z,
    double end_x, double end_z,
    double margin
) {
    int half = POINTS/2;
    double dx1 = inter_x - start_x,
           dz1 = inter_z - start_z;
    double fac1 = margin / fast_hypot(dx1, dz1);
    v2d start = {start_x + dx1 * fac1, start_z + dz1 * fac1};
    v2d mid = {inter_x, inter_z};
    double dx2 = end_x - inter_x,
           dz2 = end_z - inter_z;
    double fac2 = margin / fast_hypot(dx2, dz2);
    v2d end = {end_x - dx2 * fac2, end_z - dz2 * fac2};

    init_path(d, start, mid, 0, half);
    init_path(d, mid, end, half, POINTS-1);
}

void push_out_of_hitboxes(
    struct data *d,
    struct hitboxes *h
) {
    int n = h->num_hb; 
    for (int j = 0; j < n; j++) {
        struct hitbox hb = h->hb[j];
        for (int i = 0; i < POINTS; i++) {
            v2d p = d->points[i].v.pos;
            v2d diff = p - hb.pos;
            // branchless
            double scaling = fmax(hb.radius / fast_hypot_v(diff) - 1.0, 0.0);
            v2d sb = {scaling, scaling};
            d->points[i].v.pos = p + diff * sb;
        }
    }
}

void descend(struct data *d, v2d *delta) {
    d = __builtin_assume_aligned(d, 16);
    d->total_lagr_sum = 0.0;

    struct pt_vel *cur_pt = &d->points[1].pv;
    double cur_part_xp = lagr_partial_xp_with_side_eff(cur_pt, &d->total_lagr_sum);
    double cur_part_zp = lagr_partial_zp(cur_pt);

    for (int i = 0; i < POINTS-2; i++) {
        struct pt_vel *next_pt = &d->points[i+2].pv;
        double next_part_xp = lagr_partial_xp_with_side_eff(next_pt, &d->total_lagr_sum),
               next_part_zp = lagr_partial_zp(next_pt);
        delta[i][0] = lagr_partial_x(cur_pt) - (next_part_xp - cur_part_xp) * POINTS;
        delta[i][1] = lagr_partial_z(cur_pt) - (next_part_zp - cur_part_zp) * POINTS;
        cur_pt = next_pt;
        cur_part_xp = next_part_xp;
        cur_part_zp = next_part_zp;
    }
}

void renormalize(struct data *d, v2d *renorm_w) {
    d = __builtin_assume_aligned(d, 32);
    union point *points = d->points;
    // arclength renormalization, prevents wonky convergence behavior

    double arclength[POINTS-1],
           tot_arclength = 0.0;

    v2d last_pt = points[0].v.pos;
    for (int i = 1; i < POINTS; i++) {
        v2d current_pt = points[i].v.pos;
        tot_arclength += (
            arclength[i-1] = fast_hypot_v(current_pt - last_pt)
        );
        last_pt = current_pt;
    }

    double arclength_so_far = 0.0f;
    int i = 0; // proportion of arclength covered
    double ref = 0.0;

    for (int j = 0; j < POINTS-1; j++) {
        v2d w = points[j].v.pos, wn = points[j+1].v.pos;
        double scale = 1.0 / arclength[j];
        double along = (ref - arclength_so_far) * scale;
        while (along <= 1.0) {
            v2d u = {1 - along, 1 - along};
            v2d v = {along, along};
            // we need to add a point at the next step
            renorm_w[i] = u * w + v * wn;
            i++;
            double inc = tot_arclength / (POINTS-1);
            ref += inc;
            along += inc * scale;
        }

        arclength_so_far += arclength[j];
    }

    if (i == POINTS-1) {
        // add last point if necessary
        renorm_w[i] = points[i].v.pos;
    } else if (i != POINTS) {
        printf("Renorm only added %d points!\n", i);
        exit(1);
    }

    // update p
    for (int i = 0; i < POINTS; i++) {
        points[i].v.pos = renorm_w[i];
    }
}

void optimize(
    struct data *d,
    struct hitboxes *active_hitboxes,
    double min_diff,
    int max_iters,
    int max_iters_without_change
) {
    struct data best = *d;
    struct momentum momentum = {};
    struct memory mem = init_memory(BASE_MEM_SIZE);
    v2d renorm_w[POINTS];
    v2d delta[POINTS-2];

    double obj = objective(d),
          best_obj = obj,
          diff = HUGE_VAL,
          eps = 0.1;
    int iters = 0;
    int iters_since_last_best = 0;

    while (iters <= max_iters && diff > min_diff && iters_since_last_best < max_iters_without_change) {
        if (iters % MEM_STORE_RATE == 0) {
            printf("%d: %f (current), %f (best)\n", iters, obj, best_obj);
            store_into_memory(d, &mem);
        }

        descend(d, delta);

        if (iters++ == 0) {
            for (int i = 0; i < POINTS-2; i++) {
                struct mom_point init = {delta[i], delta[i]};
                momentum.points[i] = init;
            }
        }

        update_and_apply_momentum(d, &momentum, delta, eps);

        push_out_of_hitboxes(d, active_hitboxes);

        renormalize(d, renorm_w);

        recompute_dependent(d);

        obj = objective(d);
        diff = best_obj - obj;

        if (diff > 0) {
            iters_since_last_best = 0;
            best = *d;
            best_obj = obj;
        } else {
            iters_since_last_best++;
            diff = HUGE_VAL;
        }
    }

    if (best_obj <= obj) {
        *d = best;
        obj = best_obj;
    }

    if (iters > max_iters) {
        printf("Max iters reached: %d\n", iters);
    } else if (diff <= min_diff) {
        printf("Min diff attained: %e\n", diff);
    } else if (iters_since_last_best > max_iters_without_change) {
        printf("Ran %d iterations without improving the optimum\n", max_iters_without_change);
    }

    printf(
        "Optimized path down to %f after %d iterations\nPrinting path to file.\n",
        obj,
        iters
    );

    FILE *f = fopen(".\\path.txt", "w");
    if (!f) {
        puts("Couldn't open file for path!");
        exit(1);
    }

    for (int j = 0; j < mem.next; j++) {
        for (int i = 0; i < POINTS; i++) {
            struct pt p = mem.pts[j * POINTS + i];
            fprintf(f, "%f,%f\n", p.x, p.z);
        }
        fprintf(f, "\n");
    }

    free_memory(&mem);

    for (int i = 0; i < POINTS; i++) {
        struct pt p = d->points[i].p;
        fprintf(f, "%f,%f\n", p.x, p.z);
    }

    fclose(f);
}

int main(void) {

    struct data d;

    /* Chests:
     *  0: -1326/1198
     *  1: 1374/948
     *  2: -1326/-1202
     *  3: 774/23
     * Star: 1274/-1502
     */
    // Direct/between whirlpool and chest 3
    initialize_path(&d, 1374, 948, -1326, -1202, 150);
    // Right of chest 3 
    // initialize_path_inter(&d, 1374, 948, 1000, -1000, -1326, -1202, 150);
    // Left of whirlpool
    // initialize_path_inter(&d, 1374, 948, -400, 400, -1326, -1202, 150);

    struct hitbox hb[] = {
        {{0, 0}, 26},
        {{1374, 948}, 150},
        {{-1326, -1202}, 150},
        {{774, 23}, 150}
    };

    struct hitboxes *hitboxes = malloc(sizeof(struct hitboxes) + sizeof(hb));

    if (hitboxes == NULL) {
        puts("Couldn't allocate hitboxes!");
        exit(1);
    }

    int len = sizeof(hb) / sizeof(struct hitbox);
    hitboxes->num_hb = len;
    memcpy(hitboxes->hb, hb, sizeof(hb));

    push_out_of_hitboxes(&d, hitboxes);

    optimize(&d, hitboxes, 1e-9f, INT_MAX, 50000);

    puts("Writing debug data...");

    compute_output_resampled(&d);

    puts("Done");

    free(hitboxes);

    return 0;
}