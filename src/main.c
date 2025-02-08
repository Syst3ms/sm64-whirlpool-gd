#define _USE_MATH_DEFINES

#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <memory.h>
#include <mm_malloc.h>
#include <stdio.h>
#include <stdlib.h>

#include "autodiff.h"
#include "math_funcs.h"
#include "parameters.h"
#include "renormalization.h"
#include "util.h"

void init_path(
    struct data *d,
    v2d start, v2d end,
    size_t start_i, size_t end_i
) {
    v2d denom = {end_i - start_i, end_i - start_i};
    v2d step = (end - start) / denom;
    v2d vel = step * POINTS;
    v2d zero = {0.0, 0.0};

    v2d pos = start;
    for (size_t i = start_i; i <= end_i; i++) {
        union point *p = &d->points[i];
        p->pos = pos;
        p->vel = vel;
        p->acc = zero;
        p->time_integrand = time_integrand_alone(pos, vel);
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
    size_t n = h->num_hb;
    for (size_t j = 0; j < n; j++) {
        struct hitbox hb = h->hb[j];
        for (size_t i = 0; i < POINTS; i++) {
            v2d p = d->points[i].pos;
            v2d diff = p - hb.pos;
            // branchless
            double scaling = fmax(hb.radius / fast_hypot_v(diff) - 1.0, 0.0);
            d->points[i].pos = p + diff * scaling;
        }
    }
}

void compute_gradient(struct data *d, struct penalty_data *pdata, v2d *delta) {
    d = __builtin_assume_aligned(d, 16);

    union point *cur_pt = &d->points[1];
    double pfac = pdata->rho / 2.0;
    double cur_part_xp = lagr_partial_xp(cur_pt, pfac);
    double cur_part_zp = lagr_partial_zp(cur_pt, pfac);

    for (size_t i = 0; i < POINTS-3; i++) {
        union point *next_pt = &d->points[i+2];
        double next_part_xp = lagr_partial_xp(next_pt, pfac),
               next_part_zp = lagr_partial_zp(next_pt, pfac);
        delta[i][0] = lagr_partial_x(cur_pt, pfac) - (next_part_xp - cur_part_xp) * POINTS;
        delta[i][1] = lagr_partial_z(cur_pt, pfac) - (next_part_zp - cur_part_zp) * POINTS;
        cur_pt = next_pt;
        cur_part_xp = next_part_xp;
        cur_part_zp = next_part_zp;
    }

    union point *next_pt = &d->points[POINTS-1];
    double next_part_xp = lagr_partial_xp(next_pt, pfac),
            next_part_zp = lagr_partial_zp(next_pt, pfac);
    delta[POINTS-3][0] = lagr_partial_x(cur_pt, pfac) - (next_part_xp - cur_part_xp) * POINTS;
    delta[POINTS-3][1] = lagr_partial_z(cur_pt, pfac) - (next_part_zp - cur_part_zp) * POINTS;
}

void renormalize(struct data *d, v2d *renorm_w) {
    d = __builtin_assume_aligned(d, 32);
    union point *points = d->points;
    // arclength renormalization, prevents wonky convergence behavior

    double arclength[POINTS-1],
           tot_arclength = 0.0;

    v2d last_pt = points[0].pos;
    for (size_t i = 1; i < POINTS; i++) {
        v2d current_pt = points[i].pos;
        tot_arclength += (
            arclength[i-1] = fast_hypot_v(current_pt - last_pt)
        );
        last_pt = current_pt;
    }

    double arclength_so_far = 0.0f;
    size_t i = 0; // proportion of arclength covered
    double ref = 0.0;

    for (size_t j = 0; j < POINTS-1; j++) {
        v2d w = points[j].pos, wn = points[j+1].pos;
        double scale = 1.0 / arclength[j];
        double along = (ref - arclength_so_far) * scale;
        while (along <= 1.0) {
            // we need to add a point at the next step
            renorm_w[i] = (1 - along) * w + along * wn;
            i++;
            double inc = tot_arclength / (POINTS-1);
            ref += inc;
            along += inc * scale;
        }

        arclength_so_far += arclength[j];
    }

    if (i == POINTS-1) {
        // add last point if necessary
        renorm_w[i] = points[i].pos;
    } else if (i != POINTS) {
        printf("Renorm only added %d points!\n", i);
        exit(1);
    }

    // update p
    for (size_t i = 0; i < POINTS; i++) {
        points[i].pos = renorm_w[i];
    }
}

// finds an (approximate) minimizer for the unconstrained problem
// given by adding a penalty parameter to the lagrangian (factor is d->penalty_fac)
void optimize_unconstrained(
    struct data *d,
    struct penalty_data *pdata,
    struct history *hist,
    struct hitboxes *active_hitboxes,
    double eps,
    size_t max_iters,
    size_t max_iters_without_change
) {
    struct data best = *d;
    struct momentum mom = {};
    v2d renorm_w[POINTS];
    v2d grad[POINTS-2];

    double obj = objective(d, pdata),
          best_obj = obj;
    size_t iters = 0;
    size_t iters_since_last_best = 0;

    while (iters <= max_iters && iters_since_last_best < max_iters_without_change) {
        if (iters++ % MEM_STORE_RATE == 0) {
            //printf("%d: %f (current), %f (best)\n", iters, obj, best_obj);
            store_into_history(d, hist);
        }

        compute_gradient(d, pdata, grad);

        update_and_apply_momentum(d, &mom, grad, eps);

        push_out_of_hitboxes(d, active_hitboxes);

        renormalize(d, renorm_w);

        recompute_dependent(d);

        obj = objective(d, pdata);

        if (obj < best_obj) {
            iters_since_last_best = 0;
            best = *d;
            best_obj = obj;
        } else {
            iters_since_last_best++;
        }
    }

    if (best_obj <= obj) {
        *d = best;
        obj = best_obj;
    }

    printf("Optimized path down to %f after %d iterations.\n", obj, iters);
}

int main(void) {

    struct data d;
    puts("Initializing path");
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

    
    puts("Initializing hitboxes");
    struct hitbox hb[] = {
        {{0, 0}, 26},
        {{1374, 948}, 150},
        {{-1326, -1202}, 150},
        {{774, 23}, 150}
    };

    struct hitboxes *hitboxes = _mm_malloc(sizeof(struct hitboxes) + sizeof(hb), 16);
    if (hitboxes == NULL) {
        puts("Couldn't allocate hitboxes!");
        exit(1);
    }
    size_t len = sizeof(hb) / sizeof(struct hitbox);
    hitboxes->num_hb = len;
    memcpy(hitboxes->hb, hb, sizeof(hb));
    push_out_of_hitboxes(&d, hitboxes);

    puts("Starting optimization loops");
    struct history hist = init_history(1000);

    struct penalty_data pdata = {};
    double rho = 1.0;

    while (rho < 10000.0) {
        pdata.rho = rho;
        printf("Rho: %f\n", pdata.rho);
        optimize_unconstrained(&d, &pdata, &hist, hitboxes, 0.1, INT_MAX, 15000);
        rho *= 1.2;
    }

    puts("Writing to file...");

    FILE *f = fopen("path.txt", "w");
    if (f == NULL) {
        puts("Couldn't open path file!");
        exit(1);
    }

    for (size_t i = 0; i < hist.next; i++) {
        for (size_t j = 0; j < POINTS; j++) {
            v2d pt = hist.pts[i * POINTS + j];
            fprintf(f, "%f,%f\n", pt[0], pt[1]);
        }
        fprintf(f, "\n");
    }

    for (size_t i = 0; i < POINTS; i++) {
        v2d pt = d.points[i].pos;
        fprintf(f, "%f,%f\n", pt[0], pt[1]);
    }

    fclose(f);

    puts("Writing debug data...");

    compute_output_resampled(&d);

    puts("Done");

    free(hitboxes);
    free_history(hist);

    return 0;
}