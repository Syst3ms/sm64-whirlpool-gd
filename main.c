#define _USE_MATH_DEFINES
#define _POSIX_C_SOURCE 200809L

#include <float.h>
#include <limits.h>
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>

#include "parameters.h"
#include "util.h"
#include "math_funcs.h"
#include "renormalization.h"

void output_debug(struct data *d) {
    FILE *f = fopen("debug.txt", "w");

    fprintf(f, "x,z,theta_p,time_int,ratio\n");
    fprintf(f, "%f,%f\n", d->p.x[0], d->p.z[0]);
    for (int i = 0; i < POINTS-1; i++) {
        double tp = d->theta_p[i], ti = d->time_int[i];
        fprintf(f, "%f,%f,%f,%f,%f\n", d->p.x[i+1], d->p.z[i+1], tp, ti, fabs(tp/ti));
    }

    fclose(f);
}

// construct linear path from start to end, leaving a margin of length `margin`
// on either end
struct path initialize_path(
    double start_x, 
    double start_z, 
    double end_x, 
    double end_z,
    double margin
) {
    struct path p;
    double dx = end_x - start_x,
           dz = end_z - start_z;
    double fac = margin / hypot(dx, dz);
    start_x += dx * fac;
    start_z += dz * fac;
    end_x -= dx * fac;
    end_z -= dz * fac;
    lerp(p.x, start_x, end_x, 0, POINTS-1);
    lerp(p.z, start_z, end_z, 0, POINTS-1);
    return p;
}

struct path initialize_path_inter(
    double start_x, 
    double start_z,
    double inter_x,
    double inter_z,
    double end_x, 
    double end_z,
    double margin
) {
    int half = POINTS/2;
    struct path p;
    double dx1 = inter_x - start_x,
           dz1 = inter_z - start_z;
    double fac1 = margin / hypot(dx1, dz1);
    start_x += dx1 * fac1;
    start_z += dz1 * fac1;
    lerp(p.x, start_x, inter_x, 0, half);
    lerp(p.z, start_z, inter_z, 0, half);

    double dx2 = end_x - inter_x,
           dz2 = end_z - inter_z;
    double fac2 = margin / hypot(dx2, dz2);
    end_x -= dx2 * fac2;
    end_z -= dz2 * fac2;
    lerp(p.x, inter_x, end_x, half, POINTS-1);
    lerp(p.z, inter_z, end_z, half, POINTS-1);

    return p;
}

void push_out_of_hitboxes(
    struct path *p,
    struct hitboxes *h
) {
    int n = h->num_hb; 
    for (int j = 0; j < n; j++) {
        struct hitbox hb = h->hb[j];
        for (int i = 0; i < POINTS; i++) {
            double dx = p->x[i] - hb.x,
                   dz = p->z[i] - hb.z;
            double scaling = hb.radius / hypot(dx, dz) - 1;
            if (scaling > 0) {
                p->x[i] += dx * scaling;
                p->z[i] += dz * scaling;
            }
        }
    }
}

void descend_and_renormalize(
    struct data *d,
    struct hitboxes *hb,
    struct momentum *momentum,
    double eps,
    char first
) {
    double x2[POINTS-1], z2[POINTS-1];
    struct path *p = &d->p;

    for (int i = 0; i < POINTS-1; i++) {
        double x = p->x[i+1],
               xp = d->xp[i], 
               xpp = d->xpp[i],
               z = p->z[i+1], 
               zp = d->zp[i],
               zpp = d->zpp[i],
               theta_p = d->theta_p[i],
               partial_theta = d->partial_theta[i];
        x2[i] = -lagr_partial_xp(x, xp, xpp, z, zp, zpp, theta_p, partial_theta);
        z2[i] = -lagr_partial_zp(x, xp, xpp, z, zp, zpp, theta_p, partial_theta);
    }

    double delta_x[POINTS-2], delta_z[POINTS-2];
    derivative(delta_x, x2, POINTS-1);
    derivative(delta_z, z2, POINTS-1);

    for (int i = 0; i < POINTS-2; i++) {
        double x = p->x[i+1],
               xp = d->xp[i],
               xpp = d->xpp[i],
               z = p->z[i+1], 
               zp = d->zp[i],
               zpp = d->zpp[i],
               theta_p = d->theta_p[i],
               partial_theta = d->partial_theta[i];
        delta_x[i] += lagr_partial_x(x, xp, xpp, z, zp, zpp, theta_p, partial_theta);
        delta_z[i] += lagr_partial_z(x, xp, xpp, z, zp, zpp, theta_p, partial_theta);
    }

    update_and_apply_momentum(p, momentum, delta_x, delta_z, eps, first);

    push_out_of_hitboxes(p, hb);

    // arclength renormalization, prevents wonky convergence behavior

    double arclength[POINTS-1],
           tot_arclength = 0.0;

    for (int i = 0; i < POINTS-1; i++) {
        tot_arclength += (arclength[i] = hypot(p->x[i+1] - p->x[i], p->z[i+1] - p->z[i]));
    }

    double renorm_x[POINTS], renorm_z[POINTS];
    double arclength_so_far = 0.0f;
    int i = 0; // proportion of arclength covered

    for (int j = 0; j < POINTS-1; j++) {
        while (arclength_so_far + arclength[j] >= tot_arclength * i / (POINTS-1)) {
            // we need to add a point at the next step
            double along = (tot_arclength * i / (POINTS-1) - arclength_so_far) / arclength[j];
            renorm_x[i] = (1 - along) * p->x[j]  + along * p->x[j+1];
            renorm_z[i] = (1 - along) * p->z[j]  + along * p->z[j+1];
            i++;
        }

        arclength_so_far += arclength[j];
    }

    if (i == POINTS-1) {
        // add last point if necessary
        renorm_x[i] = p->x[i];
        renorm_z[i] = p->z[i];
    } else if (i != POINTS) {
        printf("Renorm only added %d points!\n", i);
        exit(1);
    }

    // update p
    array_copy(p->x, renorm_x, POINTS);
    array_copy(p->z, renorm_z, POINTS);
    recompute_dependent(d);
}

void optimize(
    struct data *d,
    struct hitboxes *active_hitboxes,
    double min_diff,
    int max_iters
) {
    struct path *p = &d->p;
    struct data best;
    make_copy(&best, d);

    struct memory mem = init_memory(100);

    double obj = objective(d),
          best_obj = obj,
          diff = HUGE_VAL,
          eps = 1.0;
    struct momentum momentum = {{0.0}, {0.0}, {0.0}, {0.0}};
    int iters = 0;
    while (iters <= max_iters && diff > min_diff) {
        if (iters % 100 == 0) {
            printf("%d: %f (current), %f (best)\n", iters, obj, best_obj);
            store_into_memory(&d->p, &mem);
        } 
        iters++;

        descend_and_renormalize(d, active_hitboxes, &momentum, eps, iters==1);

        obj = objective(d);
        diff = best_obj - obj;

        if (diff > 0) {
            make_copy(&best, d);
            best_obj = obj;
        } else {
            diff = HUGE_VAL;
        }
    }

    if (best_obj < obj) {
        make_copy(d, &best);
        obj = best_obj;
    }

    if (iters > max_iters) {
        printf("Max iters reached: %d\n", iters);
    } else {
        printf("Min diff attained: %e\n", diff);
    }

    printf(
        "Optimized path down to %f after %d iterations\nPrinting path to file.\n",
        obj,
        iters
    );

    FILE *f = fopen("path.txt", "w");

    for (int j = 0; j < mem.next; j++) {
        for (int i = 0; i < POINTS; i++) {
            fprintf(f, "%f,%f\n", mem.x[j * POINTS + i], mem.z[j * POINTS + i]);
        }
        fprintf(f, "\n");
    }

    free_mem(&mem);

    for (int i = 0; i < POINTS; i++) {
        fprintf(f, "%f,%f\n", p->x[i], p->z[i]);
    }

    fclose(f);
}

int main(void) {
    /* Chests:
     *  0: -1326/1198
     *  1: 1374/948
     *  2: -1326/-1202
     *  3: 774/23
     * Star: 1274/-1502
     */
    // Direct/between whirlpool and chest 3
    // struct path p = initialize_path(1374, 948, -1326, -1202, 150);
    // Right of chest 3 
    struct path p = initialize_path_inter(1374, 948, 1000, -1000, -1326, -1202, 150);
    // Left of whirlpool
    // struct path p = initialize_path_inter(1374, 948, -400, 400, -1326, -1202, 150);

    struct hitbox hb[] = {
        {0, 0, 26},
        {1374, 948, 150},
        {-1326, -1202, 150},
        {774, 23, 150}
    };

    struct hitboxes *hitboxes = malloc(sizeof(struct hitboxes) + sizeof(hb));

    if (hitboxes == NULL) {
        puts("Couldn't allocate hitboxes!");
        exit(1);
    }

    int len = sizeof(hb) / sizeof(struct hitbox);
    hitboxes->num_hb = len;
    memcpy(hitboxes->hb, hb, sizeof(hb));

    push_out_of_hitboxes(&p, hitboxes);

    struct data d;
    d.p = p;

    recompute_dependent(&d);

    optimize(&d, hitboxes, 1e-9f, INT_MAX);

    puts("Writing debug data...");

    compute_output_resampled(&d);
    output_debug(&d);

    puts("Done");

    free(hitboxes);

    return 0;
}