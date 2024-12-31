#define _USE_MATH_DEFINES
#define _POSIX_C_SOURCE 200809L

#include <float.h>
#include <limits.h>
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>

#include "util.h"
#include "math_funcs.h"

#define S 28.0
#define EPS_INCREASE_INTERVAL 2000

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
    double eps
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

    double max_delta_norm = 0.0f;

    for (int i = 0; i < POINTS-2; i++) {
        double x = p->x[i+1],
               xp = d->xp[i],
               xpp = d->xpp[i],
               z = p->z[i+1], 
               zp = d->zp[i],
               zpp = d->zpp[i],
               theta_p = d->theta_p[i],
               partial_theta = d->partial_theta[i];
        max_delta_norm = fmax(
            max_delta_norm,
            hypot(
                delta_x[i] += lagr_partial_x(x, xp, xpp, z, zp, zpp, theta_p, partial_theta),
                delta_z[i] += lagr_partial_z(x, xp, xpp, z, zp, zpp, theta_p, partial_theta)
            )
        );
    }

    for (int i = 1; i < POINTS-1; i++) {
        p->x[i] += MOMENTUM_PARAM * momentum->x[i-1] - delta_x[i-1] / max_delta_norm * eps;
        p->z[i] += MOMENTUM_PARAM * momentum->z[i-1] - delta_z[i-1] / max_delta_norm * eps; 
    }

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
    double threshold,
    double min_diff,
    int max_iters
) {
    struct path *p = &d->p;
    struct data backup;
    make_copy(&backup, d);

    struct memory mem = init_memory(100);

    double obj = objective(d),
          prev_obj = obj,
          diff = HUGE_VAL,
          eps = 1.0;
    struct momentum momentum = {{0.0}, {0.0}};
    int iters = 0;
    int streak = 0;
    while (eps >= threshold && iters <= max_iters && diff > min_diff) {
        if (iters % 100 == 0) {
            printf("%d: %f\n", iters, obj);
            store_into_memory(&d->p, &mem);
        } 
        iters++;

        descend_and_renormalize(d, active_hitboxes, &momentum, eps);

        obj = objective(d);
        diff = prev_obj - obj;

        if (diff < 0) { // worsening
            streak = 0;
            eps /= 2;
            printf("eps = %f\n", eps);

            for (int i = 1; i < POINTS-1; i++) {
                momentum.x[i-1] = momentum.z[i-1] = 0;
            }

            make_copy(d, &backup);
            diff = HUGE_VALF;
        } else {
            streak++;
            if (streak % EPS_INCREASE_INTERVAL == 0) {
                eps *= 2;
                printf("eps = %f\n", eps);
            }

            for (int i = 1; i < POINTS-1; i++) {
                momentum.x[i-1] = p->x[i] - backup.p.x[i];
                momentum.z[i-1] = p->z[i] - backup.p.z[i];
            }
            make_copy(&backup, d);
            prev_obj = obj;
        }
    }

    if (prev_obj < obj) {
        make_copy(d, &backup);
        obj = prev_obj;
    }

    if (eps < threshold) {
        printf("Eps broke threshold: %e\n", threshold);
    } else if (iters > max_iters) {
        printf("Max iters reached: %d\n", max_iters);
    } else {
        printf("Min diff attained: %e\n", min_diff);
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

    puts("Done");

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

unsigned short radians_to_au(double rad) {
    unsigned short a = (unsigned short) ((rad / M_2_PI) * 0x10000);
    return a & 0xFFFFFFF0;
}

// finds the point between (x0,z0) and (x1,z1) that is target_length away
// from (ref_x, ref_z), and stores it in (res_x,res_z)
void find_resampled(
    double ref_x, double ref_z,
    double x0, double z0, double x1, double z1,
    double target_length,
    double *res_x, double *res_z
) {
    if (x0 == ref_x && z0 == ref_z) {
        // direct lerp
        double fac = target_length / hypot(x1 - x0, z1 - z0);
        *res_x = (1 - fac) * x0 + fac * x1;
        *res_z = (1 - fac) * z0 + fac * z1;
        return; 
    }

    double vx = x1 - x0, vz = z1 - z0, wx = x0 - ref_x, wz = z0 - ref_z;

    double snv = vx * vx + vz * vz;
    double dp = vx * wx + vz * wz;
    double snw = wx * wx + wz * wz;

    if (dp < 0) {
        puts("Path segments went opposite ways, increase sample count");
        exit(1);
    }

    double fac = (-dp + sqrt(dp * dp + snv * (target_length * target_length - snw))) / snv;

    *res_x = (1 - fac) * x0 + fac * x1;
    *res_z = (1 - fac) * z0 + fac * z1;  
}

void output_debug(
    int length,
    double *renorm_x,
    double *renorm_z,
    double *yaws
) {
    FILE *f = fopen("debug.txt", "w");

    unsigned short prev_au;
    for (int i = 0; i < length; i++) {
        unsigned short cur = radians_to_au(yaws[i]);
        if (i == 0) {
            fprintf(f, "%f,%f,%f,%d\n", renorm_x[i], renorm_z[i], yaws[i], cur);
        } else {
            short diff = cur - prev_au;
            fprintf(
                f, "%f,%f,%f,%d,%d%s\n",
                renorm_x[i], renorm_z[i], yaws[i], cur, diff,
                abs(diff) > 640 ? "!" : ""
            );
        }
        prev_au = cur;
    }

    fclose(f);
}

void inc_and_realloc_if_necessary(int *i, int *arr_siz, double **renorm_x, double **renorm_z) {
    (*i)++;

    if (*i == *arr_siz) {
        *arr_siz += POINTS;
        *renorm_x = realloc(*renorm_x, *arr_siz * sizeof(double));
        *renorm_z = realloc(*renorm_z, *arr_siz * sizeof(double));
        if (renorm_x == NULL || renorm_z == NULL) {
            puts("Couldn't reallocate renorm arrays!");
            exit(1);
        }
    }
}

void compute_and_output_debug(struct path *p) {
    double cum_arclength[POINTS],
           tot_arclength = 0.0;

    cum_arclength[0] = 0.0;
    for (int i = 1; i < POINTS; i++) {
        tot_arclength += hypot(p->x[i] - p->x[i-1], p->z[i] - p->z[i-1]);  
        cum_arclength[i] = tot_arclength;
    }

    // resample to have arclength S (i.e 1 frame) between every step

    int arr_siz = POINTS;
    double *renorm_x = malloc(arr_siz * sizeof(double)),
           *renorm_z = malloc(arr_siz * sizeof(double));
    if (renorm_x == NULL || renorm_z == NULL) {
        puts("Couldn't allocate renorm arrays!");
        exit(1);
    }

    renorm_x[0] = p->x[0];
    renorm_z[0] = p->z[0];
    int i = 0;
    for (int j = 0; j < POINTS-1; j++) {
        while (cum_arclength[j+1] >= (i+1) * S) {
            find_resampled(
                renorm_x[i], renorm_z[i],
                p->x[j], p->z[j], p->x[j+1], p->z[j+1],
                S,
                &renorm_x[i+1], &renorm_z[i+1]
            );

            inc_and_realloc_if_necessary(&i, &arr_siz, &renorm_x, &renorm_z);
        }
    }

    if (cum_arclength[POINTS-1] > i * S) {
        // add one last short segment
        inc_and_realloc_if_necessary(&i, &arr_siz, &renorm_x, &renorm_z);

        renorm_x[i] = p->x[POINTS-1];
        renorm_z[i] = p->z[POINTS-1];
    }

    int len_renorm = i+1;

    double *renorm_xp = malloc((len_renorm - 1) * sizeof(double)),
           *renorm_zp = malloc((len_renorm - 1) * sizeof(double));
    derivative(renorm_xp, renorm_x, len_renorm);
    derivative(renorm_zp, renorm_z, len_renorm);

    double *yaws = malloc(len_renorm * sizeof(double));

    // now that we have proper step sizes, we can find out our prescribed yaws
    yaws[0] = theta(renorm_x[0], renorm_xp[0], renorm_z[0], renorm_zp[0]);

    for (int i = 1; i < len_renorm; i++) {
        yaws[i] = theta(renorm_x[i], renorm_xp[i-1], renorm_z[i], renorm_zp[i-1]);
    }
 
    output_debug(len_renorm, renorm_x, renorm_z, yaws);

    free(yaws);
    free(renorm_x);
    free(renorm_z);
    free(renorm_xp);
    free(renorm_zp);
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

    optimize(&d, hitboxes, 1e-5f, 1e-9f, INT_MAX);

    compute_and_output_debug(&d.p);

    free(hitboxes);

    return 0;
}