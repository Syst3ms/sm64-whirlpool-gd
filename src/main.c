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
#include "util.h"
#include "renormalization.h"

void debug_print_path(FILE *f, struct data *d) {
    // void* since the offsets vary
    for (size_t i = 0; i < POINTS; i++) {
        struct pt_vel p = d->points[i].pv;
        fprintf(f, "%f,%f,%f\n", p.x, p.z, p.theta);
    }
    fprintf(f, "\n");

    fflush(f);
}

void print_path(FILE *f, void * path, size_t elt_size, size_t count) {
    // void* since the offsets vary
    for (size_t i = 0; i < count; i++) {
        struct pt *p = path + i * elt_size;
        fprintf(f, "%f,%f\n", p->x, p->z);
    }
    fprintf(f, "\n");

    fflush(f);
}

// initialize the path while respecting yaw velocity constraint
void init_path_smooth_yaw(
    struct data *d,
    v2d start, v2d end,
    int start_i, int end_i
) {
    v2d prev_pos;
    double prev_theta;

    {
        prev_pos = d->points[start_i].v.pos = start;
        d->points[start_i].v.vel = (end - start) / (double) (end_i - start_i) * POINTS;
        compute_lagrangian_with_intermediates(&d->points[start_i].pv);
        prev_theta = d->points[start_i].pv.theta;
    }

    for (int i = start_i+1; i < end_i; i++) {
        // printf("%d\n",i);
        v2d step = (end - prev_pos) / (double) (end_i - i + 1);
        // try to go straight toward `end`
        v2d pos = prev_pos + step;
        // second-order backward difference
        v2d vel = step * POINTS;

        double cur_theta;
        double max_range = OVERSMOOTH_FACTOR * theta_maxrange(pos, vel, &cur_theta);
        double angle_diff = remainder(prev_theta - cur_theta, 2*M_PI);

        if (fabs(angle_diff) < max_range) {
            goto found_theta;
        }

        /*
         * max yaw speed exceeded, now try to rotate `pos` around `prev_pos`
         * as little as possible to satisfy the angle constraint
         */
        double r = length_v((end - prev_pos) / (double) (end_i - i + 1));
        // we want to target the side of the admissible range closest to our initial value
        // this should only be used for its sign
        double target_side = -angle_diff;

        // false position method, find phi_a and phi_b that give angles on
        // either side of the target
        double phi = atan2((end - prev_pos)[0], (end - prev_pos)[1]),
               diff = remainder(prev_theta - cur_theta + copysign(max_range, target_side), 2 * M_PI);
        double phi_a, diff_a;
        double phi_b = phi, diff_b = diff;
        double angle_step = 2*diff;

        int j = 0;
        do {
            phi_a = phi_b;
            diff_a = diff_b;
            phi_b += angle_step;

            v2d off_b = {sin(phi_b), cos(phi_b)};
            off_b *= r;
            v2d pos_b = prev_pos + off_b;
            v2d vel_b = off_b * POINTS;

            double theta_b;
            double maxrange_b = OVERSMOOTH_FACTOR * theta_maxrange(pos_b, vel_b, &theta_b);

            diff_b = remainder(prev_theta - theta_b + copysign(maxrange_b, target_side), 2*M_PI);
            j++;

            if (fabs(diff_b) < ONE_HAU) {
                pos = pos_b;
                vel = vel_b;
                cur_theta = theta_b;
                goto found_theta;
            }
        } while (is_sign_same(diff_a, diff_b));

        // false position method, run actual search
        do {
            phi = (phi_a * diff_b - phi_b * diff_a) / (diff_b - diff_a);
            v2d off = {sin(phi), cos(phi)};
            off *= r;

            pos = prev_pos + off;
            vel = off * POINTS;

            max_range = OVERSMOOTH_FACTOR * theta_maxrange(pos, vel, &cur_theta);

            diff = remainder(prev_theta - cur_theta + copysign(max_range, target_side), 2*M_PI);

            if (is_sign_same(diff, diff_a)) {
                phi_a = phi;
                diff_a = diff;
            } else {
                phi_b = phi;
                diff_b = diff;
            }
        } while (fabs(diff) >= ONE_HAU);

        found_theta:
        d->points[i].v.pos = pos;
        d->points[i].v.vel = vel;
        compute_lagrangian_with_intermediates(&d->points[i].pv);

        prev_theta = cur_theta;
        prev_pos = pos;
    }

    d->points[end_i].v.pos = end;
    d->points[end_i].v.vel = (end - prev_pos) * POINTS;
    compute_lagrangian_with_intermediates(&d->points[end_i].pv);
}

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
    double fac = margin / length(dx, dz);
    v2d start = {start_x + dx * fac, start_z + dz * fac};
    v2d end = {end_x - dx * fac, end_z - dz * fac};
    init_path_smooth_yaw(d, start, end, 0, POINTS-1);
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
    double fac1 = margin / length(dx1, dz1);
    v2d start = {start_x + dx1 * fac1, start_z + dz1 * fac1};
    v2d mid = {inter_x, inter_z};
    double dx2 = end_x - inter_x,
           dz2 = end_z - inter_z;
    double fac2 = margin / length(dx2, dz2);
    v2d end = {end_x - dx2 * fac2, end_z - dz2 * fac2};

    init_path(d, start, mid, 0, half);
    init_path(d, mid, end, half, POINTS-1);
}

double objective(struct data *d) {
    double sum = 0.0;
    for (int i = 1; i < POINTS; i++) {
        compute_lagrangian_with_intermediates(&d->points[i].pv);
        sum += d->points[i].pv.time_integrand;
    }
    return sum / (POINTS-1);
}

void descend(struct data *d, v2d *delta) {
    d = __builtin_assume_aligned(d, 16);

    struct pt_vel *cur_pt = &d->points[1].pv;
    double cur_part_xp = lagr_partial_xp(cur_pt);
    double cur_part_zp = lagr_partial_zp(cur_pt);

    for (int i = 1; i < POINTS-1; i++) {
        struct pt_vel *next_pt = &d->points[i+1].pv;
        double next_part_xp = lagr_partial_xp(next_pt),
               next_part_zp = lagr_partial_zp(next_pt);
        delta[i-1][0] = lagr_partial_x(cur_pt) - (next_part_xp - cur_part_xp) * POINTS;
        delta[i-1][1] = lagr_partial_z(cur_pt) - (next_part_zp - cur_part_zp) * POINTS;
        cur_pt = next_pt;
        cur_part_xp = next_part_xp;
        cur_part_zp = next_part_zp;
    }
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
            double scaling = fmax(hb.radius / length_v(diff) - 1.0, 0.0);
            d->points[i].v.pos = p + diff * scaling;
        }
    }
}

// tweak the path so that the yaw velocity doesn't exceed its maximum
void smooth_out(struct data *d) {
    v2d end = d->points[POINTS-1].v.pos;
    v2d prev_pos = d->points[0].v.pos;
    double prev_theta = d->points[0].pv.theta;
    v2d pos = d->points[1].v.pos;
    v2d vel = d->points[1].v.vel;

    for (int i = 1; i < POINTS-1; i++) {
        printf("%d\n", i);

        double cur_theta;
        double max_range = OVERSMOOTH_FACTOR * theta_maxrange(pos, vel, &cur_theta);
        double angle_diff = remainder_2pi(prev_theta - cur_theta);

        if (fabs(angle_diff) < max_range) {
            goto found_theta;
        }

        /*
         * max yaw speed exceeded, now try to move 'pos' along a perpendicular
         * bisector between 'prev_pos' and 'next_pos'
         */
        v2d towards_end = (end - prev_pos) / (double) (POINTS - i);
        v2d normal = ortho(towards_end);

        double target_side;
        v2d pos_a, pos_b = prev_pos + towards_end;
        double diff_a, diff_b;

        {
            v2d vel_b = towards_end * POINTS;

            double theta_b;
            double maxrange_b = OVERSMOOTH_FACTOR * theta_maxrange(pos_b, vel_b, &theta_b);

            double tmpdiff = remainder_2pi(prev_theta - theta_b);
            diff_b = tmpdiff - copysign(maxrange_b, tmpdiff);
            target_side = -tmpdiff;
            
            if (fabs(diff_b) < ONE_HAU) {
                pos = pos_b;
                cur_theta = theta_b;
                goto found_theta;
            }
        }

        // we want to target the side of the admissible range closest to our initial value
        // this should only be used for its sign
        // 'normal' is oriented such that this points the right way
        v2d bracket_step = normal * diff_b;

        do {
            pos_a = pos_b;
            diff_a = diff_b;
            pos_b += bracket_step;

            v2d vel_b = (pos_b - prev_pos) * POINTS;

            double theta_b;
            double maxrange_b = OVERSMOOTH_FACTOR * theta_maxrange(pos_b, vel_b, &theta_b);

            diff_b = remainder_2pi(prev_theta - theta_b + copysign(maxrange_b, target_side));

            if (fabs(diff_b) < ONE_HAU) {
                pos = pos_b;
                cur_theta = theta_b;
                goto found_theta;
            }
        } while (is_sign_same(diff_a, diff_b));

        double diff;

        // false position method, run actual search
        do {
            pos = (pos_a * diff_b - pos_b * diff_a) / (diff_b - diff_a);
            vel = (pos - prev_pos) * POINTS;

            max_range = OVERSMOOTH_FACTOR * theta_maxrange(pos, vel, &cur_theta);

            diff = remainder(prev_theta - cur_theta + copysign(max_range, target_side), 2*M_PI);

            if (is_sign_same(diff, diff_a)) {
                pos_a = pos;
                diff_a = diff;
            } else {
                pos_b = pos;
                diff_b = diff;
            }
        } while (fabs(diff) >= ONE_HAU);

        found_theta:
        d->points[i].v.pos = pos;
        // don't adjust anything else as we don't need it for the rest of the iteration
        // and it will be recomputed in due time from the positions

        prev_pos = pos;
        prev_theta = cur_theta;
        pos = d->points[i+1].v.pos;
        vel = (pos - prev_pos) * POINTS;
    }
}

void renormalize(struct data *d, v2d *renorm_w) {
    d = __builtin_assume_aligned(d, 16);
    union point *points = d->points;
    // arclength renormalization, prevents wonky convergence behavior

    double arclength[POINTS-1],
           tot_arclength = 0.0;

    v2d last_pt = points[0].v.pos;
    for (int i = 1; i < POINTS; i++) {
        v2d current_pt = points[i].v.pos;
        tot_arclength += (
            arclength[i-1] = length_v(current_pt - last_pt)
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
          eps = 1.0;
    int iters = 0;
    int iters_since_last_best = 0;

    while (iters < max_iters && diff > min_diff && iters_since_last_best < max_iters_without_change) {
        printf("Iter: %d\n", iters);
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

        recompute_dependent(d);

        if (iters == max_iters) {
            FILE *f = fopen("debug.txt", "w");
            debug_print_path(f, d);
            puts("Smoothing");

            smooth_out(d);
            recompute_speeds(d);

            debug_print_path(f, d);
            fclose(f);
        } else {
            smooth_out(d); 
        }

        renormalize(d, renorm_w);

        recompute_speeds(d);

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
        print_path(f, mem.pts + j + POINTS, sizeof(struct pt), POINTS);
    }

    free_memory(&mem);

    print_path(f, d->points, sizeof(d->points[0]), POINTS);

    fclose(f);
}

int main(void) {
    puts("Initializing path");
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

    struct hitboxes *hitboxes = _mm_malloc(sizeof(struct hitboxes) + sizeof(hb), 16);

    if (hitboxes == NULL) {
        puts("Couldn't allocate hitboxes!");
        exit(1);
    }

    int len = sizeof(hb) / sizeof(struct hitbox);
    hitboxes->num_hb = len;
    memcpy(hitboxes->hb, hb, sizeof(hb));

    push_out_of_hitboxes(&d, hitboxes);

    puts("Starting descent");

    optimize(&d, hitboxes, 1e-9f, 8, 50000);

    puts("Writing debug data...");

    compute_output_resampled(&d);

    puts("Done");

    _mm_free(hitboxes);

    return 0;
}