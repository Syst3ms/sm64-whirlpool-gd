#define _USE_MATH_DEFINES
#define _POSIX_C_SOURCE 200809L

#include <math.h>
#include <stdio.h>
#include <memory.h>
#include <stdlib.h>
#include <float.h>
#include <limits.h>
#include <time.h>

#define POINTS 128
#define S 28.0
#define MOMENTUM_PARAM 0.8
#define BASE_MEM_SIZE 100
#define MARIO_HITBOX_SIZE 100

// TODO: SIMD
struct path {
    double x[POINTS], z[POINTS];
    double xp[POINTS-1], zp[POINTS-1];
    double diff_eps[POINTS-1];
};

struct momentum {
    double x[POINTS-2], z[POINTS-2];
};

struct hitbox {
    double x, z, radius;
};

struct hitboxes {
    int num_hb;
    struct hitbox hb[];
};

struct memory {
    int size, next;
    double *x, *z;
};

void array_copy(double *dst, double *src, int length) {
    memcpy(dst, src, length * sizeof(double));
}

void derivative(double *dst, double *src, int srclen) {
    for (int i = 0; i < srclen-1; i++) {
        dst[i] = (src[i+1] - src[i]) * POINTS;
    }
}

// lerp between start and end on [start_i..end_i] (inclusive)
void lerp(double *dst, double start, double end, int start_i, int end_i) {
    double step = (end - start) / (end_i - start_i);
    for (int i = 0; i <= end_i - start_i; i++) {
        dst[start_i + i] = start + step * i;
    }
}

float opt_numdiff_eps(double x, double xp, double z, double zp) {
    return sqrt(DBL_EPSILON) * (fabs(x) + fabs(xp) + fabs(z) + fabs(zp) + sqrt(DBL_EPSILON));
}

void recompute_dependent(struct path *p) {
    derivative(p->xp, p->x, POINTS);
    derivative(p->zp, p->z, POINTS);

    for (int i = 1; i < POINTS; i++) {
        p->diff_eps[i-1] = opt_numdiff_eps(p->x[i], p->xp[i-1], p->z[i], p->zp[i-1]);
    }
}

double prescribed_angle(double x, double z, double xp, double zp) {
    double norm = hypot(x, z);
    double yaw_offset = - M_PI * 250 / (norm + 1000);
    double sin_o = sin(yaw_offset),
        cos_o = cos(yaw_offset),
        fac = -20 * (1/norm - 1/2000.0f);
    double current_x = fac * (cos_o * x + sin_o * z),
        current_z = fac * (-sin_o * x + cos_o * z);
    double a = S * zp,
        b = -S * xp,
        c = current_x * zp - current_z * xp;
    double theta_star = 2 * atan((b + sqrt(a*a + b*b - c*c)) / (a - c));
    return theta_star;
}

// time taken a.k.a the integrand a.k.a the lagrangian
double time_integrand(double x, double xp, double z, double zp) {
    double norm = hypot(x, z);
    double yaw_offset = - M_PI * 250 / (norm + 1000);
    double sin_o = sin(yaw_offset),
          cos_o = cos(yaw_offset),
          fac = -20 * (1/norm - 1/2000.0f);
    double current_x = fac * (cos_o * x + sin_o * z),
          current_z = fac * (-sin_o * x + cos_o * z);
    double a = S * zp,
          b = -S * xp,
          c = current_x * zp - current_z * xp;
    double inner = (b + sqrt(a*a + b*b - c*c)) / (a - c);
    double cos_theta_star = 2 / (1 + inner * inner) - 1,
          sin_theta_star = 2 * inner / (1 + inner * inner);
    double i1 = hypot(xp, zp);
    double i2 = S * cos_theta_star + current_x;
    double i3 = S * sin_theta_star + current_z;
    double i4 = hypot(i2, i3);
    return i1 / i4;
}

double total_time(struct path *p) {
    double total = 0.0f;
    double lengths[POINTS-1];

    for (int i = 1; i < POINTS; i++) {
        total += lengths[i-1] = time_integrand(
            p->x[i], p->xp[i-1], p->z[i], p->zp[i-1]
        );
    }

    return total / POINTS;
}

double partial_x(double x, double xp, double z, double zp, double diff_eps) {
    return (
        time_integrand(x + diff_eps, xp, z, zp)
        - time_integrand(x - diff_eps, xp, z, zp)
    ) / (2 * diff_eps);
}

double partial_xp(double x, double xp, double z, double zp, double diff_eps) {
    return (
        time_integrand(x, xp + diff_eps, z, zp)
        - time_integrand(x, xp - diff_eps, z, zp)
    ) / (2 * diff_eps);
}

double partial_z(double x, double xp, double z, double zp, double diff_eps) {
    return (
        time_integrand(x, xp, z + diff_eps, zp)
        - time_integrand(x, xp, z - diff_eps, zp)
    ) / (2 * diff_eps);
}

double partial_zp(double x, double xp, double z, double zp, double diff_eps) {
    return (
        time_integrand(x, xp, z, zp + diff_eps)
        - time_integrand(x, xp, z, zp - diff_eps)
    ) / (2 * diff_eps);
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

void remove_equal_points(struct path *p) {
    int last_uneq_i = 0;
    double last_uneq_x = p->x[0],
           last_uneq_z = p->z[0];
    for (int i = 1; i < POINTS; i++) {
        double x = p->x[i], z = p->z[i];
        if (x == last_uneq_x && z == last_uneq_z && i < POINTS-1) {
            // do lerp step unconditionally if last index has been reached
            // in case the array ends with a streak of equal values
            continue;
        }

        if (i - last_uneq_i > 1) {
            int len = i - last_uneq_i;
            for (int j = 1; j < len-1; j++) {
                double fac = ((double) j) / len;
                p->x[last_uneq_i + j] = last_uneq_x * fac + x * (1-fac);
                p->z[last_uneq_i + j] = last_uneq_z * fac + z * (1-fac);
            }
        }

        last_uneq_x = x;
        last_uneq_z = z;
        last_uneq_i = i; 
    }
}

void compute_arclength(struct path *p, double *arclengths, double *tot_arclength) {
    for (int i = 0; i < POINTS-1; i++) {
        *tot_arclength += (arclengths[i] = hypot(p->x[i+1] - p->x[i], p->z[i+1] - p->z[i]));
    }
}

void descend_and_renormalize(
    struct path *p,
    struct hitboxes *hb,
    struct momentum *momentum,
    double eps
) {
    double x2[POINTS-1], z2[POINTS-1];

    for (int i = 1; i < POINTS; i++) {
        double x = p->x[i],
               xp = p->xp[i-1], 
               z = p->z[i], 
               zp = p->zp[i-1],
               diff_eps = p->diff_eps[i-1];
        x2[i-1] = -partial_xp(x, xp, z, zp, diff_eps);
        z2[i-1] = -partial_zp(x, xp, z, zp, diff_eps);
    }

    double delta_x[POINTS-2], delta_z[POINTS-2];
    derivative(delta_x, x2, POINTS-1);
    derivative(delta_z, z2, POINTS-1);

    double max_delta_norm = 0.0f;

    for (int i = 1; i < POINTS-1; i++) {
        double x = p->x[i],
               xp = p->xp[i-1], 
               z = p->z[i], 
               zp = p->zp[i-1],
               diff_eps = p->diff_eps[i-1];
        delta_x[i-1] += partial_x(x, xp, z, zp, diff_eps);
        delta_z[i-1] += partial_z(x, xp, z, zp, diff_eps);
        max_delta_norm = fmax(max_delta_norm, hypot(delta_x[i-1], delta_z[i-1]));
    }

    for (int i = 1; i < POINTS-1; i++) {
        p->x[i] += MOMENTUM_PARAM * momentum->x[i-1] - delta_x[i-1] / max_delta_norm * eps;
        p->z[i] += MOMENTUM_PARAM * momentum->z[i-1] - delta_z[i-1] / max_delta_norm * eps; 
    }

    push_out_of_hitboxes(p, hb);
    // remove_equal_points(p);

    // compute derivative and arclength
    double arclength[POINTS-1],
           tot_arclength = 0.0;

    compute_arclength(p, arclength, &tot_arclength);

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
    recompute_dependent(p);
}

void make_copy(struct path *dst, struct path *src) {
    array_copy(dst->x, src->x, POINTS);
    array_copy(dst->z, src->z, POINTS);
    array_copy(dst->xp, src->xp, POINTS-1);
    array_copy(dst->zp, src->zp, POINTS-1);
    array_copy(dst->diff_eps, src->diff_eps, POINTS-1);
}

struct memory init_memory(int initial_size) {
    struct memory mem = {
        initial_size,
        0,
        malloc(initial_size * POINTS * sizeof(double)),
        malloc(initial_size * POINTS * sizeof(double))
    };

    if (mem.x == NULL || mem.z == NULL) {
        puts("Could not allocate arrays for memory!");
        exit(1);
    }
    
    return mem;
}

void store_into_memory(struct path *p, struct memory *mem) {
    array_copy(mem->x + POINTS * mem->next, p->x, POINTS);
    array_copy(mem->z + POINTS * mem->next, p->z, POINTS);

    if (++mem->next >= mem->size) {
        mem->size += 100;
        printf("Memory exceeded, reallocing to %d\n", mem->size);
        mem->x = realloc(mem->x, mem->size * POINTS * sizeof(double));
        mem->z = realloc(mem->z, mem->size * POINTS * sizeof(double));
        if (mem->x == NULL || mem->z == NULL) {
            puts("Could not realloc!");
            exit(1);
        }
    }
}

void free_mem(struct memory *mem) {
    free(mem->x);
    free(mem->z);
}

void optimize(
    struct path *p,
    struct hitboxes *active_hitboxes,
    double threshold,
    double min_diff,
    int max_iters
) {
    struct path backup;
    make_copy(&backup, p);

    struct memory mem = init_memory(100);

    double obj = total_time(p),
          prev_obj = obj,
          diff = HUGE_VAL,
          eps = 1.0;
    struct momentum momentum = {{0.0}, {0.0}};
    int iters = 0;
    while (eps >= threshold && iters <= max_iters && diff > min_diff) {
        make_copy(&backup, p);
        if (iters % 100 == 0) {
            store_into_memory(p, &mem);
            printf("%d: %f\n", iters, obj);
        } 
        iters++;

        descend_and_renormalize(p, active_hitboxes, &momentum, eps);

        obj = total_time(p);
        diff = prev_obj - obj;

        if (diff < 0) { // worsening
            eps /= 2;
            printf("eps = %f\n", eps);

            for (int i = 1; i < POINTS-1; i++) {
                momentum.x[i-1] = momentum.z[i-1] = 0;
            }

            make_copy(p, &backup);
            diff = HUGE_VALF;
        } else {
            for (int i = 1; i < POINTS-1; i++) {
                momentum.x[i-1] = p->x[i] - backup.x[i];
                momentum.z[i-1] = p->z[i] - backup.z[i];
            }
            make_copy(&backup, p);
            prev_obj = obj;
        }
    }

    if (prev_obj < obj) {
        make_copy(p, &backup);
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

    double lo = 0.0, hi = 1.0, cur;
    const double eps = 1e-4;
    // binary search
    while (hi - lo > eps) {
        cur = (lo + hi) / 2;
        double dist = hypot(
            (1 - cur) * x0 + cur * x1 - ref_x,
            (1 - cur) * z0 + cur * z1 - ref_z
        );
        if (dist < target_length - eps) {
            lo = cur;
        } else if (dist > target_length + eps) {
            hi = cur;
        } else { // eps away
            break;
        }
    }

    *res_x = (1 - cur) * x0 + cur * x1;
    *res_z = (1 - cur) * z0 + cur * z1;  
}

void output_yaws(
    int length,
    double *renorm_x,
    double *renorm_z,
    double *yaws
) {
    FILE *f = fopen("yaws.txt", "w");

    unsigned short prev_au;
    for (int i = 0; i < length; i++) {
        unsigned short cur = radians_to_au(yaws[i]);
        if (i == 0) {
            fprintf(f, "%f,%f,%f,%d\n", renorm_x[i], renorm_z[i], yaws[i], cur);
        } else {
            fprintf(f, "%f,%f,%f,%d,%d\n", renorm_x[i], renorm_z[i], yaws[i], cur, cur - prev_au);
        }
        prev_au = cur;
    }

    fclose(f);
}

void compute_and_output_yaws(struct path *p) {
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
            i++;

            if (i == arr_siz) {
                arr_siz += POINTS;
                renorm_x = realloc(renorm_x, arr_siz * sizeof(double));
                renorm_z = realloc(renorm_z, arr_siz * sizeof(double));
                if (renorm_x == NULL || renorm_z == NULL) {
                    puts("Couldn't reallocate renorm arrays!");
                    exit(1);
                }
            }
        }
    }

    if (cum_arclength[POINTS-1] > i * S) {
        // add one last short segment
        i++;

        if (i == arr_siz) {
            arr_siz++;
            renorm_x = realloc(renorm_x, arr_siz * sizeof(double));
            renorm_z = realloc(renorm_z, arr_siz * sizeof(double));
            if (renorm_x == NULL || renorm_z == NULL) {
                puts("Couldn't reallocate renorm arrays!");
                exit(1);
            }
        }

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
    yaws[0] = prescribed_angle(renorm_x[0], renorm_z[0], renorm_xp[0], renorm_zp[0]);

    for (int i = 1; i < len_renorm; i++) {
        yaws[i] = prescribed_angle(renorm_x[i], renorm_z[i], renorm_xp[i-1], renorm_zp[i-1]);
    }

    output_yaws(len_renorm, renorm_x, renorm_z, yaws);

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
    struct path p = initialize_path(1374, 948, -1326, -1202, 150);
    // Right of chest 3 
    // struct path p = initialize_path_inter(1374, 948, 1000, -1000, -1326, -1202, 150);
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
    recompute_dependent(&p);

    optimize(&p, hitboxes, 1e-5f, 1e-9f, INT_MAX);

    compute_and_output_yaws(&p);

    free(hitboxes);

    return 0;
}