#include <float.h>
#include <math.h>
#include <memory.h>
#include <stdlib.h>
#include <stdio.h>
#include <immintrin.h>
#include <stddef.h>

#include "util.h"
#include "autodiff.h"

#define M_PI		3.14159265358979323846

struct history init_history(size_t initial_size) {
    struct history mem = {
        initial_size,
        0,
        malloc(initial_size * POINTS * sizeof(v2d))
    };

    if (mem.pts == NULL) {
        puts("Could not allocate arrays for history!");
        exit(1);
    }
    
    return mem;
}

void store_into_history(struct data *d, struct history *mem) {
    size_t offset = POINTS * mem->next;
    for (size_t i = 0; i < POINTS; i++) {
        mem->pts[offset + i] = d->points[i].pos;
    }

    if (++mem->next >= mem->size) {
        mem->size += 1000;
        printf("Memory exceeded, reallocing to %d\n", mem->size);
        mem->pts = realloc(mem->pts, mem->size * POINTS * sizeof(v2d));
        if (mem->pts == NULL) {
            puts("Could not realloc!");
            exit(1);
        }
    }
}

inline void free_history(struct history mem) {
    free(mem.pts);
}

unsigned short radians_to_au(double rad) {
    unsigned short a = (unsigned short) ((rad / M_PI) * 0x8000);
    return a & 0xFFF0;
}

inline double fast_hypot(double dx, double dz) {
    return sqrt(dx * dx + dz * dz);
}

double fast_hypot_v(v2d v) {
    return sqrt(v[0] * v[0] + v[1] * v[1]);
}

v2d vabs(v2d x) {
    v2d mask = {-0.0, -0.0};
    return _mm_andnot_pd(mask, x);
}

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

