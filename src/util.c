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

struct memory init_memory(int initial_size) {
    struct memory mem = {
        initial_size,
        0,
        malloc(initial_size * POINTS * sizeof(struct pt))
    };

    if (mem.pts == NULL) {
        puts("Could not allocate arrays for memory!");
        exit(1);
    }
    
    return mem;
}

void store_into_memory(struct data *d, struct memory *mem) {
    int offset = POINTS * mem->next;
    for (int i = 0; i < POINTS; i++) {
        mem->pts[offset + i] = d->points[i].p;
    }

    if (++mem->next >= mem->size) {
        mem->size += 1000;
        printf("Memory exceeded, reallocing to %d\n", mem->size);
        mem->pts = realloc(mem->pts, mem->size * POINTS * sizeof(struct pt));
        if (mem->pts == NULL) {
            puts("Could not realloc!");
            exit(1);
        }
    }
}

inline void free_memory(struct memory *mem) {
    free(mem->pts);
}

unsigned short radians_to_au(double rad) {
    unsigned short a = (unsigned short) ((rad / M_PI) * 0x8000);
    return a & 0xFFF0;
}

inline double length(double dx, double dz) {
    return sqrt(dx * dx + dz * dz);
}

double length_v(v2d v) {
    return sqrt(v[0] * v[0] + v[1] * v[1]);
}

v2d vabs(v2d x) {
    v2d mask = {-0.0, -0.0};
    return _mm_andnot_pd(mask, x);
}

// returns a vector orthogonal to a
v2d ortho(v2d a) {
    v2d ret = {a[1], a[0]};
    v2d mask = {0.0, -0.0};
    return _mm_xor_pd(ret, mask);
}

double remainder_2pi(double x) {
    v2d y = _mm_load_sd(&x);
    v2d div = y / (2 * M_PI);
    return (y - _mm_round_pd(div, 0) * 2 * M_PI)[0];
}

// checks if the sign of x and y is different using bitwise manipulation
int is_sign_same(double x, double y) {
    union {
        unsigned long long l;
        double d;
    } a, b;
    a.d = x;
    b.d = y;
    return (a.l ^ b.l) >> 63 == 0;
}

double move_towards(double base, double direction, double amount) {
    return base + copysign(amount, direction - base);
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
        d->points[i+1].v.pos -= eps * m / u;
    }
}

void recompute_speeds(struct data *d) {
    d->points[0].v.vel = (d->points[1].v.pos - d->points[0].v.pos) * POINTS;
    for (int i = 1; i < POINTS; i++) {
        d->points[i].v.vel = (d->points[i].v.pos - d->points[i-1].v.pos) * POINTS;
    }
} 

void recompute_dependent(struct data *d) {
    recompute_speeds(d);

    for (int i = 0; i < POINTS-1; i++) {
        compute_lagrangian_with_intermediates(&d->points[i].pv);
    }
}

