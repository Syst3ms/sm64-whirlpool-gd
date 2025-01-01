#include "util.h"

#define _USE_MATH_DEFINES
#include <float.h>
#include <math.h>
#include <memory.h>
#include <stdlib.h>
#include <stdio.h>

inline void array_copy(double *dst, double *src, int length) {
    memcpy(dst, src, length * sizeof(double));
}

inline void derivative(double *dst, double *src, int srclen) {
    for (int i = 0; i < srclen-1; i++) {
        dst[i] = (src[i+1] - src[i]) * POINTS;
    }
}

inline double lerp_factor(double start, double end, double target) {
    return (target - start) / (start - end);
}

// lerp between start and end on [start_i..end_i] (inclusive)
void lerp(double *dst, double start, double end, int start_i, int end_i) {
    double step = (end - start) / (end_i - start_i);
    for (int i = 0; i <= end_i - start_i; i++) {
        dst[start_i + i] = start + step * i;
    }
}

inline void make_copy(struct data *dst, struct data *src) {
    *dst = *src;
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

inline void free_mem(struct memory *mem) {
    free(mem->x);
    free(mem->z);
}

inline unsigned short radians_to_au(double rad) {
    unsigned short a = (unsigned short) ((rad / M_PI) * 0x8000);
    return a & 0xFFFFFFF0;
}

void update_and_apply_momentum(
    struct path *p,
    struct momentum *mom,
    double *delta_x,
    double *delta_z,
    double eps,
    char first
) {
    if (first) {
        for (int i = 0; i < POINTS-2; i++) {
            mom->x[i] = mom->ut_x[i] = delta_x[i];
            mom->z[i] = mom->ut_z[i] = delta_z[i];
        }
    }
    double mx, mz, ux, uz;
    for (int i = 0; i < POINTS-2; i++) {
        mx = mom->x[i] = BETA_1 * mom->x[i] + (1 - BETA_1) * delta_x[i];
        mz = mom->z[i] = BETA_1 * mom->z[i] + (1 - BETA_1) * delta_z[i];
        ux = mom->ut_x[i] = fmax(BETA_2 * mom->ut_x[i], fabs(delta_x[i]));
        uz = mom->ut_z[i] = fmax(BETA_2 * mom->ut_z[i], fabs(delta_z[i]));
        p->x[i+1] -= eps * mx / ux;
        p->z[i+1] -= eps * mz / uz;
    }
}