#include "util.h"

#include <float.h>
#include <math.h>
#include <memory.h>
#include <stdlib.h>
#include <stdio.h>

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

void make_copy(struct data *dst, struct data *src) {
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

void free_mem(struct memory *mem) {
    free(mem->x);
    free(mem->z);
}