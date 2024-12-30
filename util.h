#pragma once

#define POINTS 128
#define MOMENTUM_PARAM 0.8
#define BASE_MEM_SIZE 100

// TYPES

struct path {
    double x[POINTS], z[POINTS];
};

struct intermediates {
    double theta, current_x;
};

struct data {
    struct path p;
    // refers to points 1 through POINTS-1
    double xp[POINTS-1], zp[POINTS-1];
    // refers to points 1 through POINTS-2, padded with last value
    double xpp[POINTS-1], zpp[POINTS-1];
    double theta[POINTS-1], partial_theta[POINTS-1];
    double lagrangian[POINTS-1];
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

// FUNCTIONS

void array_copy(double *dst, double *src, int length);
void derivative(double *dst, double *src, int srclen);
void lerp(double *dst, double start, double end, int start_i, int end_i);
float opt_numdiff_eps(double x, double xp, double z, double zp);
double theta(double x, double z, double xp, double zp);
void make_copy(struct data *dst, struct data *src);
struct memory init_memory(int initial_size);
void store_into_memory(struct path *p, struct memory *mem);
void free_mem(struct memory *mem);