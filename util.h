#pragma once

#define POINTS 128
#define S 28.0
#define MOMENTUM_PARAM 0.8
#define BASE_MEM_SIZE 100

// TYPES

// TODO
/*
 * Refactor path to only have x and z, and move all of the dependents
 * to a new struct called 'data'. This struct will contain: xp, zp, theta_star_p,
 * diff_eps, lagrangian, intermediates
 * 
 * 'intermediates' is a new struct that stores theta_star and current_x
 * 
 * When recomputing dependents, first run the first half of the lagrangian
 * to update intermediates, then take the derivative, then finish
 * computing the lagrangian using intermediates
 * 
 * The direct lagrangian for partials is separate from the two-step one above.
 * 
 * Objective function now incurs no computational cost
 */

struct path {
    double x[POINTS], z[POINTS];
    // correspond to []
    double xp[POINTS-1], zp[POINTS-1];
    // correspond to positions [1..POINTS-1] and [1..POINTS-2] respectively
    double theta_star[POINTS-1], theta_star_p[POINTS-2];
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

// FUNCTIONS

void array_copy(double *dst, double *src, int length);
void derivative(double *dst, double *src, int srclen);
void lerp(double *dst, double start, double end, int start_i, int end_i);
float opt_numdiff_eps(double x, double xp, double z, double zp);
void recompute_dependent(struct path *p);
double theta_star(double x, double z, double xp, double zp);
void make_copy(struct path *dst, struct path *src);
struct memory init_memory(int initial_size);
void store_into_memory(struct path *p, struct memory *mem);
void free_mem(struct memory *mem);