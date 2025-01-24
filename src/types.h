#pragma once

#include <stdalign.h>

#include "parameters.h"

struct pt_vel {
    double x, z, xp, zp, xpp, zpp, time_integrand, theta;
};

struct pt {
    double x, z;
};

typedef double v2d __attribute__((vector_size (16)));
typedef double v4d __attribute__((vector_size (32)));

struct vec_pt {
    _Alignas(16) v2d pos, vel, acc, extra;
};

union point {
    struct pt p;
    struct vec_pt v;
    struct pt_vel pv;
};

struct data {
    union point points[POINTS];
};

struct mom_point {
    _Alignas(16) v2d xz, ut;
};

struct momentum {
    struct mom_point points[POINTS-2];
};

struct hitbox {
    _Alignas(16) v2d pos;
    double radius;
};

struct hitboxes {
    int num_hb;
    struct hitbox hb[];
};

struct memory {
    int size, next;
    struct pt *pts;
};

struct ad1d {
    double v, d;
};