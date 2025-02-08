#pragma once

#include <stdalign.h>
#include <stddef.h>

#include "parameters.h"

typedef double v2d __attribute__((vector_size (16)));

union point {
    struct {
        double x, z, xp, zp, xpp, zpp, time_integrand, lagrangian; 
    };
    struct {
        _Alignas(16) v2d pos, vel, acc, extra;
    };
};

struct penalty_data {
    double shift[POINTS-2];
    double rho;
};

struct data {
    union point points[POINTS];
    double constraint[POINTS-2];
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
    size_t num_hb;
    struct hitbox hb[];
};

struct history {
    size_t size, next;
    v2d *pts;
};