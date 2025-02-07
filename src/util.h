#pragma once

#include "types.h"

struct memory init_memory(int initial_size);
void store_into_memory(struct data *d, struct memory *mem);
void free_memory(struct memory *mem);

unsigned short radians_to_au(double rad);
double fast_hypot(double dx, double dz);
double fast_hypot_v(v2d v);
void update_and_apply_momentum(
    struct data *d,
    struct momentum *mom,
    v2d *delta,
    double eps
);

double objective(struct data *d);
void recompute_dependent(struct data *d);