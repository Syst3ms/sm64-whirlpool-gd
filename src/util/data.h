#pragma once

#include "../types.h"

struct history init_history(size_t initial_size);
void store_into_history(struct data *d, struct history *mem);
void free_history(struct history mem);

unsigned short radians_to_au(double rad);
double fast_hypot(double dx, double dz);
double fast_hypot_v(v2d v);
void update_and_apply_momentum(
    struct data *d,
    struct momentum *mom,
    v2d *delta,
    double eps
);

double objective(struct data *d, struct penalty_data *pdata);
void recompute_dependent(struct data *d);