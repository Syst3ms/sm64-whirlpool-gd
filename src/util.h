#pragma once

#include "types.h"

struct memory init_memory(int initial_size);
void store_into_memory(struct data *d, struct memory *mem);
void free_memory(struct memory *mem);

unsigned short radians_to_au(double rad);
double length(double dx, double dz);
double length_v(v2d v);
v2d ortho(v2d a);
double remainder_2pi(double x);
int is_sign_same(double x, double y);
double move_towards(double base, double direction, double amount);
void update_and_apply_momentum(
    struct data *d,
    struct momentum *mom,
    v2d *delta,
    double eps
);

void recompute_speeds(struct data *d);
void recompute_dependent(struct data *d);