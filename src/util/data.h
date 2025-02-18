#pragma once

#include "../types.h"

void update_and_apply_momentum(
    struct data *d,
    struct momentum *mom,
    v2d *delta,
    double eps
);

void recompute_dependent(struct data *d);