#include "util.h"

#include <float.h>
#include <math.h>
#include <memory.h>

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

float opt_numdiff_eps(double x, double xp, double z, double zp) {
    return sqrt(DBL_EPSILON) * (fabs(x) + fabs(xp) + fabs(z) + fabs(zp) + sqrt(DBL_EPSILON));
}

void recompute_dependent(struct path *p) {
    derivative(p->xp, p->x, POINTS);
    derivative(p->zp, p->z, POINTS);

    for (int i = 1; i < POINTS; i++) {
        p->diff_eps[i-1] = opt_numdiff_eps(p->x[i], p->xp[i-1], p->z[i], p->zp[i-1]);
    }
}