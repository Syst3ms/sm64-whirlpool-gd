#include <immintrin.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include "./math.h"

unsigned short radians_to_au(double rad) {
    unsigned short a = (unsigned short) ((rad / M_PI) * 0x8000);
    return a & 0xFFF0;
}

inline double fast_hypot(double dx, double dz) {
    return sqrt(dx * dx + dz * dz);
}

double fast_hypot_v(v2d v) {
    return sqrt(v[0] * v[0] + v[1] * v[1]);
}

v2d vabs(v2d x) {
    v2d mask = {-0.0, -0.0};
    return _mm_andnot_pd(mask, x);
}

double square(double x) {
    return x * x;
}