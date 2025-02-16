#pragma once

#define LAGR_PARTIAL_HEADER(var)\
    double lagr_partial_##var(union point *pt, double penalty_fac, double shift)

#define LAGR_PARTIAL(var) \
LAGR_PARTIAL_HEADER(var) {\
    double orig = pt->var;\
    pt->var *= D_FAC_UP;\
    double lagr_up = lagrangian(pt, penalty_fac, shift);\
    pt->var = orig;\
    pt->var *= D_FAC_DOWN;\
    double lagr_down = lagrangian(pt, penalty_fac, shift);\
    pt->var = orig;\
    return (lagr_up - lagr_down) / (2 * D_EPS * orig);\
}

#define LAGR_PARTIAL_CHECK_ZERO(var) \
LAGR_PARTIAL_HEADER(var) {\
    double orig = pt->var;\
    if (orig != 0.0) {\
        double orig = pt->var;\
        pt->var *= D_FAC_UP;\
        double lagr_up = lagrangian(pt, penalty_fac, shift);\
        pt->var = orig;\
        pt->var *= D_FAC_DOWN;\
        double lagr_down = lagrangian(pt, penalty_fac, shift);\
        pt->var = orig;\
        return (lagr_up - lagr_down) / (2 * D_EPS * orig);\
    } else {\
        pt->var = D_EPS;\
        double lagr_up = lagrangian(pt, penalty_fac, shift);\
        pt->var = -D_EPS;\
        double lagr_down = lagrangian(pt, penalty_fac, shift);\
        pt->var = 0.0;\
        return (lagr_up - lagr_down) / (2 * D_EPS);\
    }\
}
