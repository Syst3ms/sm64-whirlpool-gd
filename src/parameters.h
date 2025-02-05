#pragma once

#define POINTS 128
#define S 28.0 // mario swimming speed
#define PENALTY_FACTOR 1000
#define ONE_HAU (M_PI/2048)
#define MAX_YAW_SPEED (5*M_PI/256)
#define OVERSMOOTH_FACTOR 0.25
#define MOMENTUM_PARAM 0.8
#define BASE_MEM_SIZE 100
#define MEM_STORE_RATE 100
#define BETA_1 0.9
#define BETA_2 0.999

#include <math.h>
#include <float.h>

#define D_EPS cbrt(DBL_EPSILON)
#define D_FAC_UP (1 + D_EPS)
#define D_FAC_DOWN (1 - D_EPS)