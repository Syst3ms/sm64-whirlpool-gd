#pragma once

#define POINTS 128
#define S 28.0 // mario swimming speed
#define PENALTY_FACTOR 1000
#define OVERSMOOTH 0.5
#define MAX_YAW_SPEED (OVERSMOOTH*(5*M_PI/256))
#define MAX_YAW_SPEED_SQUARED (MAX_YAW_SPEED * MAX_YAW_SPEED)
#define MOMENTUM_PARAM 0.8
#define BASE_MEM_SIZE 100
#define MEM_STORE_RATE 1000
#define BETA_1 0.9
#define BETA_2 0.999
#define MIN_CONSTRAINT_PROGRESS 0.99
#define CONSTRAINT_UPSCALE 1.5
#define EPS_DOWNSCALE 0.5
#define MAX_SHIFT 100.0

#include <math.h>
#include <float.h>

#define D_EPS cbrt(DBL_EPSILON)
#define D_FAC_UP (1 + D_EPS)
#define D_FAC_DOWN (1 - D_EPS)