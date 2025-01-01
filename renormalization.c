#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "parameters.h"
#include "util.h"
#include "math_funcs.h"

#define TIME_PER_FRAME 1

void output_resampled(
    int length,
    double *x,
    double *z,
    double *yaws
) {
    FILE *f = fopen("resampled.txt", "w");

    fprintf(f, "x,z,yaw,yaw_diff\n");
    unsigned short prev_au = radians_to_au(yaws[0]);
    fprintf(f, "%f,%f,%d\n", x[0], z[0], prev_au);

    for (int i = 1; i < length; i++) {
        unsigned short cur = radians_to_au(yaws[i]);
        short diff = cur - prev_au;
        fprintf(
            f, "%f,%f,%d,%d%s\n",
            x[i], z[i], cur, diff, abs(diff) > 640 ? "!" : ""
        );
        prev_au = cur;
    }

    fclose(f);
}

// finds the point between (x0,z0) and (x1,z1) that is target_length away
// from (ref_x, ref_z), and stores it in (res_x,res_z)
void find_resampled(
    double ref_x, double ref_z,
    double x0, double z0, double x1, double z1,
    double target_length,
    double *res_x, double *res_z
) {
    if (x0 == ref_x && z0 == ref_z) {
        // direct lerp
        double fac = target_length / hypot(x1 - x0, z1 - z0);
        *res_x = (1 - fac) * x0 + fac * x1;
        *res_z = (1 - fac) * z0 + fac * z1;
        return; 
    }

    double vx = x1 - x0, vz = z1 - z0, wx = x0 - ref_x, wz = z0 - ref_z;

    double snv = vx * vx + vz * vz;
    double dp = vx * wx + vz * wz;
    double snw = wx * wx + wz * wz;

    double fac = (-dp + sqrt(dp * dp + snv * (target_length * target_length - snw))) / snv;

    *res_x = (1 - fac) * x0 + fac * x1;
    *res_z = (1 - fac) * z0 + fac * z1;  
}

void inc_and_realloc_if_necessary(int *i, int *cap, double **x, double **z) {
    (*i)++;

    if (*i == *cap) {
        *cap += POINTS;
        *x = realloc(*x, *cap * sizeof(double));
        *z = realloc(*z, *cap * sizeof(double));
        if (x == NULL || z == NULL) {
            puts("Couldn't reallocate renorm arrays!");
            exit(1);
        }
    }
}

void resample_time_to_frame(struct data *d, int *res_len, double **res_x, double **res_z) {
    // resample the path so that each segment takes exactly one frame
    struct path *p = &d->p;

    int capacity = POINTS;
    double *new_x = malloc(POINTS * sizeof(double)),
           *new_z = malloc(POINTS * sizeof(double));
    double time_so_far = 0.0;
    double prev_x = new_x[0] = p->x[0];
    double prev_z = new_z[0] = p->z[0];
    int i = 1; // next point to place
    int last_sample = 0; // k such that the last sampled is between k and k+1
    double outgoing_speed = real_speed_norm(new_x[0], d->xp[0], new_z[0], d->zp[0]);

    for (int j = 0; j < POINTS-1; j++) {
        if (last_sample == POINTS-2) {
            break;
        }
        time_so_far += d->time_int[j] / (POINTS-1);

        while (time_so_far >= i * TIME_PER_FRAME) {
            int k = last_sample;
            while (k < POINTS-2 && hypot(p->x[k+1] - prev_x, p->z[k+1] - prev_z) < outgoing_speed * TIME_PER_FRAME)
            {
                k++;
            }
            
            double nx, nz;
            find_resampled(
                prev_x, prev_z, p->x[k], p->z[k], p->x[k+1], p->z[k+1],
                outgoing_speed * TIME_PER_FRAME,
                &nx, &nz
            );
            outgoing_speed = real_speed_norm(
                nx,
                (nx - prev_x) * POINTS,
                nz,
                (nz - prev_z) * POINTS
            );
            prev_x = new_x[i] = nx;
            prev_z = new_z[i] = nz;
            last_sample = k;
            inc_and_realloc_if_necessary(&i, &capacity, &new_x, &new_z);
        }
    }

    // don't add final point to ensure the end is the same, it's not that important
    *res_len = i;
    *res_x = new_x;
    *res_z = new_z;
}

void compute_output_resampled(struct data *d) {
    int len;
    double *new_x, *new_z;

    resample_time_to_frame(d, &len, &new_x, &new_z);

    double *xp = malloc((len-1) * sizeof(double)),
           *zp = malloc((len-1) * sizeof(double)),
           *yaws = malloc(len * sizeof(double));
    if (xp == NULL || zp == NULL || yaws == NULL) {
        puts("Couldn't allocate xp, zp or yaws!");
        exit(1);
    }

    derivative(xp, new_x, len);
    derivative(zp, new_z, len);

    // now that we have proper step sizes, we can find out our prescribed yaws
    yaws[0] = theta(new_x[0], xp[0], new_z[0], zp[0]);

    for (int i = 1; i < len; i++) {
        yaws[i] = theta(new_x[i], xp[i-1], new_z[i], zp[i-1]);;
    }
 
    output_resampled(len, new_x, new_z, yaws);

    free(yaws);
    free(new_x);
    free(new_z);
    free(xp);
    free(zp);
}