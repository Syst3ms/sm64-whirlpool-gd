#include <stdlib.h>
#include <mm_malloc.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <errno.h>

#include "parameters.h"
#include "util.h"
#include "math_funcs.h"

#define TIME_PER_FRAME 1

void output_resampled(int length, v2d *p, double *yaws) {
    FILE *f = fopen(".\\resampled.txt", "w");
    if (f == NULL) {
        printf("Couldn't open file for resampled: %d\n", errno);
        exit(1);
    }

    fprintf(f, "x,z,yaw,yaw_diff\n");
    unsigned short prev_au = radians_to_au(yaws[0]);
    for (int i = 0; i < length; i++) {
        unsigned short cur = radians_to_au(yaws[i]);
        short diff = cur - prev_au;
        fprintf(
            f, "%f,%f,%d,%d%s\n",
            p[i][0], p[i][1], cur, diff, abs(diff) > 640 ? "!" : ""
        );
        prev_au = cur;
    }

    fclose(f);
}

// finds the point between `first` and `second` that is `target_length` away
// from `ref`, and stores it in `res`
void find_resampled(v2d ref, v2d first, v2d second, double target_length, v2d *res) {
    if (ref[0] == first[0] && ref[1] == first[1]) {
        // direct lerp
        double fac = target_length / length_v(second - first);
        *res = (1-fac) * first + fac * second;
        return; 
    }

    v2d v = second - first, w = first - ref;
    v2d a = {v[0], w[0]}, b = {v[1], w[1]};
    v2d sums = a * a + b * b;

    double snv = sums[0];
    double snw = sums[1];
    double dp = v[0] * w[0] + v[1] * w[1];

    double fac = (-dp + sqrt(dp * dp + snv * (target_length * target_length - snw))) / snv;
    *res = (1-fac) * first + fac * second;
}

void * _mm_realloc(void *aligned_ptr, size_t prev_size, size_t size, size_t align) {
    void *new = _mm_malloc(size, align);
    memcpy(new, aligned_ptr, prev_size);
    _mm_free(aligned_ptr);
    return new;
}

v2d * inc_and_realloc_if_necessary(int *i, int *cap, v2d **arr) {
    (*i)++;

    if (*i == *cap) {
        *cap += POINTS;
        v2d *new_ptr = _mm_realloc(*arr, *i * sizeof(v2d), *cap * sizeof(v2d), 16);
        if ((*arr = new_ptr) == NULL) {
            puts("Couldn't reallocate renorm arrays!");
            exit(1);
        }
    }

    return *arr;
}

void resample_time_to_frame(struct data *d, int *res_len, v2d **new_array) {
    // resample the path so that each segment takes exactly one frame

    int capacity = POINTS;
    v2d *new = *new_array;

    double time_so_far = 0.0;
    new[0] = d->points[0].v.pos;
    v2d last_renorm = d->points[0].v.pos;
    int i = 1; // next point to place
    int last_sample = 0; // k such that the last sampled is between k and k+1
    double outgoing_speed = real_speed_norm(d->points[0].v.pos, d->points[0].v.vel);

    for (int j = 0; j < POINTS-1; j++) {
        if (last_sample == POINTS-2) {
            break;
        }
        time_so_far += d->points[j].pv.time_integrand / (POINTS-1);

        while (time_so_far >= i * TIME_PER_FRAME) {
            int k = last_sample;
            v2d prev_pos = d->points[k].v.pos;
            v2d next_pos = d->points[k+1].v.pos;
            while (k < POINTS-2 && length_v(next_pos - last_renorm) < outgoing_speed * TIME_PER_FRAME)
            {
                k++;
                prev_pos = next_pos;
                next_pos = d->points[k+1].v.pos;
            }
            
            v2d new_p;
            find_resampled(last_renorm, prev_pos, next_pos, outgoing_speed * TIME_PER_FRAME, &new_p);
            outgoing_speed = real_speed_norm(
                new_p,
                (new_p - last_renorm) * POINTS
            );
            last_renorm = new[i] = new_p;
            last_sample = k;
            new = inc_and_realloc_if_necessary(&i, &capacity, new_array);
        }
    }

    // don't add final point to ensure the end is the same, it's not that important
    *res_len = i;
}

void compute_output_resampled(struct data *d) {
    int len;
    v2d *new = _mm_malloc(POINTS * sizeof(v2d), 16);

    resample_time_to_frame(d, &len, &new);

    v2d *vel = _mm_malloc(len * sizeof(v2d), 16);
    double *yaws = malloc(len * sizeof(double));
    if (vel == NULL || yaws == NULL) {
        puts("Couldn't allocate xp, zp or yaws!");
        exit(1);
    }

    for (int i = 1; i < len; i++) {
        vel[i] = (new[i] - new[i-1]) * (double) len;
    }
    vel[0] = vel[1];

    // now that we have proper step sizes, we can find out our prescribed yaws
    for (int i = 0; i < len; i++) {
        yaws[i] = theta(new[i], vel[i]);
    }

    puts("Writing resampled");
 
    output_resampled(len, new, yaws);

    puts("Done with resampled");

    free(yaws);
    _mm_free(new);
    _mm_free(vel);
}