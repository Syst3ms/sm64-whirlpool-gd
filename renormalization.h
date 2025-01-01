#pragma once

void compute_output_resampled(struct data *d);
void resample_time_to_frame(struct data *d, int *res_len, double **res_x, double **res_z);