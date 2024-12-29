#define POINTS 128
#define S 28.0
#define MOMENTUM_PARAM 0.8
#define BASE_MEM_SIZE 100
#define MARIO_HITBOX_SIZE 100

struct path {
    double x[POINTS], z[POINTS];
    double xp[POINTS-1], zp[POINTS-1];
    double diff_eps[POINTS-1];
};

struct momentum {
    double x[POINTS-2], z[POINTS-2];
};

struct hitbox {
    double x, z, radius;
};

struct hitboxes {
    int num_hb;
    struct hitbox hb[];
};

struct memory {
    int size, next;
    double *x, *z;
};

// functions

void array_copy(double *dst, double *src, int length);
void derivative(double *dst, double *src, int srclen);
void lerp(double *dst, double start, double end, int start_i, int end_i);
float opt_numdiff_eps(double x, double xp, double z, double zp);
void recompute_dependent(struct path *p);