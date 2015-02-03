#include "state.h"

#define NUM_MOMENTS 2

enum
{
    M2 = 0,
    C,
    NUM_AVERAGES
};

typedef struct { double v[NUM_MOMENTS][NUM_AVERAGES]; } meas_basic;

void meas_basic_reset(meas_basic *m);
void meas_basic_accum(meas_basic *m, state const *s);
void meas_basic_average(meas_basic *p, unsigned long num_updates);

