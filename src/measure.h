#include "ifm.h"
#include "state.h"
#include "wolff.h"
#include <stdint.h>

#define NUM_MOMENTS 2

enum
{
    M2 = 0,
    C,
    E,
    EM2,

    NUM_AVERAGES
};

typedef struct { double v[NUM_MOMENTS][NUM_AVERAGES]; } meas;

void meas_reset(meas *m);
void meas_accum(meas *sum, wolff const *w, state const *s);
void meas_average(meas *p, uint64_t t_meas);

