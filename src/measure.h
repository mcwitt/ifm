#include "state.h"

#ifndef LOMEM
#include "wolff.h"
#else
#include "wolff_lomem.h"
typedef wolff_lomem wolff;
#endif

#define NUM_MOMENTS 2

enum
{
    M2 = 0,
    C,
#ifndef LOMEM
    E,
    EM2,
#endif

    NUM_AVERAGES
};

typedef struct { double v[NUM_MOMENTS][NUM_AVERAGES]; } meas;

void meas_reset(meas *m);
void meas_accum(meas *sum, wolff const *w, state const *s);
void meas_average(meas *p, unsigned int t_meas);

