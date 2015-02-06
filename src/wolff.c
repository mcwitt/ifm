#include "lattice.h"
#include "rng.h"
#include "spin.h"
#include "state.h"
#include "update.h"
#include "wolff.h"
#include "meas_basic.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* probability to add a site to the cluster at a given temperature T */
#define WOLFF_P_ADD(T)  (1. - exp(-2./(T)))

/* Monte Carlo step (one unit of time) */
void mc_step(lattice *l, double p, rng *rand, state *s)
{
    int i;
    for (i = 0; i < UPDATES_PER_STEP; i++) update(l, p, rand, s);
}

void thermalize(lattice *l, double p, int t_therm,
                rng *rand, state *s)
{
    int i;
    for (i = 0; i < t_therm; i++) mc_step(l, p, rand, s);
}

void measure(lattice *l, double p, int t_meas,
             rng *rand, state *s, meas_basic *meas)
{
    int i;

    meas_basic_reset(meas);

    for (i = 0; i < t_meas; i++)
    {
        mc_step(l, p, rand, s);
        meas_basic_accum(meas, s);
    }

    meas_basic_average(meas, t_meas);
}

int main(int argc, char *argv[])
{
    lattice l;
    state s;
    meas_basic meas;
    rng *rand;
    double T;
    int i, seed;

    seed = (argc > 1) ? atoi(argv[1]) : time(NULL);
    lattice_init(&l);

    state_alloc(&s);
    for (i = 0; i < LT_N; i++) s.spin[i] = SPIN_UP;
    s.magnetization = LT_N;

    rand = RNG_ALLOC();
    RNG_SEED(rand, seed);

#ifndef DEBUG
    printf("%9s %6s " \
           "%12s %12s " \
           "%12s %12s\n",
           "Temp", "Stage",
           "M2", "M4",
           "C", "C2");
#endif

    while (fscanf(stdin, "%lf", &T) != EOF)
    {
        double p = WOLFF_P_ADD(T);
        int dt = pow(2, LOG_TIME_THERM);

        thermalize(&l, p, dt, rand, &s);

        for (i = LOG_TIME_THERM; i < LOG_TIME_MEAS; i++)
        {
            dt = pow(2, i);
            measure(&l, p, dt, rand, &s, &meas);

#ifndef DEBUG
            printf("%9g %6d " \
                   "%12g %12g " \
                   "%12g %12g\n",
                    T, i,
                    meas.v[0][M2], meas.v[1][M2],
                    meas.v[0][C],  meas.v[1][C]);

            fflush(stdout);
#endif
        }
    }

    state_free(&s);
    RNG_FREE(rand);

    return EXIT_SUCCESS;
}
