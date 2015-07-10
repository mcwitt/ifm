#include "ifm.h"
#include "lattice.h"
#include "measure.h"
#include "rng.h"
#include "spin.h"
#include "state.h"
#include "wolff.h"
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* probability to add a site to the cluster at a given temperature T */
#define WOLFF_P_ADD(T)  (1. - exp(-2./(T)))

/* do Monte Carlo step and update number of spin flips */
void mc_step(lattice *l,
             double p,
             rng_state *rng,
             wolff *w,
             state *s,
             uint64_t *num_flips)
{
    int i;

    for (i = 0; i < UPDATES_PER_STEP; i++)
    {
        wolff_update(l, p, rng, w, s);
        *num_flips += w->cluster_size;
    }
}

/* get initial estimate of equilibrium cluster size for sweep calibration */
double bootstrap(lattice *l, double p, uint64_t num_steps,
                 rng_state *rng, wolff *w, state *s)
{
    uint64_t i, num_flips = 0;

    for (i = 0; i < num_steps; i++)
        mc_step(l, p, rng, w, s, &num_flips);

    /* return mean steps per sweep */
    return num_steps / ((double) num_flips/LT_N);
}

/* update and measure */
double measure(lattice *l, double p, uint64_t num_steps,
               rng_state *rng, wolff *w, state *s, meas *ms)
{
    uint64_t i, num_flips = 0;

    meas_reset(ms);

    for (i = 0; i < num_steps; i++)
    {
        mc_step(l, p, rng, w, s, &num_flips);
        meas_accum(ms, w, s);
    }

    meas_average(ms, num_steps);

    /* return mean steps per sweep */
    return num_steps / ((double) num_flips/LT_N);
}

int main(int argc, char *argv[])
{
    lattice l;
    state s;
    wolff w;
    meas ms;
    rng_state *rng;
    double T;
    int i, seed, shape[] = {[0 ... LT_D-1] = LT_LMAX};

    seed = (argc > 1) ? atoi(argv[1]) : time(NULL);
    lattice_init(&l, shape);

    if (state_alloc(&s) | wolff_alloc(&w))
    {
        int err = errno;
        fprintf(stderr, "%s: %s\n", argv[0], strerror(err));
        return EXIT_FAILURE;
    }

    for (i = 0; i < LT_N; i++) s.spin[i] = SPIN_UP;
    s.magnetization = LT_N;
    w.energy =  LT_E1;

    rng = RNG_ALLOC();
    RNG_SEED(rng, seed);

    printf("%9s %12s " \
           "%12s %12s " \
           "%12s %12s " \
           "%12s %12s " \
           "%12s %12s\n",
           "Temp", "LogSweeps",
           "M2", "M4",
           "C", "C2",
           "U", "U2",
           "UM2", "U2M4");

    while (fscanf(stdin, "%lf", &T) != EOF) /* loop over temperatures */
    {
        double steps_per_sweep, p = WOLFF_P_ADD(T);
        uint64_t bootstrap_steps = pow(2, LOG_BOOTSTRAP_STEPS);

        /* estimate equilibrium cluster size */
        steps_per_sweep = \
            bootstrap(&l, p, bootstrap_steps, rng, &w, &s);

        /* update/measure loop */
        for (i = LOG_BOOTSTRAP_STEPS; i <= LOG_MEAS_TIME; i++)
        {
            uint64_t meas_sweeps = pow(2, i),
                     meas_steps = meas_sweeps * steps_per_sweep;

            steps_per_sweep = \
                measure(&l, p, meas_steps, rng, &w, &s, &ms);

            printf("%9g %12d " \
                   "%12g %12g " \
                   "%12g %12g " \
                   "%12g %12g " \
                   "%12g %12g\n",
                    T, i,
                    ms.v[0][M2],  ms.v[1][M2],
                    ms.v[0][C],   ms.v[1][C],
                    ms.v[0][U],   ms.v[1][U],
                    ms.v[0][UM2], ms.v[1][UM2]);

            fflush(stdout);
        }
    }

    state_free(&s);
    wolff_free(&w);
    RNG_FREE(rng);

    return EXIT_SUCCESS;
}
