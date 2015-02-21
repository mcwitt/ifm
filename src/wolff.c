#include "wolff.h"
#include "spin.h"
#include <stdlib.h>

void wolff_init(wolff *w)
{
    w->stack = malloc(MAX_STACK_SIZE * sizeof *w->stack);
    w->flagged = calloc(LT_N, sizeof *w->flagged);
}

void wolff_free(wolff *w)
{
    free(w->stack);
    free(w->flagged);
}

void wolff_update(lattice const *l, double p,
                  rng_state *rng, wolff *w, state *s)
{
    site_label site, *top = w->stack;
    ising_spin cluster_spin;

    /* pick a random site to seed the cluster */
    site = (site_label) (LT_N * RNG_RAND_UNIFORM(rng));
    cluster_spin = s->spin[site];
    w->cluster_size = 0;

    while (1)
    {
        neighbor_iter iter;
        int i;

        SPIN_FLIP(s->spin[site]);
        w->flagged[site] = 0;
        w->cluster_size++;
        neighbor_iter_init(site, &iter);

        /* loop over neighbors of site */
        for (i = 0; i < LT_Z; i++)
        {
            *top = neighbor_iter_next(l, &iter);

            if (s->spin[*top] == cluster_spin)
            {
                w->energy += 2;

                if ((! w->flagged[*top]) && RNG_COIN_TOSS(rng, p))
                {
                    /*
                     * add neighbor to the cluster and push it onto the stack so
                     * that its neighbors will be considered on the next iteration
                     */

                    w->flagged[*top] = 1;
                    top++;
                }
            }
            else w->energy -= 2;
        }

        if (top == w->stack) break;    /* done? */
        site = *(--top);
    }

    s->magnetization -= 2 * cluster_spin * w->cluster_size;
}
