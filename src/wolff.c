#include "wolff.h"
#include "spin.h"
#include <string.h>
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
    lattice_site site, *top = w->stack;
    ising_spin cluster_spin;

    /* pick a random site to seed the cluster */
    lattice_random_site(l, rng, &site); /* site <- random site */
    cluster_spin = s->spin[site.i];
    w->cluster_size = 0;

    while (1)
    {
        int n;

        SPIN_FLIP(s->spin[site.i]);
        w->flagged[site.i] = 0;
        w->cluster_size++;

        /* loop over neighbors of site */
        for (n = 0; n < LT_Z; n++)
        {
            /*
             * TODO: It should be more efficient to split lattice_neighbor
             * into 2 operations: one that returns just the index, and
             * another that computes the coordinates. The latter is only
             * necessary if the neighbor site is to be added to the cluster.
             */
            if (! lattice_neighbor(l, n, &site, top)) continue;

            if (s->spin[top->i] == cluster_spin)
            {
                w->energy += 2;

                if ((! w->flagged[top->i]) && RNG_COIN_TOSS(rng, p))
                {
                    /*
                     * add neighbor to the cluster and push it onto the stack so
                     * that its neighbors will be considered on the next iteration
                     */

                    w->flagged[top->i] = 1;
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
