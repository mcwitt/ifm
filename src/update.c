#include "update.h"
#include "spin.h"

void update(lattice const *l, double p, rng *rand, state *s)
{
    lattice_site site, *top = s->stack;

    /* pick a random site to seed the cluster */
    lattice_random_site(l, rand, &site); /* site <- random site */
    s->cluster_spin = s->spin[site.i];

    SPIN_FLIP(s->spin[site.i]);
    s->cluster_size = 1;

    while (1)
    {
        int n;

        /* loop over neighbors of site */
        for (n = 0; n < LT_Z; n++)
        {
            /* TODO make this more efficient:
             * neighbor coords aren't needed until we add to cluster */
            /* top <- nth neighbor */
            int is_neighbor = lattice_neighbor(l, n, &site, top);
            if (! is_neighbor) continue;

            if ((s->spin[top->i] == s->cluster_spin) \
                && (RNG_COIN_TOSS(rand, p)))
            {
                /*
                 * add neighbor to the cluster and push it onto the stack so
                 * that its neighbors will be considered in the next iteration
                 */

                SPIN_FLIP(s->spin[top->i]);
                s->cluster_size++;
                top++;
            }
        }

        if (top == s->stack) break;    /* done? */
        site = *(--top);
    }

    s->magnetization -= 2 * s->cluster_spin * s->cluster_size;
}
