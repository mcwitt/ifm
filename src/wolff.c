#include "wolff.h"
#include "spin.h"
#include <stdlib.h>

int wolff_alloc(wolff *w)
{
    if ((w->stack = malloc(MAX_STACK_SIZE * sizeof *w->stack)) &&
        (w->flagged = calloc(LT_N, sizeof *w->flagged)))
        return 0;

    return 1;   /* oops, couldn't allocate memory... */
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

        SPIN_FLIP(s->spin[site]);
        w->flagged[site] = 0;
        w->cluster_size++;
        neighbor_iter_init(site, &iter);

        /* loop over neighbors */
        while ((*top = neighbor_iter_next(l, &iter)) != -1)
        {
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
