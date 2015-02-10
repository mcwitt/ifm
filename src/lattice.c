#include "lattice.h"
#include "neighbors.h"
#include <string.h>

/* define dimensions of the lattice */

/* default: L x L x ... x L hypercube */
int const lt_shape[LT_D] = {[0 ... LT_D-1] = LT_LMAX};

static void init_strides(int *stride)
{
    int d;

    stride[LT_D-1] = 1;

    for (d = LT_D-1; d > 0; d--)
        stride[d-1] = lt_shape[d] * stride[d];
}

static void index2coords(int i, int const *stride, int *x)
{
    int d;

    for (d = 0; d < LT_D-1; d++)
    {
        /* TODO check whether division and remainder are optimized to one
         * instruction */
        x[d] = i / stride[d];
        i %= stride[d];
    }

    x[LT_D-1] = i;
}

void lattice_init(lattice *l)
{
    init_strides(l->stride);
    init_neighbors(l->axis, l->diff);
}

void lattice_coords(lattice const *l, int i, int *x)
{
    index2coords(i, l->stride, x);
}

void lattice_random_site(lattice const *l, rng *rand, lattice_site *site)
{
    site->i = (int) (LT_N * RNG_RAND_UNIFORM(rand));
    lattice_coords(l, site->i, site->x);
}

int lattice_neighbor(lattice const *l, int n,
                     lattice_site const *restrict site,
                     lattice_site *restrict neighbor)
{
    int axis = l->axis[n],
        diff = l->diff[n][site->x[axis]];

#if (BC != BC_PERIODIC)
    if (diff == 0) return 0;   /* no neighbor */
#endif
    neighbor->i = site->i + diff * l->stride[axis];
    memcpy(neighbor->x, site->x, sizeof site->x);
    neighbor->x[axis] += diff;
    return 1;
}

