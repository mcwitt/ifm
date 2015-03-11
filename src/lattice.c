#include "lattice.h"

#if (BC == BC_PERIODIC)
#define DELTA(d) (shape[d] - 1)
#else
#define DELTA(d) 0
#endif

static void init_strides(int const *shape, site_label *stride)
{
    int d;

    stride[LT_D-1] = 1;

    for (d = LT_D-1; d > 0; d--)
        stride[d-1] = shape[d] * stride[d];
}

static void init_neighbors(int const *shape, int diff[][LT_LMAX])
{
    int d, l;

    /* forward neighbors */

    for (d = 0; d < LT_D; d++)
    {
        for (l = 0; l < shape[d]-1; l++) diff[d][l] = 1;
        diff[d][shape[d]-1] = -DELTA(d);
    }

    /* backward neighbors */

    for (d = 0; d < LT_D; d++)
    {
        int i = d + LT_D;
        diff[i][0] = DELTA(d);
        for (l = 1; l < shape[d]; l++) diff[i][l] = -1;
    }
}

void lattice_init(lattice *l, int *shape)
{
    init_strides(shape, l->stride);
    init_neighbors(shape, l->diff);
}

void neighbor_iter_init(site_label site, neighbor_iter *iter)
{
    iter->i = 0;
    iter->site = site;
    iter->r = site;
}

site_label neighbor_iter_next(lattice const *l, neighbor_iter *iter)
{
    int i, dx;

    while (iter->i < LT_D)
    {
        i = iter->i++;
        iter->x[i] = iter->r / l->stride[i];
        iter->r %= l->stride[i];
        dx = l->diff[i][iter->x[i]];
        if (dx != 0) return iter->site + dx * l->stride[i];
    }

    while (iter->i < LT_Z)
    {
        i = iter->i++;
        dx = l->diff[i][iter->x[i-LT_D]];
        if (dx != 0) return iter->site + dx * l->stride[i-LT_D];
    }

    return -1;
}

