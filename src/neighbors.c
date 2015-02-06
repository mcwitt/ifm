#include "neighbors.h"

#if (BC == BC_PERIODIC)
#define DELTA(d) (lt_shape[d] - 1)
#elif (BC == BC_FREE)
#define DELTA(d) 0
#else
#error invalid boundary condition
#endif

void init_neighbors(int *axis, int diff[][LT_LMAX])
{
    int d, l, end_diff;

    /* forward neighbors */

    for (d = 0; d < LT_D; d++)
    {
        axis[d] = d;
        for (l = 0; l < lt_shape[d]-1; l++) diff[d][l] = 1;
        diff[d][lt_shape[d]-1] = -DELTA(d);
    }

    /* backward neighbors */

    for (d = 0; d < LT_D; d++)
    {
        int n = d + LT_D;
        axis[n] = d;
        diff[n][0] = DELTA(d);
        for (l = 1; l < lt_shape[d]; l++) diff[n][l] = -1;
    }
}

