#ifndef STATE_H
#define STATE_H

#include "lattice.h"
#include "spin.h"

#define MAX_STACK_SIZE LT_N     /* conservative value */

typedef struct
{
    ising_spin *spin;
    lattice_site *stack;
    uint_ord_N cluster_size;
    ising_spin cluster_spin;
    int_ord_N magnetization;
}
state;

void state_alloc(state *s);

void state_free(state *s);

#endif
