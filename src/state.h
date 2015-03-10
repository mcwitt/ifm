#ifndef STATE_H
#define STATE_H

#include "lattice.h"
#include "spin.h"

typedef struct
{
    ising_spin *spin;
    int64_t magnetization;
}
state;

int state_alloc(state *s);
void state_free(state *s);

#endif
