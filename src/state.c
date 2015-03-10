#include "state.h"
#include <stdlib.h>

int state_alloc(state *s)
{
    if ((s->spin = malloc(LT_N * sizeof *s->spin)))
        return 0;

    return 1;
}

void state_free(state *s) { free(s->spin); }
