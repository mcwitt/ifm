#include "state.h"

void state_alloc(state *s)
{
    s->spin = malloc(LT_N * sizeof *s->spin);
    s->stack = malloc(MAX_STACK_SIZE * sizeof *s->stack);
}

void state_free(state *s)
{
    free(s->spin);
    free(s->stack);
}
