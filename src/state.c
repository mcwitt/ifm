#include "state.h"
#include <stdlib.h>

void state_alloc(state *s) { s->spin = malloc(LT_N * sizeof *s->spin); }
void state_free(state *s) { free(s->spin); }
