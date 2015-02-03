#ifndef SPIN_H
#define SPIN_H

typedef int ising_spin;

/* define representation of Ising spins */

#define SPIN_UP 1
#define SPIN_DN -1
#define SPIN_MOMENT(s) (s)
#define SPIN_FLIP(s) ((s) *= -1)

#endif
