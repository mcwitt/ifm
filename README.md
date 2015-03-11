ifm: Monte Carlo simulation of Ising ferromagnets
-------------------------------------------------------------------------------

Basic usage
-------------------------------------------------------------------------------

Simulation parameters are set before compilation to allow for the most
effective optimization. Use the script `setup.sh` to prepare a build for a
given set of parameters. Binaries with different sets of parameters should be
built in separate directories.

For example, to build a simulation of a 2-d model on a square lattice with
periodic boundary conditions and linear size 64:

```bash
mkdir -p build/d2-L64
cd build/d2-L64
BC=BC_PERIODIC LT_D=2 LT_LMAX=64 ../../setup.sh
make
```

The `ifm` binary takes a single optional argument, which is used to seed the
random number generator. Temperatures are read from standard input. Typical
usage looks like

```bash
./ifm 123 < temps.txt > d2-L64-123.txt
```

Parameters
-------------------------------------------------------------------------------

The following environment variables are understood by `setup.sh`:

__Basic parameters__

 Parameter          | Type                      | Description
 ------------------ | ------------------------- | -----------
 `BC`               | {`BC_PERIODIC`,`BC_FREE`} | boundary conditions
 `LT_D`             | integer                   | dimension of lattice
 `LT_LMAX`          | integer                   | linear size of lattice
 `LOG_TIME_THERM`   | integer                   | log time to thermalize
 `LOG_TIME_MEAS`    | integer                   | log time to measure

__Advanced parameters__ (usually don't need adjustment)

 Parameter          | Type                      | Description
 ------------------ | ------------------------- | -----------
 `RNG`              | [see `rng.h.in` for list] | type of random number generator to use
 `UPDATES_PER_STEP` | integer                   | number of cluster updates per unit time
 `MAX_STACK_SIZE`   | integer                   | size of the stack used for cluster updates

