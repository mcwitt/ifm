#!/bin/bash

# Script to generate source files from templates ending with `.in`.
# Uses `sed` to replace the pattern `@variable@` with variable's value.

srcdir="$(dirname "$0")/src"
cmd=sed

# Set parameter value (appends to sed command)
function setp {
    eval $1=$2
    cmd="$cmd -e s|@$1@|$2|"
}

# Use parameter value from environment or default if not set.
# (appends to sed command)
function getp {
    eval $1=\${$1:-$2}
    val=$(eval echo \$$1)
    cmd="$cmd -e s|@$1@|$val|"
}

# PARAMETERS AND DEFAULT VALUES

getp srcdir $srcdir
getp BC BC_PERIODIC
getp LT_LMAX 64
getp LT_D 2
getp LOG_BOOTSTRAP_STEPS 3
getp LOG_MEAS_TIME 8
getp UPDATES_PER_STEP 1
getp MAX_STACK_SIZE $((LT_LMAX**LT_D))
getp RNG RNG_DSFMT

# DERIVED CONSTANTS

setp LT_N $((LT_LMAX**LT_D))

for f in $srcdir/*.in; do
    $cmd $f > $(basename "${f%.in}")
done

