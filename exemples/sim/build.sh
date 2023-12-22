#!/bin/sh
CC="gcc"

CFLAGS="-Wall -O3 -g -fopenmp"
LIB="-lm -lSDL2"
IDIR="../../."

# sim.c
SRC1="./sim.c"
out1="sim.out"
# sim_slice1d.c
SRC2="./sim_slice1d.c"
out2="slice1d.out"

set -xe
$CC $CFLAGS -I$IDIR -o $out1 $SRC1 $LIB
# $CC $CFLAGS -I$IDIR -o $out2 $SRC2 $LIB

