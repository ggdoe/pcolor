#!/bin/sh
CC="gcc"

CFLAGS="-Wall -O3 -g -fopenmp"
LIB="-lm -lSDL2"
IDIR="../../."

# animate.c
SRC="./sod.c"
out="sod.out"

set -xe
$CC $CFLAGS -I$IDIR -o $out $SRC $LIB

