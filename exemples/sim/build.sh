#!/bin/sh
CC="gcc"

CFLAGS="-Wall -O3 -g"
LIB="-lm -lSDL2"
IDIR="../../."

# animate.c
SRC="./sim.c"
out="sim.out"

set -xe
$CC $CFLAGS -I$IDIR -o $out $SRC $LIB -Wno-unknown-pragmas

