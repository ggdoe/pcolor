#!/bin/sh

SRC="./main.c"
out_name="a.out"

CFLAGS="-Wall -O3"
LFLAGS="-fopenmp -flto"
LIB="-lm -lSDL2"

CC="gcc"

set -xe
$CC $CFLAGS $LFLAGS -o $out_name $SRC $LIB


