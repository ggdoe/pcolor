#!/bin/sh
CC="gcc"

CFLAGS="-Wall -O3 -g"
LIB="-lm -lSDL2"
IDIR="../."

# plot.c
SRC1="./plot.c"
out1="plot.out"
# pcolor_nostate.c
SRC2="./pcolor_nostate.c"
out2="pcolor_nostate.out"
# pcolor.c
SRC3="./pcolor.c"
out3="pcolor.out"
# animate.c
SRC4="./animate.c"
out4="animate.out"

set -xe
$CC $CFLAGS -I$IDIR -o $out1 $SRC1 $LIB
$CC $CFLAGS -I$IDIR -o $out2 $SRC2 $LIB
$CC $CFLAGS -I$IDIR -o $out3 $SRC3 $LIB
$CC $CFLAGS -I$IDIR -o $out4 $SRC4 $LIB


