#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <time.h>
#include <math.h>

#include "show.h"
#include "plot.h"
#include "cmap.h"
#include "utils.h"

#include "pcolor.h"

#define IMG_FACTOR 60
#define IMG_WIDTH  (16 * IMG_FACTOR)
#define IMG_HEIGHT (12  * IMG_FACTOR)

int main(int argc, char **argv)
{
  uint32_t *pixels = malloc(IMG_WIDTH * IMG_HEIGHT * sizeof(uint32_t));
  memset(pixels, 0xFFFFFFFF, IMG_WIDTH * IMG_HEIGHT * sizeof(uint32_t));
  
  // setup grid
  int w = 88;
  int h = 19;
  double x[w * h];
  double y[w * h];
  uint32_t C[(w-1) * (h-1)];
  
  meshgrid(x, y, w, h, 0.0, 1.0, 0.0, 1.0);

  // remap grid
  ring_map(x, y, w, h);
  colella_map(x, y, w, h);

  for (int i = 0; i < (w-1) * (h-1); i++){
      C[i] = cmap_nipy_spectral((double)rand()/RAND_MAX);
  }

  struct pcolor_state pstate = pcolor_state_alloc(pixels + 100, IMG_HEIGHT, IMG_HEIGHT);

  // shifting with leading dimension parameter
  pstate.ld_pixels = IMG_WIDTH;
  
  pcolor_fill_state(&pstate, x, y, w, h);
  pstate.show_edge = false;

  pcolor(&pstate, C);
  
  // pcolor_nostate(pixels, IMG_WIDTH, IMG_HEIGHT, x, y, w, h, C, true, 0xFFFFFFFF, 0xFF111111);

  show(pixels, IMG_WIDTH, IMG_HEIGHT);

  free(pixels);
  pcolor_free(&pstate);
  return 0;
}
