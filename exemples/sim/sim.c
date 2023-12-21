#include "sim.h"
#include "show.h"
#include "plot.h"
#include "pcolor.h"

#define IMG_FACTOR 60
#define IMG_WIDTH  (16 * IMG_FACTOR)
#define IMG_HEIGHT (12  * IMG_FACTOR)
#define NB_POINTS 128

// void print_state(struct sim *sim){
//   for(u32 i=sim->grid.imin ; i < sim->grid.imax; i++)
//     printf("%.3f\t%.3f\t%.3f\t%.3f\n", sim->grid.cellcenter[i], sim->pstate.rho[i], sim->pstate.u[i], sim->pstate.p[i]);
// }

int main(int argc, char** argv)
{
  uint32_t *pixels = malloc(IMG_WIDTH * IMG_HEIGHT * sizeof(uint32_t));
  fill(pixels, IMG_WIDTH, IMG_HEIGHT, 255, 255, 255);


  struct sim sim = init_sim(NB_POINTS, NB_POINTS);
  init_state(&sim);
  // run(&sim, 0.2);

  free_sim(&sim);
  return 0;
}
