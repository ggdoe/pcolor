#include "sim1d.h"
#include "show.h"
#include "plot.h"
#include "pcolor.h"

#define IMG_FACTOR 64
#define IMG_WIDTH  (16 * IMG_FACTOR)
#define IMG_HEIGHT (12  * IMG_FACTOR)
#define NB_POINTS 256

int main(int argc, char** argv)
{
  uint32_t *pixels = malloc(IMG_WIDTH * IMG_HEIGHT * sizeof(uint32_t));
  fill(pixels, IMG_WIDTH, IMG_HEIGHT, 255, 255, 255);

  struct sim sim = init_sim(NB_POINTS);
  init_state(&sim);
  run(&sim, 0.2);
  
  const real_t *rho = sim.pstate.rho + sim.grid.gx;
  const real_t *  u = sim.pstate.u   + sim.grid.gx;
  const real_t *  p = sim.pstate.p   + sim.grid.gx;

  struct lim xlim = {0.0, 1.0};
  struct lim ylim = compute_lim(rho, sim.grid.Nx, NULL);
             ylim = compute_lim(  u, sim.grid.Nx, &ylim);
             ylim = compute_lim(  p, sim.grid.Nx, &ylim);
  ylim.min -= 0.02;
  ylim.max += 0.02;

  plot(pixels, IMG_WIDTH, IMG_HEIGHT, rho,  sim.grid.Nx, &ylim, 0xFFFF0000);
  plot(pixels, IMG_WIDTH, IMG_HEIGHT,   u,  sim.grid.Nx, &ylim, 0xFF00FF00);
  plot(pixels, IMG_WIDTH, IMG_HEIGHT,   p,  sim.grid.Nx, &ylim, 0xFF0000FF);
  fill_grid(pixels, IMG_WIDTH, IMG_HEIGHT, &xlim, &ylim, 0xFFAAAAAA);
  show(pixels, IMG_WIDTH, IMG_HEIGHT);

  free_sim(&sim);
  return 0;
}
