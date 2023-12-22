#include "sim.h"
#include "show.h"
#include "plot.h"
#include "pcolor.h"
#include "cmap.h"

#define IMG_FACTOR 64
#define IMG_WIDTH  (16 * IMG_FACTOR)
#define IMG_HEIGHT (12  * IMG_FACTOR)
#define NB_POINTS_X 256
#define NB_POINTS_Y 32

int main(int argc, char** argv)
{
  uint32_t *pixels = malloc(IMG_WIDTH * IMG_HEIGHT * sizeof(uint32_t));
  fill(pixels, IMG_WIDTH, IMG_HEIGHT, 255, 255, 255);

  struct sim sim = init_sim(NB_POINTS_X, NB_POINTS_Y);
  init_state(&sim);

  run(&sim, 0.2);

  // plot 1d
  const real_t *rho_slice = sim.pstate.rho + sim.grid.Nx_tot*sim.grid.Ny_tot/2 + sim.grid.gx;
  const real_t *  u_slice = sim.pstate.u   + sim.grid.Nx_tot*sim.grid.Ny_tot/2 + sim.grid.gx;
  const real_t *  v_slice = sim.pstate.v   + sim.grid.Nx_tot*sim.grid.Ny_tot/2 + sim.grid.gx;
  const real_t *  p_slice = sim.pstate.p   + sim.grid.Nx_tot*sim.grid.Ny_tot/2 + sim.grid.gx;
  
  struct lim ylim = compute_lim(rho_slice, sim.grid.Nx, NULL);
             ylim = compute_lim(  u_slice, sim.grid.Nx, &ylim);
             ylim = compute_lim(  v_slice, sim.grid.Nx, &ylim);
             ylim = compute_lim(  p_slice, sim.grid.Nx, &ylim);
  ylim.min -= 0.02;
  ylim.max += 0.02;
  plot(pixels, IMG_WIDTH, IMG_HEIGHT, rho_slice, sim.grid.Nx, &ylim, 0xFFFF0000);
  plot(pixels, IMG_WIDTH, IMG_HEIGHT,   u_slice, sim.grid.Nx, &ylim, 0xFF00FF00);
  plot(pixels, IMG_WIDTH, IMG_HEIGHT,   v_slice, sim.grid.Nx, &ylim, 0xFF00FFFF);
  plot(pixels, IMG_WIDTH, IMG_HEIGHT,   p_slice, sim.grid.Nx, &ylim, 0xFF0000FF);
  struct lim xlim = {0.0, 1.0};
  fill_grid(pixels, IMG_WIDTH, IMG_HEIGHT, &xlim, &ylim, 0xFFAAAAAA);
  show(pixels, IMG_WIDTH, IMG_HEIGHT);

  free_sim(&sim);
  return 0;
}
