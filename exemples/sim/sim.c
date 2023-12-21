#include "sim.h"
#include "show.h"
#include "plot.h"
#include "pcolor.h"

#define IMG_FACTOR 60
#define IMG_WIDTH  (16 * IMG_FACTOR)
#define IMG_HEIGHT (12  * IMG_FACTOR)
#define NB_POINTS 256

// void print_state(struct sim *sim){
//   for(u32 i=sim->grid.imin ; i < sim->grid.imax; i++)
//     printf("%.3f\t%.3f\t%.3f\t%.3f\n", sim->grid.cellcenter[i], sim->pstate.rho[i], sim->pstate.u[i], sim->pstate.p[i]);
// }

int main(int argc, char** argv)
{
  uint32_t *pixels = malloc(IMG_WIDTH * IMG_HEIGHT * sizeof(uint32_t));
  fill(pixels, IMG_WIDTH, IMG_HEIGHT, 255, 255, 255);


  struct sim sim = init_sim(NB_POINTS, 64);
  init_state(&sim);
  run(&sim, 0.2);
  
  const real_t *rho_slice = sim.pstate.rho + sim.grid.Nx_tot * sim.grid.Ny_tot / 2;
  // const real_t *  u_slice = sim.pstate.u   + sim.grid.Nx_tot * sim.grid.Ny_tot / 2;
  // const real_t *  v_slice = sim.pstate.v   + sim.grid.Nx_tot * sim.grid.Ny_tot / 2;
  // const real_t *  p_slice = sim.pstate.p   + sim.grid.Nx_tot * sim.grid.Ny_tot / 2;
  
  // plot
  struct lim xlim = {0.0, 1.0};
  struct lim ylim = compute_lim(rho_slice, NB_POINTS, NULL);
            //  ylim = compute_lim(  u_slice, NB_POINTS, &ylim);
            //  ylim = compute_lim(  v_slice, NB_POINTS, &ylim);
            //  ylim = compute_lim(  p_slice, NB_POINTS, &ylim);
  ylim.min -= 0.2;
  ylim.max += 0.2;
  fill_grid(pixels, IMG_WIDTH, IMG_HEIGHT, &xlim, &ylim, 0xFFAAAAAA);
  plot(pixels, IMG_WIDTH, IMG_HEIGHT, rho_slice, NB_POINTS, &ylim, 0xFFFF0000);
  // plot(pixels, IMG_WIDTH, IMG_HEIGHT,   u_slice, NB_POINTS, &ylim, 0xFF00FF00);
  // plot(pixels, IMG_WIDTH, IMG_HEIGHT,   v_slice, NB_POINTS, &ylim, 0xFF0000FF);
  // plot(pixels, IMG_WIDTH, IMG_HEIGHT,   p_slice, NB_POINTS, &ylim, 0xFF00FFFF);
  show(pixels, IMG_WIDTH, IMG_HEIGHT);

  free_sim(&sim);
  return 0;
}
