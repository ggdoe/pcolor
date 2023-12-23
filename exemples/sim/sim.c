#include "sim.h"
#include "show.h"
#include "plot.h"
#include "pcolor.h"
#include "cmap.h"

#define IMG_FACTOR 64
#define IMG_WIDTH  (10 * IMG_FACTOR)
#define IMG_HEIGHT (10  * IMG_FACTOR)
#define NB_POINTS 10

void compute_minmax(real_t *out_min, real_t *out_max, struct grid *grid, real_t *values)
{
  real_t min =  FLT_MAX;
  real_t max = -FLT_MAX;

  #pragma omp parallel for collapse(2) reduction(min:min) reduction(max:max)
  for(u64 i=grid->imin  ; i < grid->imax  ; i++){
    for(u64 j=grid->jmin  ; j < grid->jmax  ; j++)
    {
      const u64 id = (i*grid->Nx_tot + j);
      min = (min < values[id]) ? min : values[id];
      max = (max > values[id]) ? max : values[id];
    }
  }

  *out_min = min;
  *out_max = max;
}

void fill_pixels(uint32_t *pixels, struct grid *grid, real_t *values)
{
  real_t min = 0.0, max = 1.0;
  compute_minmax(&min, &max, grid, values);

  #pragma omp parallel for collapse(2)
  for(u64 i=grid->imin  ; i < grid->imax  ; i++)
    for(u64 j=grid->jmin  ; j < grid->jmax  ; j++){
      const u64 id  = (i*grid->Nx_tot + j);
      const u64 id_pixel = ((i-grid->imin)*grid->Nx + (j-grid->jmin));
      pixels[id_pixel] = cmap_nipy_spectral((values[id]-min)/(max-min));
    }
}

struct callback_args{
  struct pcolor_state *pcolor_state;
  struct sim *sim;
  u64 offset;
};

void callback_dummy(void *args){}

void callback_update(void *args)
{
  struct callback_args *s = args;
  run(s->sim, s->sim->t + 0.005);
  real_t min, max;
  compute_minmax(&min, &max, &s->sim->grid, s->sim->pstate.rho);
  // pcolor_real(s->pcolor_state, s->sim->pstate.rho+s->offset, 0.5, 2.6);
  pcolor_real(s->pcolor_state, s->sim->pstate.rho+s->offset, min, max);
}

int main(int argc, char** argv)
{
  uint32_t *pixels = malloc(IMG_WIDTH * IMG_HEIGHT * sizeof(uint32_t));
  fill(pixels, IMG_WIDTH, IMG_HEIGHT, 255, 255, 255);

  struct sim sim = init_sim(NB_POINTS, NB_POINTS);
  init_state(&sim);

  struct pcolor_state pcolor_state = pcolor_state_alloc(pixels, IMG_WIDTH, IMG_HEIGHT);
  // pcolor_state.show_edge = false;

  const u64 offset = 0; sim.grid.gy*sim.grid.Nx_tot + sim.grid.gx;
  struct callback_args callback_args = {.pcolor_state = &pcolor_state, .sim=&sim, .offset=offset};
  struct custom_keyevent kevent = {.key=SDLK_SPACE, .callback=callback_update, .callback_args=&callback_args};
  pcolor_fill_state(&pcolor_state, sim.grid.vertex_x+offset, sim.grid.vertex_y+offset, sim.grid.Nx_tot, sim.grid.Ny_tot, sim.grid.Nx_tot);

  real_t min, max;
  compute_minmax(&min, &max, &sim.grid, sim.pstate.rho);
  pcolor_real(&pcolor_state, sim.pstate.rho+offset, min, max);
  
  animate(pixels, IMG_WIDTH, IMG_HEIGHT, 0.0, callback_dummy, &callback_args, &kevent, 1);

  free_sim(&sim);
  return 0;
}
