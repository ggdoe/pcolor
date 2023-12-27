#include "sim.h"
#include "show.h"
#include "plot.h"
#include "pcolor.h"
#include "cmap.h"

// set initial window size
#define IMG_FACTOR 64
#define IMG_WIDTH  (10 * IMG_FACTOR)
#define IMG_HEIGHT (10  * IMG_FACTOR)

// set number cell of the simulation in each direction
#define NB_POINTS 128


// compute minmax values that are not ghost cells
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

// arguments common to all callback
struct callback_args{
  struct pcolor_config *pcolor_config;
  struct sim *sim;
  real_t *state_to_draw;
  u64 offset;
  real_t min, max;
  cmap_function_t cmap;
};

// do nothing
void callback_dummy(void *args){}

// compute minmax for the cmap at the current timestep
void callback_minmax(void *args)
{
  struct callback_args *s = args;
  compute_minmax(&s->min, &s->max, &s->sim->grid, s->state_to_draw);
  printf("t: %lf\tmin: %lf\t%lf\n", s->sim->t, s->min, s->max);
  pcolor_real(s->pcolor_config, s->state_to_draw+s->offset, s->min, s->max, s->cmap);
}

// cycle state to draw (rho, p, u, v, E)
void callback_cycle_drawstate(void *args)
{
  struct callback_args *s = args;
  static u64 i=0;
  switch(++i%5)
  {
    case 0:
      printf("state rho\n");
      s->state_to_draw = s->sim->pstate.rho;
      break;
    case 1:
      printf("state p\n");
      s->state_to_draw = s->sim->pstate.p;
      break;
    case 2:
      printf("state u\n");
      s->state_to_draw = s->sim->pstate.u;
      break;
    case 3:
      printf("state v\n");
      s->state_to_draw = s->sim->pstate.v;
      break;
    case 4:
      printf("state E\n");
      s->state_to_draw = s->sim->cstate.E;
      break;
  }
  callback_minmax(args);
}

// callback that is call at each frame
void callback_update(void *args)
{
  struct callback_args *s = args;
  run(s->sim, s->sim->t + 0.01);
  printf("t: %lf\n", s->sim->t);
  pcolor_real(s->pcolor_config, s->state_to_draw+s->offset, s->min, s->max, s->cmap);
}

// restart the simulation
void callback_init(void *args)
{
  struct callback_args *s = args;
  init_state(s->sim);
  callback_minmax(args);
}

// --- MAIN ---

int main(int argc, char** argv)
{
  // create pixel array
  uint32_t *pixels = malloc(IMG_WIDTH * IMG_HEIGHT * sizeof(uint32_t));

  // init fluid simulation on cartesian grid
  struct sim sim = init_sim(NB_POINTS, NB_POINTS);
  init_state(&sim);

  // // offset to avoid ghost cells
  const u64 offset =  sim.grid.gy*sim.grid.Nx_tot + sim.grid.gx;

  // // create pcolor config
  struct pcolor_config pcolor_config = pcolor_config_alloc(pixels, IMG_WIDTH, IMG_HEIGHT);
  pcolor_config.show_edge = false;
  pcolor_config_fill(&pcolor_config, sim.grid.vertex_x+offset, sim.grid.vertex_y+offset, sim.grid.Nx, sim.grid.Ny, sim.grid.Nx_tot);

  // callback functions arguments
  struct callback_args callback_args = {.state_to_draw=sim.pstate.rho, .pcolor_config = &pcolor_config, .sim=&sim, .offset=offset, .cmap=cmap_cool_warm};

  // custom key event
  struct custom_keyevent kevent[4]; const int nb_custom_event = sizeof(kevent)/sizeof(struct custom_keyevent);
  kevent[0] = (struct custom_keyevent){.key=SDLK_SPACE, .callback=callback_update,          .callback_args=&callback_args};
  kevent[1] = (struct custom_keyevent){.key=SDLK_s,     .callback=callback_cycle_drawstate, .callback_args=&callback_args};
  kevent[2] = (struct custom_keyevent){.key=SDLK_r,     .callback=callback_init,            .callback_args=&callback_args};
  kevent[3] = (struct custom_keyevent){.key=SDLK_m,     .callback=callback_minmax,          .callback_args=&callback_args};

  // compute minmax for the cmap function
  callback_minmax(&callback_args);

  // start window and animate the simulation
  animate(pixels, IMG_WIDTH, IMG_HEIGHT, 60.0, callback_update, &callback_args, kevent, nb_custom_event);
  // animate(pixels, IMG_WIDTH, IMG_HEIGHT, 30.0, callback_dummy, &callback_args, kevent, nb_custom_event);

  // freeing everything
  free(pixels);
  free_sim(&sim);
  pcolor_free(&pcolor_config);
  return 0;
}
