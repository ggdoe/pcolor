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
  const char *initial_condition_name;
};

// do nothing
void callback_dummy(void *args){}

// compute minmax for the cmap at the current timestep
void callback_minmax(void *args)
{
  struct callback_args *s = args;
  compute_minmax(&s->min, &s->max, &s->sim->grid, s->state_to_draw);
  printf("\tminmax: %lf, %lf\n", s->min, s->max);
  pcolor(s->pcolor_config, s->state_to_draw+s->offset, s->min, s->max, s->cmap);
}

// cycle state to draw (rho, u, v, p, E)
void callback_cycle_drawstate(void *args)
{
  struct callback_args *s = args;
  static u64 i=0;
  printf("\x1b[2K"); // clear line
  switch(i++)
  {
    case 0:
      printf("\rstate rho");
      s->state_to_draw = s->sim->pstate.rho;
      break;
    case 1:
      printf("\rstate u  ");
      s->state_to_draw = s->sim->pstate.u;
      break;
    case 2:
      printf("\rstate v  ");
      s->state_to_draw = s->sim->pstate.v;
      break;
    case 3:
      printf("\rstate p  ");
      s->state_to_draw = s->sim->pstate.p;
      break;
    case 4:
      printf("\rstate E  ");
      s->state_to_draw = s->sim->cstate.E;
      break;
  }
  callback_minmax(args);
  i = i%5;
}

// callback that is call at each frame
void callback_update(void *args)
{
  struct callback_args *s = args;
  run(s->sim, s->sim->t + 0.01);
  printf("\rt: %.4lf", s->sim->t);fflush(stdout);
  pcolor(s->pcolor_config, s->state_to_draw+s->offset, s->min, s->max, s->cmap);

  // slice_x on middle of the y axis
  // plot(s->pcolor_config->pixels, IMG_WIDTH, IMG_HEIGHT, s->state_to_draw+s->sim->grid.Ny_tot*s->sim->grid.Nx_tot/2 + s->sim->grid.gx, s->sim->grid.Nx, NULL, 0xFFFFFFFF);
}

// activate or deactivate edges
void callback_cycle_edge(void *args)
{
  struct callback_args *s = args;
  s->pcolor_config->show_edge = !s->pcolor_config->show_edge;
  pcolor(s->pcolor_config, s->state_to_draw+s->offset, s->min, s->max, s->cmap);
}

// restart the simulation
void callback_init(void *args)
{
  struct callback_args *s = args;
  initial_condition(s->sim, s->initial_condition_name);
  callback_minmax(args);
}

// --- MAIN ---

int main(int argc, char** argv)
{
  // create pixel array
  uint32_t *pixels = MALLOC(IMG_WIDTH * IMG_HEIGHT * sizeof(uint32_t));

  // init fluid simulation on cartesian grid
  struct sim sim = init_sim(NB_POINTS, NB_POINTS);
  const char* initial_condition_name = "kelvin-helmholtz";
  initial_condition(&sim, initial_condition_name);

  // // offset to avoid ghost cells
  const u64 offset =  sim.grid.gy*sim.grid.Nx_tot + sim.grid.gx;

  // // create pcolor config
  struct pcolor_config pcolor_config = pcolor_config_alloc(pixels, IMG_WIDTH, IMG_HEIGHT);
  pcolor_config.show_edge = false;
  pcolor_config_fill(&pcolor_config, sim.grid.vertex_x+offset, sim.grid.vertex_y+offset, sim.grid.Nx, sim.grid.Ny, sim.grid.Nx_tot);

  // callback functions arguments
  struct callback_args callback_args = {.state_to_draw=sim.pstate.rho, .pcolor_config = &pcolor_config, 
                                        .sim=&sim, .offset=offset, .cmap=cmap_cool_warm, .initial_condition_name=initial_condition_name};

  // custom key event
  struct custom_keyevent kevent[5]; const int nb_custom_event = sizeof(kevent)/sizeof(struct custom_keyevent);
  kevent[0] = (struct custom_keyevent){.key=SDLK_SPACE, .callback=callback_update,          .callback_args=&callback_args};
  kevent[1] = (struct custom_keyevent){.key=SDLK_s,     .callback=callback_cycle_drawstate, .callback_args=&callback_args};
  kevent[2] = (struct custom_keyevent){.key=SDLK_i,     .callback=callback_init,            .callback_args=&callback_args};
  kevent[3] = (struct custom_keyevent){.key=SDLK_m,     .callback=callback_minmax,          .callback_args=&callback_args};
  kevent[4] = (struct custom_keyevent){.key=SDLK_e,     .callback=callback_cycle_edge,      .callback_args=&callback_args};

  // compute minmax for the cmap function
  callback_cycle_drawstate(&callback_args);

  // start window and animate the simulation
  animate(pixels, IMG_WIDTH, IMG_HEIGHT, 60.0, callback_update, &callback_args, kevent, nb_custom_event);
  // animate(pixels, IMG_WIDTH, IMG_HEIGHT, 30.0, callback_dummy, &callback_args, kevent, nb_custom_event);

  // freeing everything
  free(pixels);
  free_sim(&sim);
  pcolor_free(&pcolor_config);
  return 0;
}
