#pragma once
#include "utils.h"
#include "sim.h"

void prim_to_cons(struct sim *sim);

void init_sod(struct sim *sim, enum dir dir);
void init_kelvin_helmholtz(struct sim *sim, enum dir dir);
void init_blast(struct sim *sim);

typedef void(*initial_condition_t)(struct sim*, enum dir);

void initial_condition(struct sim *sim, const char* init_name)
{
  const char *avail_init[] = {"sod", "blast", "kelvin-helmholtz"};
  initial_condition_t avail_init_fct[] = {init_sod, (initial_condition_t)init_blast, init_kelvin_helmholtz};

  for(int i=0; i < sizeof(avail_init)/sizeof(char*); i++)
    if(!strncmp(avail_init[i], init_name, 30))
    {
      avail_init_fct[i](sim, IX);
      return;
    }
  printf("unknown initial condition, available: \n\t");
  for(int i=0; i < sizeof(avail_init)/sizeof(char*); i++)
    printf("%s  ", avail_init[i]);
  printf("\n");
  exit(1);
}

void init_blast(struct sim *sim)
{
  DECLARE_PSTATE_VAR
  sim->t = 0;
  sim->boundary[IX] = REFLECTING;
  sim->boundary[IY] = REFLECTING;

  const double r0 = 0.1;

  const real_t *center_x = sim->grid.cellcenter_x;
  const real_t *center_y = sim->grid.cellcenter_y;
  const real_t mid_x = 0.5 * (sim->grid.xmin + sim->grid.xmax);
  const real_t mid_y = 0.5 * (sim->grid.ymin + sim->grid.ymax);

  const real_t thck  = 0.02;
  const real_t p_in  = 1.0;
  const real_t p_out = 0.1;

  for_each_cells_and_ghost_y(i)
    for_each_cells_and_ghost_x(j)
    {
      const u64 id = cell_id(i,j);
      const real_t r = sqrt((center_x[id] - mid_x)*(center_x[id] - mid_x) + (center_y[id] - mid_y)*(center_y[id] - mid_y));
      const real_t tr = 0.5 * (tanh((r - r0) / thck) + 1.0);
      rho[id] = 1.0;
      u[id]   = 0;
      v[id]   = 0;
      p[id]   = p_out * tr + p_in * (1.0-tr);
    }
  prim_to_cons(sim);
}

void init_sod(struct sim *sim, enum dir dir)
{
  DECLARE_PSTATE_VAR
  sim->t = 0;
  sim->boundary[IX] = ABSORBING;
  sim->boundary[IY] = ABSORBING;

  const real_t *cell_center = (dir == IX) ? sim->grid.cellcenter_x : sim->grid.cellcenter_y;
  const real_t middle = (dir == IX) ? 0.5 * (sim->grid.xmin + sim->grid.xmax) : 0.5 * (sim->grid.ymin + sim->grid.ymax);

  for_each_cells_and_ghost_y(i)
    for_each_cells_and_ghost_x(j)
    {
      const u64 id = cell_id(i,j);
      if(cell_center[id] < middle){
        rho[id] = 1.0;
        u[id]   = 0;
        v[id]   = 0;
        p[id]   = 1.0;
      }
      else{
        rho[id] = 0.125;
        u[id]   = 0;
        v[id]   = 0;
        p[id]   = 0.1;
      }
    }
  prim_to_cons(sim);
}

void init_kelvin_helmholtz(struct sim *sim, enum dir dir)
{
  DECLARE_PSTATE_VAR
  sim->t = 0;
  sim->boundary[dir]    = PERIODIC;
  sim->boundary[IY-dir] = PERIODIC;

  const real_t *cell_center = (dir == IY) ? sim->grid.cellcenter_x : sim->grid.cellcenter_y;
  const real_t middle = (dir == IX) ? 0.5 * (sim->grid.xmin + sim->grid.xmax) : 0.5 * (sim->grid.ymin + sim->grid.ymax);

  const real_t pert = 0.02;

  for_each_cells_and_ghost_y(i)
    for_each_cells_and_ghost_x(j)
    {
      const u64 id = cell_id(i,j);

      const real_t tr = sin(M_PI * cell_center[id]/middle);
      // real_t tr = sin(M_PI * cell_center[id]/middle); tr = (fabs(tr) > 0.9) ? tr : 0.0;
      // const real_t top = middle / 2.0;
      // const real_t bot = middle * 3.0 / 2.0;
      // const real_t tr = (fmin(fabs(cell_center[id]-top)/middle, fabs(cell_center[id]-bot)/middle) < 0.1) ? 1.0 : 0.0;

      if(fabs(cell_center[id] - middle) < 0.5 * middle){
        rho[id] = 2.0;
        u[id]   = ((dir == IY) ? pert * tr * ((real_t)rand()/RAND_MAX - 0.5) :  0.5);
        v[id]   = ((dir == IX) ? pert * tr * ((real_t)rand()/RAND_MAX - 0.5) :  0.5);
        p[id]   = 2.5;
      }
      else{
        rho[id] = 1.0;
        u[id]   = ((dir == IY) ? pert * tr * ((real_t)rand()/RAND_MAX - 0.5) : -0.5);
        v[id]   = ((dir == IX) ? pert * tr * ((real_t)rand()/RAND_MAX - 0.5) : -0.5);
        p[id]   = 2.5;
      }
    }
  prim_to_cons(sim);
}

