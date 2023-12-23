#pragma once
#include "utils.h"
#include "sim.h"

void prim_to_cons(struct sim *sim);

void init_sod(struct sim *sim, enum dir dir);
void init_kelvin_helmholtz(struct sim *sim, enum dir dir);

void init_state(struct sim *sim)
{
  sim->t = 0;
  
  // init_sod(sim, IY);
  init_kelvin_helmholtz(sim, IX);

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
  sim->boundary[IX] = ABSORBING;
  sim->boundary[IY] = PERIODIC;

  const real_t *cell_center = (dir == IY) ? sim->grid.cellcenter_x : sim->grid.cellcenter_y;
  const real_t middle = (dir == IX) ? 0.5 * (sim->grid.xmin + sim->grid.xmax) : 0.5 * (sim->grid.ymin + sim->grid.ymax);

  const real_t pert = 0.04;

  for_each_cells_and_ghost_y(i)
    for_each_cells_and_ghost_x(j)
    {
      const u64 id = cell_id(i,j);
      if(cell_center[id] < middle){
        rho[id] = 2.0;
        u[id]   = ((dir == IY) ? pert * ((real_t)rand()/RAND_MAX - 0.5) :  1.0);
        v[id]   = ((dir == IX) ? pert * ((real_t)rand()/RAND_MAX - 0.5) :  1.0);
        p[id]   = 1.0;
      }
      else{
        rho[id] = 1.0;
        u[id]   = ((dir == IY) ? pert * ((real_t)rand()/RAND_MAX - 0.5) : -1.0);
        v[id]   = ((dir == IX) ? pert * ((real_t)rand()/RAND_MAX - 0.5) : -1.0);
        p[id]   = 1.0;
      }
    }
  prim_to_cons(sim);
}

