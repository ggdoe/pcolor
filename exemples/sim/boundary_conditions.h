#pragma once
#include "utils.h"
#include "sim.h"

void boundary_periodic(struct sim *sim, enum dir dir);
void boundary_absorbing(struct sim *sim, enum dir dir);

void fill_boundaries(struct sim *sim)
{
  for(int dir=IX; dir<=IY; dir++)
    switch(sim->boundary[dir])
    {
      case ABSORBING:
        boundary_absorbing(sim, dir);
        break;
      case PERIODIC:
        boundary_periodic(sim, dir);
        break;
    }
}

void boundary_periodic(struct sim *sim, enum dir dir)
{
  DECLARE_CSTATE_VAR
  u64 imaxl, jmaxl, lo;
  u64 iminr, imaxr, jminr, jmaxr, hi;
  
  if(dir == IX){
    imaxl = sim->grid.Ny_tot;
    jmaxl = sim->grid.gx;
    iminr = 0;
    imaxr = sim->grid.Ny_tot;
    jminr = sim->grid.Nx_tot - sim->grid.gx;
    jmaxr = sim->grid.Nx_tot;
    lo    = jmaxl;
    hi    = jminr;
  }
  else{
    imaxl = sim->grid.gy;
    jmaxl = sim->grid.Nx_tot;
    iminr = sim->grid.Ny_tot - sim->grid.gy;
    imaxr = sim->grid.Ny_tot;
    jminr = 0;
    jmaxr = sim->grid.Nx_tot;
    lo    = imaxl;
    hi    = iminr;
  }

  // left side
  #pragma omp parallel for collapse(2)
  for(u64 i=0; i < imaxl; i++)
    for(u64 j=0; j < jmaxl; j++){
      const u64 id    = cell_id(i,j);
      const u64 id_l = (dir == IX) ? cell_id(i,hi-lo+j) : cell_id(hi-lo+i,j);
      rho  [id] = rho  [id_l];
      rho_u[id] = rho_u[id_l];
      rho_v[id] = rho_v[id_l];
      E    [id] = E    [id_l];
    }
  // right side
  #pragma omp parallel for collapse(2)
  for(u64 i=iminr; i < imaxr; i++)
    for(u64 j=jminr; j < jmaxr; j++){
      const u64 id    = cell_id(i,j);
      const u64 id_r = (dir == IX) ? cell_id(i,lo+j) : cell_id(lo+j,j);
      rho  [id] = rho  [id_r];
      rho_u[id] = rho_u[id_r];
      rho_v[id] = rho_v[id_r];
      E    [id] = E    [id_r];
    }
}

void boundary_absorbing(struct sim *sim, enum dir dir)
{
  DECLARE_CSTATE_VAR
  u64 imaxl, jmaxl, lo;
  u64 iminr, imaxr, jminr, jmaxr, hi;
  
  if(dir == IX){
    imaxl = sim->grid.Ny_tot;
    jmaxl = sim->grid.gx;
    iminr = 0;
    imaxr = sim->grid.Ny_tot;
    jminr = sim->grid.Nx_tot - sim->grid.gx;
    jmaxr = sim->grid.Nx_tot;
    lo    = jmaxl;
    hi    = jminr-1;
  }
  else{
    imaxl = sim->grid.gy;
    jmaxl = sim->grid.Nx_tot;
    iminr = sim->grid.Ny_tot - sim->grid.gy;
    imaxr = sim->grid.Ny_tot;
    jminr = 0;
    jmaxr = sim->grid.Nx_tot;
    lo    = imaxl;
    hi    = iminr-1;
  }

  // left side
  #pragma omp parallel for collapse(2)
  for(u64 i=0; i < imaxl; i++)
    for(u64 j=0; j < jmaxl; j++){
      const u64 id    = cell_id(i,j);
      const u64 id_l = (dir == IX) ? cell_id(i,lo) : cell_id(lo,j);
      rho  [id] = rho  [id_l];
      rho_u[id] = rho_u[id_l];
      rho_v[id] = rho_v[id_l];
      E    [id] = E    [id_l];
    }
  // right side
  #pragma omp parallel for collapse(2)
  for(u64 i=iminr; i < imaxr; i++)
    for(u64 j=jminr; j < jmaxr; j++){
      const u64 id    = cell_id(i,j);
      const u64 id_r = (dir == IX) ? cell_id(i,hi) : cell_id(hi,j);
      rho  [id] = rho  [id_r];
      rho_u[id] = rho_u[id_r];
      rho_v[id] = rho_v[id_r];
      E    [id] = E    [id_r];
    }
}
