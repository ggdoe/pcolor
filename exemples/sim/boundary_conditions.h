#pragma once
#include "utils.h"
#include "sim.h"

void boundary_reflecting(struct sim *sim, enum dir dir);
void boundary_periodic(struct sim *sim, enum dir dir);
void boundary_periodic_v2(struct sim *sim, enum dir dir);
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
      case REFLECTING:
        boundary_reflecting(sim, dir);
        break;
    }
}


void boundary_reflecting(struct sim *sim, enum dir dir)
{
  DECLARE_CSTATE_VAR
  
  if(dir == IX){
    const int64_t jmax = sim->grid.jmax;
    const int64_t gx   = sim->grid.gx;

    #pragma omp parallel for collapse(2)
    for(int64_t i=0; i < sim->grid.Ny_tot; i++)
      for(int64_t j=0; j < gx; j++){
        const int64_t id_l  = cell_id(i,j);
        const int64_t sym_l = cell_id(i,2*gx-j-1);
        const int64_t id_r  = cell_id(i,jmax+j);
        const int64_t sym_r = cell_id(i,jmax-j-1);

        rho  [id_l] = rho  [sym_l];
        rho_u[id_l] = -rho_u[sym_l];
        rho_v[id_l] = rho_v[sym_l];
        E    [id_l] = E    [sym_l];

        rho  [id_r] = rho  [sym_r];
        rho_u[id_r] = -rho_u[sym_r];
        rho_v[id_r] = rho_v[sym_r];
        E    [id_r] = E    [sym_r];
      }
  }
  else{
    const int64_t imax = sim->grid.imax;
    const int64_t gy   = sim->grid.gy;
    
    #pragma omp parallel for collapse(2)
    for(int64_t i=0; i < gy; i++)
      for(int64_t j=0; j < sim->grid.Nx_tot; j++){
        const int64_t id_l  = cell_id(i,j);
        const int64_t sym_l = cell_id(2*gy-i-1,j);
        const int64_t id_r  = cell_id(imax+i,j);
        const int64_t sym_r = cell_id(imax-i-1,j);

        rho  [id_l] = rho  [sym_l];
        rho_u[id_l] = rho_u[sym_l];
        rho_v[id_l] = -rho_v[sym_l];
        E    [id_l] = E    [sym_l];

        rho  [id_r] = rho  [sym_r];
        rho_u[id_r] = rho_u[sym_r];
        rho_v[id_r] = -rho_v[sym_r];
        E    [id_r] = E    [sym_r];
      }
  }
}

void boundary_periodic(struct sim *sim, enum dir dir)
{
  DECLARE_CSTATE_VAR
  
  if(dir == IX){
    const int64_t jmax = sim->grid.jmax;
    const int64_t gx   = sim->grid.gx;
    
    #pragma omp parallel for collapse(2)
    for(int64_t i=0; i < sim->grid.Ny_tot; i++)
      for(int64_t j=0; j < gx; j++){
        const int64_t id_l  = cell_id(i,j);
        const int64_t sym_l = cell_id(i,jmax-gx+j);
        const int64_t id_r  = cell_id(i,jmax+j);
        const int64_t sym_r = cell_id(i,j+gx);

        rho  [id_l] = rho  [sym_l];
        rho_u[id_l] = rho_u[sym_l];
        rho_v[id_l] = rho_v[sym_l];
        E    [id_l] = E    [sym_l];

        rho  [id_r] = rho  [sym_r];
        rho_u[id_r] = rho_u[sym_r];
        rho_v[id_r] = rho_v[sym_r];
        E    [id_r] = E    [sym_r];
      }
  }
  else{
    const int64_t imax = sim->grid.imax;
    const int64_t gy   = sim->grid.gy;

    #pragma omp parallel for collapse(2)
    for(int64_t i=0; i < gy; i++)
      for(int64_t j=0; j < sim->grid.Nx_tot; j++){
        const int64_t id_l  = cell_id(i,j);
        const int64_t sym_l = cell_id(imax-gy+i,j);
        const int64_t id_r  = cell_id(imax+i,j);
        const int64_t sym_r = cell_id(gy+i,j);

        rho  [id_l] = rho  [sym_l];
        rho_u[id_l] = rho_u[sym_l];
        rho_v[id_l] = rho_v[sym_l];
        E    [id_l] = E    [sym_l];

        rho  [id_r] = rho  [sym_r];
        rho_u[id_r] = rho_u[sym_r];
        rho_v[id_r] = rho_v[sym_r];
        E    [id_r] = E    [sym_r];
      }
  }
}

// void boundary_periodic_old(struct sim *sim, enum dir dir)
// {
//   DECLARE_CSTATE_VAR
//   int64_t imaxl, jmaxl, lo;
//   int64_t iminr, imaxr, jminr, jmaxr, hi;
  
//   if(dir == IX){
//     imaxl = sim->grid.Ny_tot;
//     jmaxl = sim->grid.gx;
//     iminr = 0;
//     imaxr = sim->grid.Ny_tot;
//     jminr = sim->grid.Nx_tot - sim->grid.gx;
//     jmaxr = sim->grid.Nx_tot;
//     lo    = jmaxl;
//     hi    = jminr;
//   }
//   else{
//     imaxl = sim->grid.gy;
//     jmaxl = sim->grid.Nx_tot;
//     iminr = sim->grid.Ny_tot - sim->grid.gy;
//     imaxr = sim->grid.Ny_tot;
//     jminr = 0;
//     jmaxr = sim->grid.Nx_tot;
//     lo    = imaxl;
//     hi    = iminr;
//   }

//   // left side
//   // #pragma omp parallel for collapse(2)
//   for(int64_t i=0; i < imaxl; i++)
//     for(int64_t j=0; j < jmaxl; j++){
//       const int64_t id    = cell_id(i,j);
//       const int64_t id_l = (dir == IX) ? cell_id(i,hi-lo+j) : cell_id(hi-lo+i,j);
//       rho  [id] = rho  [id_l];
//       rho_u[id] = rho_u[id_l];
//       rho_v[id] = rho_v[id_l];
//       E    [id] = E    [id_l];
//     }
//   // right side
//   // #pragma omp parallel for collapse(2)
//   for(int64_t i=iminr; i < imaxr; i++)
//     for(int64_t j=jminr; j < jmaxr; j++){
//       const int64_t id    = cell_id(i,j);
//       const int64_t id_r = (dir == IX) ? cell_id(i,lo-hi+j) : cell_id(lo-hi+i,j);
//       rho  [id] = rho  [id_r];
//       rho_u[id] = rho_u[id_r];
//       rho_v[id] = rho_v[id_r];
//       E    [id] = E    [id_r];
//     }
// }

void boundary_absorbing(struct sim *sim, enum dir dir)
{
  DECLARE_CSTATE_VAR
  
  if(dir == IX){
    const int64_t jmax = sim->grid.jmax;
    const int64_t gx   = sim->grid.gx;
    
    #pragma omp parallel for collapse(2)
    for(int64_t i=0; i < sim->grid.Ny_tot; i++)
      for(int64_t j=0; j < gx; j++){
        const int64_t id_l  = cell_id(i,j);
        const int64_t sym_l = cell_id(i,gx);
        const int64_t id_r  = cell_id(i,jmax+j);
        const int64_t sym_r = cell_id(i,jmax-1);

        rho  [id_l] = rho  [sym_l];
        rho_u[id_l] = rho_u[sym_l];
        rho_v[id_l] = rho_v[sym_l];
        E    [id_l] = E    [sym_l];

        rho  [id_r] = rho  [sym_r];
        rho_u[id_r] = rho_u[sym_r];
        rho_v[id_r] = rho_v[sym_r];
        E    [id_r] = E    [sym_r];
      }
  }
  else{
    const int64_t imax = sim->grid.imax;
    const int64_t gy   = sim->grid.gy;

    #pragma omp parallel for collapse(2)
    for(int64_t i=0; i < gy; i++)
      for(int64_t j=0; j < sim->grid.Nx_tot; j++){
        const int64_t id_l  = cell_id(i,j);
        const int64_t sym_l = cell_id(gy,j);
        const int64_t id_r  = cell_id(imax+i,j);
        const int64_t sym_r = cell_id(imax-1,j);

        rho  [id_l] = rho  [sym_l];
        rho_u[id_l] = rho_u[sym_l];
        rho_v[id_l] = rho_v[sym_l];
        E    [id_l] = E    [sym_l];

        rho  [id_r] = rho  [sym_r];
        rho_u[id_r] = rho_u[sym_r];
        rho_v[id_r] = rho_v[sym_r];
        E    [id_r] = E    [sym_r];
      }
  }
}

// void boundary_absorbing_old(struct sim *sim, enum dir dir)
// {
//   DECLARE_CSTATE_VAR
//   int64_t imaxl, jmaxl, lo;
//   int64_t iminr, imaxr, jminr, jmaxr, hi;
  
//   if(dir == IX){
//     imaxl = sim->grid.Ny_tot;
//     jmaxl = sim->grid.gx;
//     iminr = 0;
//     imaxr = sim->grid.Ny_tot;
//     jminr = sim->grid.Nx_tot - sim->grid.gx;
//     jmaxr = sim->grid.Nx_tot;
//     lo    = jmaxl;
//     hi    = jminr-1;
//   }
//   else{
//     imaxl = sim->grid.gy;
//     jmaxl = sim->grid.Nx_tot;
//     iminr = sim->grid.Ny_tot - sim->grid.gy;
//     imaxr = sim->grid.Ny_tot;
//     jminr = 0;
//     jmaxr = sim->grid.Nx_tot;
//     lo    = imaxl;
//     hi    = iminr-1;
//   }

//   // left side
//   #pragma omp parallel for collapse(2)
//   for(int64_t i=0; i < imaxl; i++)
//     for(int64_t j=0; j < jmaxl; j++){
//       const int64_t id    = cell_id(i,j);
//       const int64_t id_l = (dir == IX) ? cell_id(i,lo) : cell_id(lo,j);
//       rho  [id] = rho  [id_l];
//       rho_u[id] = rho_u[id_l];
//       rho_v[id] = rho_v[id_l];
//       E    [id] = E    [id_l];
//     }
//   // right side
//   #pragma omp parallel for collapse(2)
//   for(int64_t i=iminr; i < imaxr; i++)
//     for(int64_t j=jminr; j < jmaxr; j++){
//       const int64_t id    = cell_id(i,j);
//       const int64_t id_r = (dir == IX) ? cell_id(i,hi) : cell_id(hi,j);
//       rho  [id] = rho  [id_r];
//       rho_u[id] = rho_u[id_r];
//       rho_v[id] = rho_v[id_r];
//       E    [id] = E    [id_r];
//     }
// }
