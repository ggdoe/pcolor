#pragma once
#include "type.h"
#include "utils.h"
#include "initial_conditions.h"
#include "boundary_conditions.h"

void prim_to_cons(struct sim *sim)
{
  DECLARE_STATES_VAR

  #pragma omp parallel for collapse(2)
  for_each_cells_and_ghost_y(i)
    for_each_cells_and_ghost_x(j)
    {
      const u64 id = cell_id(i,j);
      real_t Ec = 0.5 * prho[id] * (u[id] * u[id] + v[id] * v[id]);
      crho[id]  = prho[id];
      rho_u[id] = prho[id] * u[id];
      rho_v[id] = prho[id] * v[id];
      E[id]     = Ec + p[id] / (sim->gamma - 1.0);
    }
}

void cons_to_prim(struct sim *sim)
{
  DECLARE_STATES_VAR

  #pragma omp parallel for collapse(2)
  for_each_cells_and_ghost_y(i)
    for_each_cells_and_ghost_x(j)
    {
      const u64 id = cell_id(i,j);
      real_t Ec = 0.5 * (rho_u[id] * rho_u[id] + rho_v[id] * rho_v[id]) / crho[id];
      prho[id] = crho[id];
      u[id]    = rho_u[id] / crho[id];
      v[id]    = rho_v[id] / crho[id];
      p[id]    = (E[id] - Ec) * (sim->gamma - 1.0);
      // if(p[i] < 0)
      //   printf("negative pressure !!\n");
    }
}

struct grid init_grid(u32 Nx, u32 Ny, u32 gx, u32 gy)
{
  const real_t xmin = 0.0, ymin = 0.0;
  const real_t xmax = 1.0, ymax = 1.0;
  const real_t dx = (xmax - xmin) / Nx;
  const real_t dy = (ymax - ymin) / Ny;

  struct grid grid;
  grid.xmin = xmin;         grid.ymin = ymin;
  grid.xmax = xmax;         grid.ymax = ymax;

  grid.Nx_tot = Nx + 2*gx;  grid.Ny_tot = Ny + 2*gy;
  grid.Nx     = Nx;         grid.Ny     = Ny;
  grid.dx     = dx;         grid.dy     = dy;
  grid.gx     = gx;         grid.gy     = gy;
  grid.jmin   = gx;         grid.imin   = gy;
  grid.jmax   = Nx + gx;    grid.imax   = Ny + gy;

  grid.cellcenter_x = (real_t*)MALLOC(grid.Nx_tot * grid.Ny_tot * sizeof(real_t));
  grid.cellcenter_y = (real_t*)MALLOC(grid.Nx_tot * grid.Ny_tot * sizeof(real_t));
  for(u64 i=0; i < grid.Ny_tot; i++){
    for(u64 j=0; j < grid.Nx_tot; j++){
      grid.cellcenter_x[i * grid.Nx_tot + j] = xmin + (0.5 - gx + j) * dx;
      grid.cellcenter_y[i * grid.Nx_tot + j] = ymin + (0.5 - gy + i) * dy;
    }
  }

  grid.vertex_x = (real_t*)MALLOC((grid.Nx_tot+1) * (grid.Ny_tot+1) * sizeof(real_t));
  grid.vertex_y = (real_t*)MALLOC((grid.Nx_tot+1) * (grid.Ny_tot+1) * sizeof(real_t));
  for(u64 i=0; i < grid.Ny_tot+1; i++){
    for(u64 j=0; j < grid.Nx_tot+1; j++){
      grid.vertex_x[i * (grid.Nx_tot+1) + j] = xmin + (j - (real_t)gx) * dx;
      grid.vertex_y[i * (grid.Nx_tot+1) + j] = ymin + (i - (real_t)gy) * dy;
    }
  }

  return grid;
}

void alloc_state(void *state, size_t n)
{
  struct pstate *s = (struct pstate*)state;
  s->rho = (real_t*)MALLOC(n * sizeof(real_t));
  s->u   = (real_t*)MALLOC(n * sizeof(real_t));
  s->v   = (real_t*)MALLOC(n * sizeof(real_t));
  s->p   = (real_t*)MALLOC(n * sizeof(real_t));
  memset(s->rho, 0, n * sizeof(real_t));
  memset(s->u  , 0, n * sizeof(real_t));
  memset(s->v  , 0, n * sizeof(real_t));
  memset(s->p  , 0, n * sizeof(real_t));
}

struct sim init_sim(u32 Nx, u32 Ny)
{
  struct sim sim;
  const u32 gx = 2, gy = 2;
  sim.grid = init_grid(Nx, Ny, gx, gy);
  sim.gamma = 1.4;
  sim.cfl = 0.1;
  sim.t = 0.0;

  alloc_state(&sim.cstate,  sim.grid.Nx_tot * sim.grid.Ny_tot);
  alloc_state(&sim.pstate,  sim.grid.Nx_tot * sim.grid.Ny_tot);
  alloc_state(&sim.slope_x, sim.grid.Nx_tot * sim.grid.Ny_tot);
  alloc_state(&sim.slope_y, sim.grid.Nx_tot * sim.grid.Ny_tot);
  alloc_state(&sim.flux_x,  sim.grid.Nx_tot * sim.grid.Ny_tot);
  alloc_state(&sim.flux_y,  sim.grid.Nx_tot * sim.grid.Ny_tot);

  return sim;
}

void free_state(void *state)
{
  struct pstate *s = (struct pstate*)state;
  free(s->rho);
  free(s->u);
  free(s->v);
  free(s->p);
}

void free_sim(struct sim *sim)
{
  free(sim->grid.cellcenter_x);
  free(sim->grid.cellcenter_y);
  free(sim->grid.vertex_x);
  free(sim->grid.vertex_y);
  free_state(&sim->cstate);
  free_state(&sim->pstate);
  free_state(&sim->slope_x);
  free_state(&sim->slope_y);
  free_state(&sim->flux_x);
  free_state(&sim->flux_y);
}

real_t compute_dt(struct sim *sim)
{
  DECLARE_PSTATE_VAR
  real_t dt = __FLT_MAX__;
  const double dx = sim->grid.dx;
  const double dy = sim->grid.dy;
  const double gamma = sim->gamma;

  #pragma omp parallel for reduction(min:dt) collapse(2)
  for_each_cells_y(i)
    for_each_cells_x(j)
    {
      const u64 id = cell_id(i,j);
      real_t cs = sqrt(gamma * p[id] / rho[id]);
      real_t dtx = dx / (cs + fabs(u[id]));
      real_t dty = dy / (cs + fabs(v[id]));
      real_t local_dt = fmin(dtx, dty);
      dt = fmin(dt, local_dt);
    }
  return sim->cfl * dt;
}

struct fcell riemann_hllc(struct pcell *restrict left, struct pcell *restrict right, real_t gamma)
{
  const real_t entho = 1.0 / (gamma - 1.0);

  real_t rl = left->rho;
  real_t ul = left->u;
  real_t vl = left->v;
  real_t pl = left->p;
  real_t ecinl = 0.5*rl*(ul*ul+vl*vl);
  real_t etotl = ecinl + pl*entho;

  real_t rr = right->rho;
  real_t ur = right->u;
  real_t vr = right->v;
  real_t pr = right->p;
  real_t ecinr = 0.5*rr*(ur*ur+vr*vr);
  real_t etotr = ecinr + pr*entho;

  real_t cl = sqrt(gamma*pl/rl);
  real_t cr = sqrt(gamma*pr/rr);

  real_t SL = fmin(ul, ur) - fmax(cl, cr);
  real_t SR = fmin(ul, ur) + fmax(cl, cr);

  // Compute lagrangian sound speed
  real_t rcl = rl*(ul-SL);
  real_t rcr = rr*(SR-ur);
    
  // Compute acoustic star state
  real_t ustar    = (rcr*ur   +rcl*ul   +  (pl-pr))/(rcr+rcl);
  real_t ptotstar = (rcr*pl+rcl*pr+rcl*rcr*(ul-ur))/(rcr+rcl);

  // Left star region variables
  real_t rstarl    = rl*(SL-ul)/(SL-ustar);
  real_t etotstarl = ((SL-ul)*etotl-pl*ul+ptotstar*ustar)/(SL-ustar);
    
  // Right star region variables
  real_t rstarr    = rr*(SR-ur)/(SR-ustar);
  real_t etotstarr = ((SR-ur)*etotr-pr*ur+ptotstar*ustar)/(SR-ustar);
    
  // Sample the solution at x/t=0
  real_t ro, uo, vo, ptoto, etoto;
  if (SL > 0) {
    ro=rl;
    uo=ul;
    vo=vl;
    ptoto=pl;
    etoto=etotl;
  } else if (ustar > 0) {
    ro=rstarl;
    uo=ustar;
    vo=vl;
    ptoto=ptotstar;
    etoto=etotstarl;
  } else if (SR > 0) {
    ro=rstarr;
    uo=ustar;
    vo=vr;
    ptoto=ptotstar;
    etoto=etotstarr;
  } else {
    ro=rr;
    uo=ur;
    vo=vr;
    ptoto=pr;
    etoto=etotr;
  }
  struct fcell flux;

  // Compute the Godunov flux
  flux.rho   = ro*uo;
  flux.rho_u = ro*uo*uo+ptoto;
  flux.rho_v = flux.rho*vo;
  flux.E = (etoto+ptoto)*uo;

  return flux;
}

real_t minmod(real_t a, real_t b)
{
    return (a * b < 0.0) ? 0.0 : (fabs(a) < fabs(b)) ? a : b;
}

real_t limited_slope(real_t *q, u64 id, const u64 off)
{
  return minmod(q[id] - q[id-off], q[id+off] - q[id]);
}

struct interface_values{
  struct pcell left;
  struct pcell right;
};

struct interface_values get_interface_values(struct sim *sim, u64 id, enum dir dir)
{
  DECLARE_PSTATE_VAR
  struct pstate *slope = (dir == IX) ? &sim->slope_x : &sim->slope_y;
  const u64        off = (dir == IX) ? 1 : sim->grid.Nx_tot;

  struct interface_values iv = {
    .left = {
      .rho = rho[id]     + 0.5 * slope->rho[id],
      .u   = u  [id]     + 0.5 * slope->u  [id],
      .v   = v  [id]     + 0.5 * slope->v  [id],
      .p   = p  [id]     + 0.5 * slope->p  [id]
    },
    .right = {
      .rho = rho[id+off] - 0.5 * slope->rho[id+off],
      .u   = u  [id+off] - 0.5 * slope->u  [id+off],
      .v   = v  [id+off] - 0.5 * slope->v  [id+off],
      .p   = p  [id+off] - 0.5 * slope->p  [id+off]
    }
  };
  
  return iv;
}

void compute_slope(struct sim *sim)
{
  DECLARE_PSTATE_VAR
  const u64 offset_x = 1;
  const u64 offset_y = sim->grid.Nx_tot;

  #pragma omp parallel for collapse(2)
  for(u64 i=sim->grid.jmin-1; i < sim->grid.imax+1  ; i++)   // +1 for muscl hancock,
    for(u64 j=sim->grid.imin-1; j < sim->grid.jmax+1  ; j++){// else just use for_each_interface
      const u64 id = cell_id(i,j);
      // x-dir
      sim->slope_x.rho[id] = limited_slope(rho, id, offset_x);
      sim->slope_x.  u[id] = limited_slope(  u, id, offset_x);
      sim->slope_x.  v[id] = limited_slope(  v, id, offset_x);
      sim->slope_x.  p[id] = limited_slope(  p, id, offset_x);

      // y-dir
      sim->slope_y.rho[id] = limited_slope(rho, id, offset_y);
      sim->slope_y.  u[id] = limited_slope(  u, id, offset_y);
      sim->slope_y.  v[id] = limited_slope(  v, id, offset_y);
      sim->slope_y.  p[id] = limited_slope(  p, id, offset_y);
    }
}

void swap_uv(struct pcell* cell)
{
  real_t swap = cell->u;
  cell->u = cell->v;
  cell->v = swap;
}

void compute_fluxes(struct sim *sim, enum dir dir)
{
  struct cstate *flux = (dir == IX) ? &sim->flux_x : &sim->flux_y;

  #pragma omp parallel for collapse(2)
  for_each_interfaces_y(i)
    for_each_interfaces_x(j)
    {
      const u64 id = cell_id(i,j);
      struct interface_values iv = get_interface_values(sim, id, dir);
      if(dir == IY){
        swap_uv(&iv.left);
        swap_uv(&iv.right);
      }

      struct fcell riemannflux = riemann_hllc(&iv.left, &iv.right, sim->gamma);

      if(dir == IY)
        swap_uv((struct pcell*)&riemannflux);

      flux->rho  [id] = riemannflux.rho;
      flux->rho_u[id] = riemannflux.rho_u;
      flux->rho_v[id] = riemannflux.rho_v;
      flux->E    [id] = riemannflux.E;
    }
}

void reconstruct_muscl_hancock(struct sim *sim, real_t dt)
{
  DECLARE_PSTATE_VAR
  const real_t dtdx = 0.5 * dt / sim->grid.dx;
  const real_t dtdy = 0.5 * dt / sim->grid.dy;

  #pragma omp parallel for collapse(2)
  for(u64 i=sim->grid.jmin-1; i < sim->grid.imax+1  ; i++)
    for(u64 j=sim->grid.imin-1; j < sim->grid.jmax+1  ; j++)
    {
      const u64 id = cell_id(i,j);
      const real_t dxr = sim->slope_x.rho[id];
      const real_t dxu = sim->slope_x.  u[id];
      const real_t dxv = sim->slope_x.  v[id];
      const real_t dxp = sim->slope_x.  p[id];

      const real_t dyr = sim->slope_y.rho[id];
      const real_t dyu = sim->slope_y.  u[id];
      const real_t dyv = sim->slope_y.  v[id];
      const real_t dyp = sim->slope_y.  p[id];

      const real_t pr  = rho[id];
      const real_t pu  =   u[id];
      const real_t pv  =   v[id];
      const real_t pp  =   p[id];
    
      rho[id] = pr + (- pu * dxr - pr * dxu) * dtdx
                   + (- pv * dyr - pr * dyv) * dtdy;
      u[id]   = pu + (- pu * dxu - dxp / pr) * dtdx
                   + (- pv * dyu) * dtdy;
      v[id]   = pv + (- pu * dxv) * dtdx
                   + (- pv * dyv - dyp / pr) * dtdy;
      p[id]   = pp + (- sim->gamma * pp * dxu - pu * dxp) * dtdx
                   + (- sim->gamma * pp * dyv - pv * dyp) * dtdy;
    }
}

void cells_update(struct sim *sim, real_t dt)
{
  DECLARE_CSTATE_VAR
  const real_t dtdx = dt / sim->grid.dx;
  const real_t dtdy = dt / sim->grid.dy;
  const u64 off_x = 1;
  const u64 off_y = sim->grid.Nx_tot;
  struct cstate *fx = &sim->flux_x;
  struct cstate *fy = &sim->flux_y;

  #pragma omp parallel for collapse(2)
  for_each_cells_y(i)
    for_each_cells_x(j){
      const u64 id = cell_id(i,j);
      rho  [id] += dtdx * (fx->rho  [id-off_x] - fx->rho  [id])
                 + dtdy * (fy->rho  [id-off_y] - fy->rho  [id]);
      rho_u[id] += dtdx * (fx->rho_u[id-off_x] - fx->rho_u[id])
                 + dtdy * (fy->rho_u[id-off_y] - fy->rho_u[id]);
      rho_v[id] += dtdx * (fx->rho_v[id-off_x] - fx->rho_v[id])
                 + dtdy * (fy->rho_v[id-off_y] - fy->rho_v[id]);
      E    [id] += dtdx * (fx->E    [id-off_x] - fx->E    [id])
                 + dtdy * (fy->E    [id-off_y] - fy->E    [id]);
  }
}

void step(struct sim *sim, real_t dt_max)
{
  real_t dt = compute_dt(sim);
  dt = (dt < dt_max) ? dt : dt_max;
  sim->t += dt;
  
  compute_slope(sim);
  reconstruct_muscl_hancock(sim, dt);
  compute_fluxes(sim, IX);
  compute_fluxes(sim, IY);
  cells_update(sim, dt);
  fill_boundaries(sim);
  cons_to_prim(sim);
  // printf("t: %.4lf\tdt: %.3e\n", sim->t, dt);
}

void run(struct sim *sim, real_t tmax)
{
  while(sim->t < tmax)
  {
    real_t dt_max = (tmax - sim->t);
    step(sim, dt_max);
  }
}
  
