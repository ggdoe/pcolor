#include "type.h"

#define for_each_cells_x(j)       for(u64 j=sim->grid.jmin; j < sim->grid.jmax; j++)
#define for_each_cells_y(i)       for(u64 i=sim->grid.imin; i < sim->grid.imax; i++)
#define for_each_cells_and_ghost_x(j)    for(u64 j=0; j < sim->grid.Nx_tot; j++)
#define for_each_cells_and_ghost_y(i)    for(u64 i=0; i < sim->grid.Ny_tot; i++)
#define for_each_interfaces_x(j)  for(u64 j=sim->grid.imin-1; j < sim->grid.jmax+1; j++)
#define for_each_interfaces_y(i)  for(u64 i=sim->grid.jmin-1; i < sim->grid.imax+1; i++)
#define cell_id(i,j) (i*sim->grid.Nx_tot + j)

#define DECLARE_PSTATE_VAR              \
        double *rho = sim->pstate.rho;  \
        double *u   = sim->pstate.u;    \
        double *v   = sim->pstate.v;    \
        double *p   = sim->pstate.p;

#define DECLARE_CSTATE_VAR                  \
        double *rho   = sim->cstate.rho;    \
        double *rho_u = sim->cstate.rho_u;  \
        double *rho_v = sim->cstate.rho_v;  \
        double *E     = sim->cstate.E;
        
#define DECLARE_STATES_VAR                  \
        double *rho   = sim->pstate.rho;    \
        double *u     = sim->pstate.u;      \
        double *v     = sim->pstate.v;      \
        double *p     = sim->pstate.p;      \
        double *rho_u = sim->cstate.rho_u;  \
        double *rho_v = sim->cstate.rho_v;  \
        double *E     = sim->cstate.E;

enum dir{
  IX, IY
};

struct cstate{
  real_t *restrict rho;
  real_t *restrict rho_u;
  real_t *restrict rho_v;
  real_t *restrict E;
};

struct pstate{
  real_t *restrict rho;
  real_t *restrict u;
  real_t *restrict v;
  real_t *restrict p;
};

struct grid {
  real_t xmin, ymin;
  real_t xmax, ymax;

  u32 Nx, Ny;
  u32 Nx_tot, Ny_tot;
  u32 gx, gy;
  real_t dx, dy;

  real_t *vertex_x, *vertex_y;
  real_t *cellcenter_x, *cellcenter_y;

  u32 jmin, imin;
  u32 jmax, imax;
};

struct sim{
  struct grid grid;
  struct cstate cstate;
  struct pstate pstate;
  struct pstate slope_x;
  struct pstate slope_y;
  struct cstate flux_x;
  struct cstate flux_y;
  real_t gamma;
  real_t cfl;
  real_t t;
};

struct pcell {
  real_t rho;
  real_t u;
  real_t v;
  real_t p;
};

struct fcell{
  real_t rho;
  real_t rho_u;
  real_t rho_v;
  real_t E;
};

void prim_to_cons(struct sim *sim)
{
  DECLARE_STATES_VAR

  //#pragma omp parallel for
  for_each_cells_and_ghost_y(i)
    for_each_cells_and_ghost_x(j)
    {
      const u64 id = cell_id(i,j);
      real_t Ec = 0.5 * rho[id] * (u[id] * u[id] + v[id] * v[id]);
      rho[id]   = rho[id];
      rho_u[id] = rho[id] * u[id];
      rho_v[id] = rho[id] * v[id];
      E[id]     = Ec + p[id] / (sim->gamma - 1.0);
    }
}

void cons_to_prim(struct sim *sim)
{
  DECLARE_STATES_VAR

  //#pragma omp parallel for
  for_each_cells_and_ghost_y(i)
    for_each_cells_and_ghost_x(j)
    {
      const u64 id = cell_id(i,j);
      real_t Ec = 0.5 * (rho_u[id] * rho_u[id] + rho_v[id] * rho_v[id]) / rho[id];
      rho[id] = rho[id];
      u[id]   = rho_u[id] / rho[id];
      v[id]   = rho_v[id] / rho[id];
      p[id]   = (E[id] - Ec) * (sim->gamma - 1.0);
      // if(p[i] < 0)
      //   printf("negative pressure !!\n");
    }
}

void init_state(struct sim *sim)
{
  DECLARE_PSTATE_VAR
  sim->t = 0;

  for_each_cells_and_ghost_y(i)
    for_each_cells_and_ghost_x(j)
    {
      const u64 id = cell_id(i,j);
      if(sim->grid.cellcenter_x[id] < 0.5 * (sim->grid.xmax + sim->grid.xmin)){
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

struct grid init_grid(u32 Nx, u32 Ny, u32 gx, u32 gy)
{
  const real_t xmin = 0.0, ymin = 0.0;
  const real_t xmax = 1.0, ymax = 1.0;
  const real_t dx = (xmax - xmin) / Nx;
  const real_t dy = (ymax - ymin) / Ny;

  struct grid grid;
  grid.xmin = xmin; grid.ymin = ymin;
  grid.xmax = xmax; grid.ymax = ymax;

  grid.Nx_tot = Nx + 2*gx;  grid.Ny_tot = Ny + 2*gy;
  grid.Nx = Nx;             grid.Ny = Ny;
  grid.dx = dx;             grid.dy = dy;
  grid.gx = gx;             grid.gy = gy;
  grid.jmin = gx;           grid.imin = gy;
  grid.jmax = Nx + gx;      grid.imax = Ny + gy;

  grid.cellcenter_x = (real_t*)malloc(grid.Nx_tot * grid.Ny_tot * sizeof(real_t));
  grid.cellcenter_y = (real_t*)malloc(grid.Nx_tot * grid.Ny_tot * sizeof(real_t));
  for(u64 i=0; i < grid.Ny_tot; i++){
    for(u64 j=0; j < grid.Nx_tot; j++){
      grid.cellcenter_x[i * grid.Nx_tot + j] = xmin + (0.5 - gx + j) * dx;
      grid.cellcenter_y[i * grid.Nx_tot + j] = ymin + (0.5 - gy + i) * dy;
    }
  }

  grid.vertex_x = (real_t*)malloc((grid.Nx_tot+1) * (grid.Ny_tot+1) * sizeof(real_t));
  grid.vertex_y = (real_t*)malloc((grid.Nx_tot+1) * (grid.Ny_tot+1) * sizeof(real_t));
  for(u64 i=0; i < grid.Ny_tot+1; i++){
    for(u64 j=0; j < grid.Nx_tot+1; j++){
      grid.vertex_x[i * (grid.Nx_tot+1) + j] = xmin + (j - gx) * dx;
      grid.vertex_y[i * (grid.Nx_tot+1) + j] = ymin + (i - gy) * dy;
    }
  }

  return grid;
}

void alloc_state(void *state, size_t n)
{
  struct pstate *s = (struct pstate*)state;
  s->rho = (real_t*)malloc(n * sizeof(real_t));
  s->u   = (real_t*)malloc(n * sizeof(real_t));
  s->v   = (real_t*)malloc(n * sizeof(real_t));
  s->p   = (real_t*)malloc(n * sizeof(real_t));
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
  sim.cfl = 0.2;
  sim.t = 0.0;

  alloc_state(&sim.cstate, sim.grid.Nx_tot * sim.grid.Ny_tot);
  alloc_state(&sim.pstate, sim.grid.Nx_tot * sim.grid.Ny_tot);
  alloc_state(&sim.slope_x,  sim.grid.Nx_tot * sim.grid.Ny_tot);
  alloc_state(&sim.slope_y,  sim.grid.Nx_tot * sim.grid.Ny_tot);
  alloc_state(&sim.flux_x, sim.grid.Nx_tot * sim.grid.Ny_tot);
  alloc_state(&sim.flux_y, sim.grid.Nx_tot * sim.grid.Ny_tot);

  return sim;
}

void free_state(void *state)
{
  struct pstate *s = (struct pstate*)state;
  free(s->rho);
  free(s->u  );
  free(s->v  );
  free(s->p  );
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

  //#pragma omp parallel for reduction(min:inv_dt)
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
  real_t ro, uo, ptoto, etoto;
  if (SL > 0) {
    ro=rl;
    uo=ul;
    ptoto=pl;
    etoto=etotl;
  } else if (ustar > 0) {
    ro=rstarl;
    uo=ustar;
    ptoto=ptotstar;
    etoto=etotstarl;
  } else if (SR > 0) {
    ro=rstarr;
    uo=ustar;
    ptoto=ptotstar;
    etoto=etotstarr;
  } else {
    ro=rr;
    uo=ur;
    ptoto=pr;
    etoto=etotr;
  }
  struct fcell flux;

  // Compute the Godunov flux
  flux.rho   = ro*uo;
  flux.rho_u = ro*uo*uo+ptoto;
  flux.rho_v = (flux.rho > 0) ? flux.rho*vl : flux.rho*vr;
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

  //#pragma omp parallel for
  for_each_cells_y(i)
    for_each_cells_x(j){
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

void swap_uv(void* state)
{
  struct pcell *s = state;
  real_t swap = s->u;
  s->u = s->v;
  s->v = swap;
}

void compute_fluxes(struct sim *sim, enum dir dir)
{
  struct cstate *flux = (dir == IX) ? &sim->flux_x : &sim->flux_y;

  //#pragma omp parallel for
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
        swap_uv(&riemannflux);
        
      flux->rho  [id] = riemannflux.rho;
      flux->rho_u[id] = riemannflux.rho_u;
      flux->rho_v[id] = riemannflux.rho_v;
      flux->E    [id] = riemannflux.E;
    }
}

void cells_update(struct sim *sim, real_t dt)
{
  DECLARE_CSTATE_VAR
  const real_t dtdx = dt / sim->grid.dx;
  const real_t dtdy = dt / sim->grid.dy;
  const u64 off_x = 1;
  const u64 off_y = sim->grid.Nx_tot;
  struct cstate *flux_x = &sim->flux_x;
  struct cstate *flux_y = &sim->flux_y;

  //#pragma omp parallel for
  for_each_cells_y(i)
    for_each_cells_x(j){
      const u64 id = cell_id(i,j);
      rho  [id] += dtdx * (flux_x->rho  [id-off_x] - flux_x->rho  [id])
                 + dtdy * (flux_y->rho  [id-off_y] - flux_y->rho  [id]);
      rho_u[id] += dtdx * (flux_x->rho_u[id-off_x] - flux_x->rho_u[id])
                 + dtdy * (flux_y->rho_u[id-off_y] - flux_y->rho_u[id]);
      rho_v[id] += dtdx * (flux_x->rho_v[id-off_x] - flux_x->rho_v[id])
                 + dtdy * (flux_y->rho_v[id-off_y] - flux_y->rho_v[id]);
      E    [id] += dtdx * (flux_x->E    [id-off_x] - flux_x->E    [id])
                 + dtdy * (flux_y->E    [id-off_y] - flux_y->E    [id]);
  }
}

void fill_boundaries_absorbing(struct sim *sim, enum dir dir)
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
  //#pragma omp parallel for
  for(u64 i=0; i < imaxl; i++)
    for(u64 j=0; j < jmaxl; j++){
      const u64 id    = cell_id(i,j);
      const u64 id_lo = cell_id(i,lo);
      rho  [id] = rho  [id_lo];
      rho_u[id] = rho_u[id_lo];
      rho_v[id] = rho_v[id_lo];
      E    [id] = E    [id_lo];
    }
  // right side
  //#pragma omp parallel for
  for(u64 i=iminr; i < imaxr; i++)
    for(u64 j=jminr; j < jmaxr; j++){
      const u64 id    = cell_id(i,j);
      const u64 id_hi = cell_id(i,hi-1);
      rho  [id] = rho  [id_hi];
      rho_u[id] = rho_u[id_hi];
      rho_v[id] = rho_v[id_hi];
      E    [id] = E    [id_hi];
    }
}

void step(struct sim *sim, real_t dt_max)
{
  real_t dt = compute_dt(sim);
  dt = (dt < dt_max) ? dt : dt_max;
  sim->t += dt;
  
  compute_slope(sim);
  compute_fluxes(sim, IX);
  compute_fluxes(sim, IY);
  cells_update(sim, dt);
  fill_boundaries_absorbing(sim, IX);
  fill_boundaries_absorbing(sim, IY);
  cons_to_prim(sim);
}

void run(struct sim *sim, real_t tmax)
{
  while(sim->t < tmax)
  {
    real_t dt_max = (tmax - sim->t);
    step(sim, dt_max);
  }
}
  
