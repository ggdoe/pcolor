#include "type.h"

#define for_each_cells(i)       for(u64 i=sim->grid.imin; i < sim->grid.imax; i++)
#define for_each_cells_and_ghost(i)    for(u64 i=0; i < sim->grid.Nx_tot; i++)
#define for_each_interfaces(i)  for(u64 i=sim->grid.imin-1; i < sim->grid.imax+1; i++)

struct cstate{
  real_t *restrict rho;
  real_t *restrict rho_u;
  real_t *restrict E;
};

struct pstate{
  real_t *restrict rho;
  real_t *restrict u;
  real_t *restrict p;
};

struct grid {
  real_t xmin;
  real_t xmax;

  u32 Nx;
  u32 Nx_tot;
  u32 gx;
  real_t dx;
  real_t *cellcenter;

  u32 imin;
  u32 imax;
};

struct sim{
  struct grid grid;
  struct cstate cstate;
  struct pstate pstate;
  struct pstate slope;
  struct cstate fluxe;
  real_t gamma;
  real_t cfl;
  real_t t;
};

struct pcell {
  real_t rho;
  real_t u;
  real_t p;
};

struct fcell{
  real_t rho;
  real_t rho_u;
  real_t E;
};

void prim_to_cons(struct sim *sim)
{
  #pragma omp parallel for
  for_each_cells_and_ghost(i)
  {
    real_t Ec = 0.5 * sim->pstate.rho[i] * (sim->pstate.u[i] * sim->pstate.u[i]);
    sim->cstate.rho[i]   = sim->pstate.rho[i];
    sim->cstate.rho_u[i] = sim->pstate.rho[i] * sim->pstate.u[i];
    sim->cstate.E[i]     = Ec + sim->pstate.p[i] / (sim->gamma - 1.0);
  }
}

void cons_to_prim(struct sim *sim)
{
  #pragma omp parallel for
  for_each_cells_and_ghost(i)
  {
    real_t Ec = 0.5 * (sim->cstate.rho_u[i] * sim->cstate.rho_u[i]) / sim->cstate.rho[i];
    sim->pstate.rho[i] = sim->cstate.rho[i];
    sim->pstate.u[i]   = sim->cstate.rho_u[i] / sim->cstate.rho[i];
    sim->pstate.p[i]   = (sim->cstate.E[i] - Ec) * (sim->gamma - 1.0);
    // if(sim->pstate.p[i] < 0)
    //   printf("negative pressure !!\n");
  }
}

void init_state(struct sim *sim)
{
  sim->t = 0;

  for_each_cells_and_ghost(i)
  {
    if(sim->grid.cellcenter[i] < 0.5 * (sim->grid.xmax + sim->grid.xmin)){
      sim->pstate.rho[i] = 1.0;
      sim->pstate.u[i]   = 0;
      sim->pstate.p[i]   = 1.0;
    }
    else{
      sim->pstate.rho[i] = 0.125;
      sim->pstate.u[i]   = 0;
      sim->pstate.p[i]   = 0.1;
    }
  }
  prim_to_cons(sim);
}

struct grid init_grid(u32 Nx, u32 gx)
{
  real_t xmin = 0.0;
  real_t xmax = 1.0;
  real_t dx = (xmax - xmin) / Nx;

  struct grid grid;
  grid.xmin = xmin;
  grid.xmax = xmax;

  grid.Nx_tot = Nx + 2*gx;
  grid.Nx = Nx;
  grid.dx = dx;
  grid.gx = gx;
  grid.imin = gx;
  grid.imax = Nx + gx;

  grid.cellcenter = (real_t*)malloc(grid.Nx_tot * sizeof(real_t));
  for(u64 i=0; i < grid.Nx_tot; i++) // todo : boudary condition for ghost cell
    grid.cellcenter[i] = xmin + (0.5 - gx + i) * dx;

  return grid;
}

void alloc_state(void *state, size_t n)
{
  struct pstate *s = (struct pstate*)state;
  s->rho = (real_t*)malloc(n * sizeof(real_t));
  s->u   = (real_t*)malloc(n * sizeof(real_t));
  s->p   = (real_t*)malloc(n * sizeof(real_t));
  memset(s->rho, 0, n * sizeof(real_t));
  memset(s->u  , 0, n * sizeof(real_t));
  memset(s->p  , 0, n * sizeof(real_t));
}

struct sim init_sim(u32 Nx)
{
  struct sim sim;
  const u32 gx = 2;
  sim.grid = init_grid(Nx, gx);
  sim.gamma = 1.4;
  sim.cfl = 0.02;
  sim.t = 0.0;
  
  alloc_state(&sim.cstate, sim.grid.Nx_tot);
  alloc_state(&sim.pstate, sim.grid.Nx_tot);
  alloc_state(&sim.slope,  sim.grid.Nx_tot);
  alloc_state(&sim.fluxe,  sim.grid.Nx_tot);

  return sim;
}

void free_state(void *state)
{
  struct pstate *s = (struct pstate*)state;
  free(s->rho);
  free(s->u  );
  free(s->p  );
}

void free_sim(struct sim *sim)
{
  free(sim->grid.cellcenter);
  free_state(&sim->cstate);
  free_state(&sim->pstate);
  free_state(&sim->slope);
  free_state(&sim->fluxe);
}

real_t compute_dt(struct sim *sim)
{
  real_t max_speed = 0.0;
  #pragma omp parallel for reduction(max:max_speed)
  for_each_cells(i)
  {
    real_t cs = sqrt(sim->pstate.p[i] * sim->gamma / sim->pstate.rho[i]);
    real_t speed = cs + fabs(sim->pstate.u[i]);
    max_speed = (speed > max_speed) ? speed : max_speed;
  }
  const real_t dt = sim->grid.dx / max_speed;
  return sim->cfl * dt;
}

struct fcell riemann_hllc(struct pcell *restrict left, struct pcell *restrict right, real_t gamma)
{
  const real_t entho = 1.0 / (gamma - 1.0);

  real_t rl = left->rho;
  real_t ul = left->u;
  real_t pl = left->p;
  real_t ecinl = 0.5*rl*(ul*ul);
  real_t etotl = ecinl + pl*entho;

  real_t rr = right->rho;
  real_t ur = right->u;
  real_t pr = right->p;
  real_t ecinr = 0.5*rr*(ur*ur);
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
  flux.E = (etoto+ptoto)*uo;

  return flux;
}

real_t minmod(real_t a, real_t b)
{
    return (a * b < 0.0) ? 0.0 : (fabs(a) < fabs(b)) ? a : b;
}

real_t limited_slope(real_t *q, u64 id)
{
  return minmod(q[id] - q[id-1], q[id+1] - q[id]);
}

struct interface_values{
  struct pcell left;
  struct pcell right;
};

struct interface_values get_interface_values(struct sim *sim, u64 id)
{
  struct interface_values iv = {
    .left = {
      .rho = sim->pstate.rho[id]   + 0.5 * sim->slope.rho[id],
      .u   = sim->pstate.u  [id]   + 0.5 * sim->slope.u  [id],
      .p   = sim->pstate.p  [id]   + 0.5 * sim->slope.p  [id]
    },
    .right = {
      .rho = sim->pstate.rho[id+1] - 0.5 * sim->slope.rho[id+1],
      .u   = sim->pstate.u  [id+1] - 0.5 * sim->slope.u  [id+1],
      .p   = sim->pstate.p  [id+1] - 0.5 * sim->slope.p  [id+1]
    }
  };
  
  return iv;
}

void compute_slope(struct sim *sim)
{
  #pragma omp parallel for
  for_each_cells(i){
    sim->slope.rho[i] = limited_slope(sim->pstate.rho, i);
    sim->slope.u[i]   = limited_slope(sim->pstate.u, i);
    sim->slope.p[i]   = limited_slope(sim->pstate.p, i);
  }
}

void compute_fluxes(struct sim *sim)
{
  #pragma omp parallel for
  for_each_interfaces(i)
  {
    struct interface_values iv = get_interface_values(sim, i);
    struct fcell flux = riemann_hllc(&iv.left, &iv.right, sim->gamma);

    sim->fluxe.rho  [i] = flux.rho;
    sim->fluxe.rho_u[i] = flux.rho_u;
    sim->fluxe.E    [i] = flux.E;
  }
}

void cells_update(struct sim *sim, real_t dt)
{
  const real_t dtdx = dt / sim->grid.dx;

  #pragma omp parallel for
  for_each_cells(i)
  {
    sim->cstate.rho  [i] += dtdx * (sim->fluxe.rho  [i-1] - sim->fluxe.rho  [i]);
    sim->cstate.rho_u[i] += dtdx * (sim->fluxe.rho_u[i-1] - sim->fluxe.rho_u[i]);
    sim->cstate.E    [i] += dtdx * (sim->fluxe.E    [i-1] - sim->fluxe.E    [i]);
  }
}

void fill_boundaries_absorbing(struct sim *sim)
{
  const u64 lo = sim->grid.gx;
  const u64 hi = sim->grid.Nx_tot - sim->grid.gx;

  #pragma omp parallel for
  for(u64 i=0; i < lo; i++){
    sim->cstate.rho  [i] = sim->cstate.rho  [lo];
    sim->cstate.rho_u[i] = sim->cstate.rho_u[lo];
    sim->cstate.E    [i] = sim->cstate.E    [lo];
  }
  #pragma omp parallel for
  for(u64 i=hi; i < sim->grid.Nx_tot; i++){
    sim->cstate.rho  [i] = sim->cstate.rho  [hi-1];
    sim->cstate.rho_u[i] = sim->cstate.rho_u[hi-1];
    sim->cstate.E    [i] = sim->cstate.E    [hi-1];
  }
}

void step(struct sim *sim, real_t dt_max)
{
  real_t dt = compute_dt(sim);
  dt = (dt < dt_max) ? dt : dt_max;
  sim->t += dt;
  
  compute_slope(sim);
  compute_fluxes(sim);
  cells_update(sim, dt);
  fill_boundaries_absorbing(sim);
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
  
