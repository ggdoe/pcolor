#pragma once
enum dir{
  IX, IY
};

enum boundary_type{
  ABSORBING, PERIODIC, REFLECTING,
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

  int64_t    Nx,     Ny;
  int64_t    Nx_tot, Ny_tot;
  int64_t    gx,     gy;
  real_t dx,     dy;

  real_t *vertex_x,     *vertex_y;
  real_t *cellcenter_x, *cellcenter_y;

  int64_t jmin, imin;
  int64_t jmax, imax;
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
  enum boundary_type boundary[2];
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

#define for_each_cells_x(j)            for(u64 j=sim->grid.jmin  ; j < sim->grid.jmax  ; j++)
#define for_each_cells_y(i)            for(u64 i=sim->grid.imin  ; i < sim->grid.imax  ; i++)
#define for_each_cells_and_ghost_x(j)  for(u64 j=0               ; j < sim->grid.Nx_tot; j++)
#define for_each_cells_and_ghost_y(i)  for(u64 i=0               ; i < sim->grid.Ny_tot; i++)
#define for_each_interfaces_x(j)       for(u64 j=sim->grid.imin-1; j < sim->grid.jmax  ; j++)
#define for_each_interfaces_y(i)       for(u64 i=sim->grid.jmin-1; i < sim->grid.imax  ; i++)
#define cell_id(i,j) ((i)*sim->grid.Nx_tot + (j))

#define DECLARE_PSTATE_VAR                  \
        double *rho   = sim->pstate.rho;    \
        double *u     = sim->pstate.u;      \
        double *v     = sim->pstate.v;      \
        double *p     = sim->pstate.p;

#define DECLARE_CSTATE_VAR                  \
        double *rho   = sim->cstate.rho;    \
        double *rho_u = sim->cstate.rho_u;  \
        double *rho_v = sim->cstate.rho_v;  \
        double *E     = sim->cstate.E;
        
#define DECLARE_STATES_VAR                  \
        double *prho  = sim->pstate.rho;    \
        double *u     = sim->pstate.u;      \
        double *v     = sim->pstate.v;      \
        double *p     = sim->pstate.p;      \
        double *crho  = sim->cstate.rho;    \
        double *rho_u = sim->cstate.rho_u;  \
        double *rho_v = sim->cstate.rho_v;  \
        double *E     = sim->cstate.E;
