#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <time.h>
#include <math.h>

#include "show.h"
#include "plot.h"
#include "cmap.h"
#include "utils.h"

#include "pcolor.h"

#define IMG_FACTOR 60
#define IMG_WIDTH  (16 * IMG_FACTOR)
#define IMG_HEIGHT (16  * IMG_FACTOR)

struct callback_args
{
  int w;
  int h;
  struct pcolor_state *pstate;
  int map;
  double deform_factor;
};

void callback_fill_color(void *args)
{
  static int iter = 0;
  struct callback_args *s = args;

  int w = s->w;
  int h = s->h;
  uint32_t C[(w-1)*(h-1)];

  for (int i = 0; i < (w-1)*(h-1); i++){
      C[i]  = cmap_nipy_spectral(0.4 + 0.5 * sin(2 * M_PI * (double)i/((w-1)*(h-1)) + (double)iter/20.))
            ^ cmap_nipy_spectral(0.2 + 0.3 * cos((double)i/8.0 + (double)iter/100.0));
  }

  pcolor(s->pstate, C);

  iter++;
  iter = iter%400;
}

void callback_cycle_mapping(void *args)
{
  struct callback_args *s = args;

  int w = s->w;
  int h = s->h;
  double x[w * h];
  double y[w * h];
  meshgrid(x, y, w, h, 0.0, 1.0, 0.0, 1.0);

  switch(s->map++ % 5){
    case 0:
      ring_map(x, y, w, h);
      break;
    case 1:
      colella_map(x, y, w, h, s->deform_factor);
      break;
    case 2:
      ring_map(x, y, w, h);
      colella_map(x, y, w, h, s->deform_factor);
      break;
    case 3:
      ring_map(x, y, w, h);
      colella_map(x, y, w, h, s->deform_factor);
      pert_map(x, y, w, h, 0.01);
      // rectify perturbation mapping
      for(int i=0; i < h; i++){
        x[i*w] = x[(i+1)*w - 1];
        y[i*w] = y[(i+1)*w - 1];
      }
      break;
  }
  pcolor_fill_state(s->pstate, x, y, w, h);
}

void callback_toggle_edge(void *args)
{
  struct pcolor_state *pstate = args;
  pstate->show_edge = !pstate->show_edge;
}
void callback_inc_h(void *args)
{
  struct callback_args *s = args;
  s->h+=5;
  s->map--;callback_cycle_mapping(args);
}
void callback_dec_h(void *args)
{
  struct callback_args *s = args;
  s->h-=5;
  if(s->h < 2) s->h=2;
  s->map--;callback_cycle_mapping(args);
}
void callback_inc_w(void *args)
{
  struct callback_args *s = args;
  s->w+=5;
  s->map--;callback_cycle_mapping(args);
}
void callback_dec_w(void *args)
{
  struct callback_args *s = args;
  s->w-=5;
  if(s->w < 2) s->w=2;
  s->map--;callback_cycle_mapping(args);
}
void callback_inc_deform(void *args)
{
  struct callback_args *s = args;
  s->deform_factor+=0.005;
  s->map--;callback_cycle_mapping(args);
}
void callback_dec_deform(void *args)
{
  struct callback_args *s = args;
  s->deform_factor-=0.005;
  s->map--;callback_cycle_mapping(args);
}


int main(int argc, char **argv)
{
  uint32_t *pixels = malloc(IMG_WIDTH * IMG_HEIGHT * sizeof(uint32_t));
  memset(pixels, 0xFFFFFFFF, IMG_WIDTH * IMG_HEIGHT * sizeof(uint32_t));

  // init pcolor state
  struct pcolor_state pstate = pcolor_state_alloc(pixels, IMG_WIDTH, IMG_HEIGHT);
  // pstate.show_edge = false;
  pstate.color_egde = 0xFF000000;

  struct callback_args args = {.h = 20, .w = 80, .map = 0, .deform_factor = 0.01, .pstate = &pstate};
  
  callback_cycle_mapping(&args);
  callback_fill_color(&args);
  
  struct custom_keyevent ckey[8]; const int number_custom_keyevent = sizeof(ckey)/sizeof(struct custom_keyevent);
  ckey[0] = (struct custom_keyevent){.key = SDLK_m, .callback = callback_cycle_mapping, .callback_args = &args};
  ckey[1] = (struct custom_keyevent){.key = SDLK_e, .callback = callback_toggle_edge, .callback_args = &pstate};
  ckey[2] = (struct custom_keyevent){.key = SDLK_UP, .callback = callback_inc_h, .callback_args = &args};
  ckey[3] = (struct custom_keyevent){.key = SDLK_DOWN, .callback = callback_dec_h, .callback_args = &args};
  ckey[4] = (struct custom_keyevent){.key = SDLK_LEFT, .callback = callback_dec_w, .callback_args = &args};
  ckey[5] = (struct custom_keyevent){.key = SDLK_RIGHT, .callback = callback_inc_w, .callback_args = &args};
  ckey[6] = (struct custom_keyevent){.key = SDLK_KP_PLUS, .callback = callback_inc_deform, .callback_args = &args};
  ckey[7] = (struct custom_keyevent){.key = SDLK_KP_MINUS, .callback = callback_dec_deform, .callback_args = &args};
  
  animate(pixels, IMG_WIDTH, IMG_HEIGHT, 60.0, callback_fill_color, &args, ckey, number_custom_keyevent);

  free(pixels);
  pcolor_free(&pstate);
  return 0;
}
