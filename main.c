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
#define IMG_HEIGHT (12  * IMG_FACTOR)
// #define IMG_HEIGHT (16  * IMG_FACTOR)

void callback(void *args)
{
    uint32_t *pixels = args;
    int w = 88;
    int h = 19;
    double x[w * h];
    double y[w * h];
    uint32_t C[(w-1) * (h-1)];
    
    meshgrid(x, y, w, h, 0.0, 1.0, 0.0, 1.0);


    ring_map(x, y, w, h);
    colella_map(x, y, w, h);
    // pert_map(x, y, w, h, 0.01);

    for (int i = 0; i < (w-1) * (h-1); i++){
        C[i] = cmap_nipy_spectral((double)rand()/RAND_MAX);
    }
    pcolor_nostate(pixels, IMG_WIDTH, IMG_HEIGHT, x, y, w, h, C, true, 0xFFFFFFFF, 0xFF111111);
    // printf("inside callback\n");
}

void callback_keyJ(void *args)
{
    printf("touche J appuyer !!\n");
}

int main(int argc, char **argv)
{
    uint32_t *pixels = malloc(IMG_WIDTH * IMG_HEIGHT * sizeof(uint32_t));
    memset(pixels, 0xFFFFFFFF, IMG_WIDTH * IMG_HEIGHT * sizeof(uint32_t));

    #if 0
        #define NB_POINTS 300
        double value[NB_POINTS];
        double value2[NB_POINTS];
        for(int i=0; i < NB_POINTS; i++)
        {
            double x = 5.*((double)i/(double)(NB_POINTS-1) - 0.5f);
            value[i] = exp(-x*x);
            value2[i] = sin(-(double)i/NB_POINTS * 9 * 3.141592);
        }

        fill(pixels, IMG_WIDTH, IMG_HEIGHT, 255, 255, 255);

        struct lim xlim = {5.0*(0.0 - 0.5f) , 5.0*(1.0 - 0.5f)};
        struct lim ylim;
        ylim = compute_lim(value,  NB_POINTS, NULL);
        ylim = compute_lim(value2, NB_POINTS, &ylim);

        fill_grid(pixels, IMG_WIDTH, IMG_HEIGHT, &xlim, &ylim, 0xFFAAAAAA);

        plot(pixels, IMG_WIDTH, IMG_HEIGHT, value,  NB_POINTS, &ylim, 0xFFFF0000);
        plot(pixels, IMG_WIDTH, IMG_HEIGHT, value2, NB_POINTS, &ylim, 0xFF0000FF);
        
    #elif 0
        for(int k=0; k < IMG_WIDTH * IMG_HEIGHT; k++)
        {
            double i = (double)(k / IMG_HEIGHT) / (double)IMG_WIDTH;
            double j = (double)(k % IMG_WIDTH) / (double)IMG_HEIGHT;
            double v = fabs(sin(4.5*M_PI*i + exp(j*(1-j)*10)));
            pixels[k] = cmap_nipy_spectral(v*(1-v)*4);
        }
    #else
    
        int w = 88;
        int h = 19;
        double x[w * h];
        double y[w * h];
        uint32_t C[(w-1) * (h-1)];
        
        meshgrid(x, y, w, h, 0.0, 1.0, 0.0, 1.0);


        ring_map(x, y, w, h);
        colella_map(x, y, w, h);
        // pert_map(x, y, w, h, 0.01);

        for (int i = 0; i < (w-1) * (h-1); i++){
            C[i] = cmap_nipy_spectral((double)rand()/RAND_MAX);
        }

        struct pcolor_state pstate = pcolor_state_alloc(pixels, IMG_HEIGHT, IMG_HEIGHT);
        pstate.ld_pixels = IMG_WIDTH;
        pcolor_fill_state(&pstate, x, y, w, h);
        pstate.show_edge = false;
        pcolor(&pstate, C);
        // pcolor_nostate(pixels, IMG_WIDTH, IMG_HEIGHT, x, y, w, h, C, true, 0xFFFFFFFF, 0xFF111111);

    #endif

    // show(pixels, IMG_WIDTH, IMG_HEIGHT);
    
    struct custom_keyevent ckey = {.key = SDLK_j, .callback = callback_keyJ, .callback_args = NULL};
    animate(pixels, IMG_WIDTH, IMG_HEIGHT, 60.0, callback, pixels, &ckey, 1);

    free(pixels);
    pcolor_free(&pstate);
    return 0;
}
