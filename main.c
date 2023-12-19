#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <time.h>
#include <math.h>

#include "show.h"
#include "plot.h"
#include "cmap.h"

#include "pcolor.h"

#define IMG_FACTOR 60
#define IMG_WIDTH  (16 * IMG_FACTOR)
#define IMG_HEIGHT (12  * IMG_FACTOR)

int main(int argc, char **argv)
{
    uint32_t *pixels = malloc(IMG_WIDTH * IMG_HEIGHT * sizeof(uint32_t));
    memset(pixels, 0xFFFFFFFF, IMG_WIDTH * IMG_HEIGHT * sizeof(uint32_t));

    #if 0
        #define NB 300
        double value[NB];
        double value2[NB];
        for(int i=0; i < NB; i++)
        {
            double x = 5.*((double)i/(double)(NB-1) - 0.5f);
            value[i] = exp(-x*x);
            value2[i] = sin(-(double)i/NB * 9 * 3.141592);
        }

        fill(pixels, IMG_WIDTH, IMG_HEIGHT, 255, 255, 255);

        struct lim ylim, xlim = {51.0*(0.0 - 0.5f) , 51.0*(1.0 - 0.5f)};
        ylim = compute_lim(value,  NB, NULL);
        ylim = compute_lim(value2, NB, &ylim);

        grid(pixels, IMG_WIDTH, IMG_HEIGHT, &xlim, &ylim, 0xFFAAAAAA);

        plot(pixels, IMG_WIDTH, IMG_HEIGHT, value,  NB, &ylim, 0xFFFF0000);
        plot(pixels, IMG_WIDTH, IMG_HEIGHT, value2, NB, &ylim, 0xFF0000FF);
        
    #elif 0
        for(int k=0; k < IMG_WIDTH * IMG_HEIGHT; k++)
        {
            double i = (double)(k / IMG_HEIGHT) / (double)IMG_WIDTH;
            double j = (double)(k % IMG_WIDTH) / (double)IMG_HEIGHT;
            double v = fabs(sin(4.5*M_PI*i + exp(j*(1-j)*10)));
            pixels[k] = cmap_nipy_spectral(v*(1-v)*4);
        }
    #else
    
        int w = 80;
        int h = 20;
        double x[w * h];
        double y[w * h];
        uint32_t C[(w-1) * (h-1)];
        
        for (int i = 0; i < (w-1) * (h-1); i++)
            C[i] = rand() | 0xFF000000;
            
        meshgrid(x, y, w, h, 0.0, 1.0, 0.0, 1.0);

        ring_map(x, y, w, h);
        pert_map(x, y, w, h, 0.008);

        pcolor(pixels, IMG_WIDTH, IMG_HEIGHT, x, y, w, h, C, 0xFFFFFFFF, 0xFF111111);

    #endif

    show(pixels, IMG_WIDTH, IMG_HEIGHT);

    free(pixels);
    return 0;
}
