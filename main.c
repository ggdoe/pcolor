#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <time.h>
#include <math.h>

#include "show.h"
#include "plot.h"
#include "cmap.h"

#define WIN_FACTOR 60
#define WIN_WIDTH  (16 * WIN_FACTOR)
#define WIN_HEIGHT (12  * WIN_FACTOR)

int main(int argc, char **argv)
{
    uint32_t *pixels = malloc(WIN_WIDTH * WIN_HEIGHT * sizeof(uint32_t));
    memset(pixels, 0xFFFFFFFF, WIN_WIDTH * WIN_HEIGHT * sizeof(uint32_t));

    #if 1
        #define NB 300
        double value[NB];
        double value2[NB];
        for(int i=0; i < NB; i++)
        {
            double x = 10.*((double)i/(double)(NB-1) - 0.5f);
            value[i] = exp(-x*x);
            value2[i] = sin(-(double)i/NB * 9 * 3.141592);
        }

        fill(pixels, WIN_WIDTH, WIN_HEIGHT, 255, 255, 255);

        struct lim ylim, xlim = {10.0*(0.0 - 0.5f) , 10.0*(1.0 - 0.5f)};
        ylim = compute_lim(value,  NB, NULL);
        ylim = compute_lim(value2, NB, &ylim);

        grid(pixels, WIN_WIDTH, WIN_HEIGHT, &xlim, &ylim, 0xFFAAAAAA);

        plot(pixels, WIN_WIDTH, WIN_HEIGHT, value,  NB, &ylim, 0xFFFF0000);
        plot(pixels, WIN_WIDTH, WIN_HEIGHT, value2, NB, &ylim, 0xFF0000FF);
        
    #else
        for(int k=0; k < WIN_WIDTH * WIN_HEIGHT; k++)
        {
            double i = (double)(k / WIN_HEIGHT) / (double)WIN_WIDTH;
            double j = (double)(k % WIN_WIDTH) / (double)WIN_HEIGHT;
            double v = fabs(sin(2*M_PI*i + exp(j*3)));
            pixels[k] = cmap_nipy_spectral(v);
        }
    #endif

    show(pixels, WIN_WIDTH, WIN_HEIGHT);

    free(pixels);
    return 0;
}
