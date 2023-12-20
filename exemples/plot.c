#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <time.h>
#include <math.h>

#include "show.h"
#include "plot.h"
#include "cmap.h"
#include "utils.h"

#define IMG_FACTOR 60
#define IMG_WIDTH  (16 * IMG_FACTOR)
#define IMG_HEIGHT (12  * IMG_FACTOR)

int main(int argc, char **argv)
{
    uint32_t *pixels = malloc(IMG_WIDTH * IMG_HEIGHT * sizeof(uint32_t));
    memset(pixels, 0xFFFFFFFF, IMG_WIDTH * IMG_HEIGHT * sizeof(uint32_t));

    #define NB_POINTS 300
    double value[NB_POINTS];
    double value2[NB_POINTS];
    for(int i=0; i < NB_POINTS; i++)
    {
        double x = 5.*((double)i/(double)(NB_POINTS-1) - 0.5f); // use linspace instead
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

    show(pixels, IMG_WIDTH, IMG_HEIGHT);

    free(pixels);
    return 0;
}
