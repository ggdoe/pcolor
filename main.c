#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <time.h>

#include "show.h"
#include "plot.h"
#include "cmap.h"

#define WIN_WIDTH  400
#define WIN_HEIGHT 600

int main(int argc, char **argv)
{
    uint32_t pixels[WIN_WIDTH * WIN_HEIGHT];
    memset(pixels, 0xFFFFFFFF, WIN_WIDTH * WIN_HEIGHT * sizeof(uint32_t));

    init_window(WIN_WIDTH, WIN_HEIGHT);

    #if 1
    #define NB 3000
    double value[NB];
    double value2[NB];
    for(int i=0; i < NB; i++)
    {
        double x = 5.*((double)i/(double)NB - 0.5f);
        value[i] = exp(-x*x);
        value2[i] = sin(-(double)i/NB * M_PI);
    }

    plot(pixels, WIN_WIDTH, WIN_HEIGHT, value, NB, 0xFFFF0000);
    plot(pixels, WIN_WIDTH, WIN_HEIGHT, value2, NB, 0xFF0000FF);
        
    #else
    for(int i=0; i < WIN_WIDTH * WIN_HEIGHT; i++)
        pixels[i] = cmap_nipy_spectral(sin((double)i / WIN_WIDTH / WIN_HEIGHT)) | 0xFF000000;
    #endif

    show(pixels, WIN_WIDTH, WIN_HEIGHT);

    close_window();
    return 0;
}
