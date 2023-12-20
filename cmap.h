#pragma once
#include <math.h>
#include <stdint.h>

uint32_t cmap_nipy_spectral(double v)
{
    static float data[22][3] = 
            {   
                {  0.,   0.,   0.},
                {119.,   0., 136.},
                {136.,   0., 153.},
                {  0.,   0., 170.},
                {  0.,   0., 221.},
                {  0., 119., 221.},
                {  0., 153., 221.},
                {  0., 170., 170.},
                {  0., 170., 136.},
                {  0., 153.,   0.},
                {  0., 187.,   0.},
                {  0., 221.,   0.},
                {  0., 255.,   0.},
                {187., 255.,   0.},
                {238., 238.,   0.},
                {255., 204.,   0.},
                {255., 153.,   0.},
                {255.,   0.,   0.},
                {221.,   0.,   0.},
                {204.,   0.,   0.},
                {204., 204., 204.},
                {204., 204., 204.}
            };

    v = (v < 0.0) ? 0.0 : (v > 1.0) ? 1.0 : v;
    double index = v * 20.;
    double frac = index - floor(index);
    size_t id = index;

    uint32_t r = ((1-frac) * data[id][0] + frac * data[id+1][0]);
    uint32_t g = ((1-frac) * data[id][1] + frac * data[id+1][1]);
    uint32_t b = ((1-frac) * data[id][2] + frac * data[id+1][2]);

    return (r<<16) | (g<<8) | b;
}
