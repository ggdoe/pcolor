#pragma once

#include <stdlib.h>
#include <math.h>
#include <stdint.h>

void meshgrid(double *x, double *y, int w, int h, double xmin, double xmax, double ymin, double ymax) 
{
    double dx = (xmax - xmin) / (w - 1);
    double dy = (ymax - ymin) / (h - 1);

    int i, j;
    for (i = 0; i < w; i++) {
        for (j = 0; j < h; j++) {
            int index = j * w + i;
            x[index] = xmin + i * dx;
            y[index] = ymin + j * dy;
        }
    }
}

void mapc2p(double *x, double *y, int w, int h)
{
    const double Lx = 1.0;
    const double Ly = 1.0;

    for (int i = 0; i < w*h; i++)
    {
        int sinsin = 0.75 * sin(M_PI * x[i] / Lx) * sin(M_PI * y[i] / Ly);
        x[i] += sinsin;
        y[i] += sinsin;
    }
}

void draw_cross(uint32_t *pixels, int w_pixels, int h_pixels, int i, int j, uint32_t color)
{
    const int w_cross = 5;

    int imin = (i - w_cross     < 0)        ?        0 : i - w_cross;
    int imax = (i + w_cross + 1 > h_pixels) ? h_pixels : i + w_cross + 1;
    int jmin = (j - w_cross     < 0)        ?        0 : j - w_cross;
    int jmax = (j + w_cross + 1 > w_pixels) ? w_pixels : j + w_cross + 1;

    for(int jj = jmin; jj < jmax; jj++)
        pixels[i * w_pixels + jj] = color;
    for(int ii = imin; ii < imax; ii++)
        pixels[ii * w_pixels + j] = color;
}

// void draw_line(uint32_t *pixels, int w_pixels, int h_pixels, int i0, int j0, int i1, int j1, uint32_t color)
// {
//     // Bresenham's line algorithm
//     const int imin = (i0 < i1) ? i0 : i1;
//     const int imax = (i0 > i1) ? i0 : i1;
//     const int jmin = (j0 < j1) ? j0 : j1;
//     const int jmax = (j0 > j1) ? j0 : j1;

//     const int di = (imax - imin);
//     const int dj = (jmax - jmin);

//     int i = imin;
//     int D = 2*di - dj;

//     for(int j = jmin; j <= jmax; j++)
//     {
//         pixels[i*w_pixels + j] = color;
//         if(D > 0){
//             i++;
//             D -= 2*dj;
//         }
//         D += 2*di;
//     }
// }
void draw_line(uint32_t *pixels, int w_pixels, int h_pixels, int i0, int j0, int i1, int j1, uint32_t color)
{
    // Bresenham's line algorithm
    const int di = -abs(i1 - i0);
    const int dj =  abs(j1 - j0);
    const int si = (i0 < i1) ? 1 : -1;
    const int sj = (j0 < j1) ? 1 : -1;
    int error = di + dj;

    while(true)
    {
        pixels[i0*w_pixels + j0] = color;
        if(i0 == i1 && j0 == j1) break;
        if(2*error >= di){
            if(j0 == j1) break;
            error += di;
            j0 += sj;
        }
        if(2*error <= dj){
            if(i0 == i1) break;
            error += dj;
            i0 += si;
        }
    }
}

// void draw_line(uint32_t *pixels, int w_pixels, int h_pixels, int y0, int x0, int y1, int x1, uint32_t color)
// {
//     int dx = abs(x1-x0), sx = x0<x1 ? 1 : -1;
//     int dy = -abs(y1-y0), sy = y0<y1 ? 1 : -1;
//     int err = dx+dy, e2; /* error value e_xy */

//     for (;;){ /* loop */
//         pixels[y0*w_pixels + x0] = color;
//         e2 = 2*err;
//         if (e2 >= dy) { /* e_xy+e_x > 0 */
//             if (x0 == x1) break;
//             err += dy; x0 += sx;
//         }
//         if (e2 <= dx) { /* e_xy+e_y < 0 */
//             if (y0 == y1) break;
//             err += dx; y0 += sy;
//         }
//     }
// }


static double _pcolor_min(const double *x, int n)
{
    double min = x[0];
    for(int i = 1; i < n; i++)
        min = (x[i] < min) ? x[i] : min;
    return min;
}
static double _pcolor_max(const double *x, int n)
{
    double max = x[0];
    for(int i = 1; i < n; i++)
        max = (x[i] > max) ? x[i] : max;
    return max;
}

static int _pcolor_min_int(const int *x, int n)
{
    int min = x[0];
    for(int i = 1; i < n; i++)
        min = (x[i] < min) ? x[i] : min;
    return min;
}
static int _pcolor_max_int(const int *x, int n)
{
    int max = x[0];
    for(int i = 1; i < n; i++)
        max = (x[i] > max) ? x[i] : max;
    return max;
}

// static double _pcolor_argmin(const double *x, int n)
// {
//     int min = x[0];
//     int arg = 0;
//     for(int i = 1; i < n; i++)
//         if(x[i] < min){
//             min = x[i];
//             arg = i;
//         }
//     return arg;
// }
// static double _pcolor_argmax(const double *x, int n)
// {
//     int max = x[0];
//     int arg = 0;
//     for(int i = 1; i < n; i++)
//         if(x[i] > max){
//             max = x[i];
//             arg = i;
//         }
//     return arg;
// }

void ring_map(double *x, double *y, int w_grid, int h_grid)
{
    // const double xmin = _pcolor_min(x, w_grid*h_grid), ymin = _pcolor_min(y, w_grid*h_grid);
    // const double xmax = _pcolor_max(x, w_grid*h_grid), ymax = _pcolor_max(y, w_grid*h_grid);

    for (int i = 0; i < h_grid * w_grid; i++) {
        double r = 0.5 * (1.0 + y[i]);
        double t = x[i] * 2 * M_PI;
        x[i] = r * cos(t);
        y[i] = r * sin(t);
    }
}

void pert_map(double *x, double *y, int w_grid, int h_grid, double amplitude)
{
    srand(0);
    for (int i = 0; i < h_grid * w_grid; i++) {
        double r1 = amplitude * 2.0*((double)rand()/(double)RAND_MAX - 0.5);
        double r2 = amplitude * 2.0*((double)rand()/(double)RAND_MAX - 0.5);
        x[i] += r1;
        y[i] += r2;
    }
}

bool is_inside_quad(int *iv, int *jv, int i, int j)
{
    // if sign is always the same --> point is inside
    bool sign  = (iv[0]-iv[3])*(j-jv[3]) >= (jv[0]-jv[3])*(i-iv[3]);

    for(int id = 0; id < 3; id++){
        int i0 = iv[id], i1 = iv[id+1];
        int j0 = jv[id], j1 = jv[id+1];

        if( sign ^ ((j-j0)*(i1-i0) >= (i-i0)*(j1-j0)))
            return false;
    }
    return true;
}

// void check_fct(uint32_t *pixels, int w_pixels, int h_pixels)
// {
    
//     int jv[] = {100, 400, 400, 100, jv[0]};
//     int iv[] = {100, 100, 400, 400, iv[0]};
//     jv[1] -= 50;
//     iv[1] -= 50;

//     for (int i = 0; i < h_pixels; i++)
//         for (int j = 0; j < w_pixels; j++)
//             pixels[i * w_pixels + j] = is_inside_quad(iv, jv, i, j) ? 0x88888888 : 0;
            
//     // printf("%d\n", );
//     for (int i = 0; i < 4; i++)
//         draw_cross(pixels, w_pixels, h_pixels, iv[i], jv[i], 0xFFFFFFFF);
// }

void pcolor(uint32_t *pixels, int w_pixels, int h_pixels, double *x, double *y, int w_grid, int h_grid, uint32_t *color_grid, uint32_t color_egde, uint32_t color_outside) 
{
    for (int i = 0; i < w_pixels * h_pixels; i++)
        pixels[i] = color_outside;

    const double margin = 0.0;
    const double xmin = _pcolor_min(x, w_grid*h_grid) - margin, ymin = _pcolor_min(y, w_grid*h_grid) - margin;
    const double xmax = _pcolor_max(x, w_grid*h_grid) + margin, ymax = _pcolor_max(y, w_grid*h_grid) + margin;

    for (int i = 0; i < h_grid - 1; i++) {
        for (int j = 0; j < w_grid - 1; j++) {

            // vertex of the cell
            int jv[] = {
                (x[ i    * w_grid + j  ] - xmin) / (xmax - xmin) * (w_pixels-1), 
                (x[ i    * w_grid + j+1] - xmin) / (xmax - xmin) * (w_pixels-1), 
                (x[(i+1) * w_grid + j+1] - xmin) / (xmax - xmin) * (w_pixels-1), 
                (x[(i+1) * w_grid + j  ] - xmin) / (xmax - xmin) * (w_pixels-1)
            };
            int iv[] = {
                (y[ i    * w_grid + j  ] - ymin) / (ymax - ymin) * (h_pixels-1), 
                (y[ i    * w_grid + j+1] - ymin) / (ymax - ymin) * (h_pixels-1), 
                (y[(i+1) * w_grid + j+1] - ymin) / (ymax - ymin) * (h_pixels-1),
                (y[(i+1) * w_grid + j  ] - ymin) / (ymax - ymin) * (h_pixels-1)
            };

            // boundary of the box
            int jmin = _pcolor_min_int(jv, 4);
            int jmax = _pcolor_max_int(jv, 4);
            int imin = _pcolor_min_int(iv, 4);
            int imax = _pcolor_max_int(iv, 4);

            for(int id = imin; id < imax; id++)
                for(int jd = jmin; jd < jmax; jd++)
                    if(is_inside_quad(iv, jv, id, jd))
                        pixels[id * w_pixels + jd] =  color_grid[i * w_grid + j];
                    
            draw_line(pixels, w_pixels, h_pixels, iv[0], jv[0], iv[1], jv[1], color_egde);
            draw_line(pixels, w_pixels, h_pixels, iv[0], jv[0], iv[3], jv[3], color_egde);
            draw_line(pixels, w_pixels, h_pixels, iv[3], jv[3], iv[2], jv[2], color_egde);
            draw_line(pixels, w_pixels, h_pixels, iv[1], jv[1], iv[2], jv[2], color_egde);
        }
    }
    // for (int i = 0; i < h_grid; i++) {
    //     for (int j = 0; j < w_grid; j++) {
    //         int ii = (y[i * w_grid + j]     - ymin) / (ymax - ymin) * (h_pixels-1);
    //         int jj = (x[i * w_grid + j]     - xmin) / (xmax - xmin) * (w_pixels-1);
    //         draw_cross(pixels, w_pixels, h_pixels, ii, jj, 0x00FF0000);
    //     }
    // }
}