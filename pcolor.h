#pragma once

#include <stdlib.h>
#include <math.h>
#include <stdint.h>

struct pcolor_state{
    int w_pixels;
    int h_pixels;
    int ld_pixels;

    uint32_t *pixels;
    uint32_t *pixel_color_id;
    enum pixel_type {
        PCOLOR_TYPE_CELL,
        PCOLOR_TYPE_EDGE,
        PCOLOR_TYPE_OUTSIDE
    } *pixel_type;

    uint32_t color_egde;
    uint32_t color_outside;
    bool show_edge;
};

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

bool is_inside_quad(int *iv, int *jv, int i, int j)
{
    bool sign  = (iv[0]-iv[3])*(j-jv[3]) > (jv[0]-jv[3])*(i-iv[3]);

    for(int id = 0; id < 3; id++){
        int i0 = iv[id], i1 = iv[id+1];
        int j0 = jv[id], j1 = jv[id+1];

        // idk, it works enough
        if( (!sign && ((j-j0)*(i1-i0) > (i-i0)*(j1-j0))) ^
            ( sign && ((j-j0)*(i1-i0) < (i-i0)*(j1-j0))) )
            return false;
    }
    return true;
}

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

struct pcolor_state pcolor_state_alloc(uint32_t *pixels, int w_pixels, int h_pixels)
{
    struct pcolor_state state = {
        .pixels = pixels,
        .w_pixels = w_pixels,
        .h_pixels = h_pixels,
        .ld_pixels = w_pixels,
        .color_egde = 0xFFFFFFFF,
        .color_outside = 0xFF111111,
        .show_edge = true
    };
    state.pixel_color_id = realloc(state.pixel_color_id, w_pixels*h_pixels*sizeof(uint32_t));
    state.pixel_type     = realloc(state.pixel_type    , w_pixels*h_pixels*sizeof(enum pixel_type));
    return state;
}

void pcolor_free(struct pcolor_state *state)
{
    free(state->pixel_type);
    free(state->pixel_color_id);
}

void pcolor_fill_state(struct pcolor_state *state, double *x, double *y, int w_grid, int h_grid)
{

    const int w_pixels  = state->w_pixels;
    const int h_pixels  = state->h_pixels;

    for (int i = 0; i < w_pixels * h_pixels; i++){
        state->pixel_color_id[i] = 0;
        state->pixel_type[i] = PCOLOR_TYPE_OUTSIDE;
    }

    const double xmin = _pcolor_min(x, w_grid*h_grid), ymin = _pcolor_min(y, w_grid*h_grid);
    const double xmax = _pcolor_max(x, w_grid*h_grid), ymax = _pcolor_max(y, w_grid*h_grid);
    
    for (int i = 0; i < h_grid - 1; i++) {
        for (int j = 0; j < w_grid - 1; j++) {

            // vertex of the current cell
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

            for(int id = imin; id <= imax; id++)
                for(int jd = jmin; jd <= jmax; jd++)
                    if(is_inside_quad(iv, jv, id, jd)){
                        state->pixel_type[id * w_pixels + jd] = PCOLOR_TYPE_CELL;
                        state->pixel_color_id[id * w_pixels + jd] = i * (w_grid-1) + j;
                    }

            draw_line(state->pixel_type, w_pixels, h_pixels, iv[0], jv[0], iv[1], jv[1], PCOLOR_TYPE_EDGE);
            draw_line(state->pixel_type, w_pixels, h_pixels, iv[0], jv[0], iv[3], jv[3], PCOLOR_TYPE_EDGE);
            draw_line(state->pixel_type, w_pixels, h_pixels, iv[3], jv[3], iv[2], jv[2], PCOLOR_TYPE_EDGE);
            draw_line(state->pixel_type, w_pixels, h_pixels, iv[1], jv[1], iv[2], jv[2], PCOLOR_TYPE_EDGE);
        }
    }
}

void pcolor(struct pcolor_state *state, uint32_t *color_grid)
{
    for(int i = 0; i < state->h_pixels; i++)
        for(int j = 0; j < state->w_pixels; j++)
        {
            const int id_pixels = i*state->ld_pixels + j;
            const int id_color  = i*state-> w_pixels + j;
            switch (state->pixel_type[id_color]){
                case PCOLOR_TYPE_CELL:
                    state->pixels[id_pixels] = color_grid[state->pixel_color_id[id_color]];
                    break;
                case PCOLOR_TYPE_EDGE:
                    state->pixels[id_pixels] = (state->show_edge) ? state->color_egde : color_grid[state->pixel_color_id[id_color]];
                    break;
                case PCOLOR_TYPE_OUTSIDE:
                    state->pixels[id_pixels] = state->color_outside;
                    break;
            }
        }
}

void pcolor_nostate(uint32_t *pixels, int w_pixels, int h_pixels, double *x, double *y, int w_grid, int h_grid, uint32_t *color_grid, bool show_edge, uint32_t color_egde, uint32_t color_outside) 
{
    for (int i = 0; i < w_pixels * h_pixels; i++)
        pixels[i] = color_outside;

    const double xmin = _pcolor_min(x, w_grid*h_grid), ymin = _pcolor_min(y, w_grid*h_grid);
    const double xmax = _pcolor_max(x, w_grid*h_grid), ymax = _pcolor_max(y, w_grid*h_grid);

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

            for(int id = imin; id <= imax; id++)
                for(int jd = jmin; jd <= jmax; jd++)
                    if(is_inside_quad(iv, jv, id, jd))
                        pixels[id * w_pixels + jd] =  color_grid[i * (w_grid-1) + j];
            
            if(show_edge){
                draw_line(pixels, w_pixels, h_pixels, iv[0], jv[0], iv[1], jv[1], color_egde);
                draw_line(pixels, w_pixels, h_pixels, iv[0], jv[0], iv[3], jv[3], color_egde);
                draw_line(pixels, w_pixels, h_pixels, iv[3], jv[3], iv[2], jv[2], color_egde);
                draw_line(pixels, w_pixels, h_pixels, iv[1], jv[1], iv[2], jv[2], color_egde);
            }
        }
    }
}