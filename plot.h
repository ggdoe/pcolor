#include <stdint.h>
#include <string.h>
#include <float.h>

#define COLOR_OUT 0xFFFFFFFF
#define COLOR(r,g,b) (0xFF<<24 | r<<16 | g<<8 | b)

struct lim{
    double min;
    double max;
};

struct lim compute_lim(const double *data, uint32_t nb_elem, const struct lim *old_lim)
{
    struct lim lim;
    if(old_lim != NULL)
        lim = *old_lim;
    else{
        lim.min = data[0];
        lim.max = data[0];
    }

    for(uint32_t k=1; k < nb_elem; k++){
        lim.min = (lim.min < data[k]) ? lim.min : data[k];
        lim.max = (lim.max > data[k]) ? lim.max : data[k];
    }

    return lim;
}

void fill(uint32_t *pixels, int width, int height, uint8_t red, uint8_t green, uint8_t blue)
{
    uint32_t color = 0xFF000000 | red<<16 | green<<8 | blue;
    memset(pixels, color, width * height * sizeof(uint32_t));
}

void plot(uint32_t *pixels, int width, int height, double *data, uint32_t nb_elem, struct lim *ylim, uint32_t color)
{
    double min, max;
    if(ylim == NULL){
        struct lim lim = compute_lim(data, nb_elem, NULL);
        min = lim.min;
        max = lim.max;
    }
    else {
        min = ylim->min;
        max = ylim->max;
    }
    
    // const double padding = 0.02; // doesnt work with grid()
    // {
    //     double diff = (max - min);
    //     min -= padding * diff;
    //     max += padding * diff;
    // }

    int i = 0;
    int j = (height - 1) * (1.0 - (data[0] - min) / (max - min));

    for(uint32_t k=1; k < nb_elem; k++){

        int new_i = k * (width - 1) / (nb_elem - 1);
        int new_j = (height - 1) * (1.0 - (data[k] - min) / (max - min));
        
        int jmin = (j < new_j) ? j : new_j;
        int jmax = (j > new_j) ? j : new_j;

            for(int jj=jmin; jj<jmax; jj++){
                int ii = new_i + (i - new_i) * (new_j - jj)/(new_j - j);
                pixels[jj * width + ii] = color;
            }
            for(int ii=i; ii<new_i; ii++){
                int jj = new_j + (j - new_j) * (new_i - ii)/(new_i - i);
                pixels[jj * width + ii] = color;
            }
        
        i = new_i;
        j = new_j;
    }
}

void grid(uint32_t *pixels, int width, int height, struct lim *xlim, struct lim *ylim, uint32_t color)
{
    const uint32_t axis_color = 0xFF000000;
    double diff_x, dx;
    double diff_y, dy;

    {
        const double count = 4.0;
        const double margin_factor = 1.3;
        diff_x = xlim->max - xlim->min;
        diff_y = ylim->max - ylim->min;
        dx = pow(10, ceil(log10(margin_factor * diff_x) - 1.0))/count;
        dy = pow(10, ceil(log10(margin_factor * diff_y) - 1.0))/count;
    }
    printf("  [grid]  xmin: %lf  xmax: %lf\n", xlim->min, xlim->max);
    printf("  [grid]  dx:   %lf    dy: %lf\n", dx, dy);

    // grid
    for(double x = ceil(xlim->min/dx)*dx; x < xlim->max; x+=dx){
        uint32_t i = width * (x - xlim->min) / diff_x;
        for(uint32_t j = 0; j < height; j++)
            pixels[j * width + i] = color;
    }
    for(double y = ceil(ylim->min/dy)*dy; y < ylim->max; y+=dy){
        uint32_t j = height * (y - ylim->min) / diff_y;
        for(uint32_t i = 0; i < width; i++)
            pixels[j * width + i] = color;
    }

    // axis
    if(xlim->min <= 0.0 && xlim->max >= 0.0){
        uint32_t i = width * (0.0 - xlim->min) / diff_x;
        for(uint32_t j = 0; j < height; j++)
            pixels[j * width + i] = axis_color;
    }
    if(ylim->min <= 0.0 && ylim->max >= 0.0){
        uint32_t j = height * (0.0 - ylim->min) / diff_y;
        for(uint32_t i = 0; i < width; i++)
            pixels[j * width + i] = axis_color;
    }

}
