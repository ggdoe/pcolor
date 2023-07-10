#include <stdint.h>
#include <string.h>
#include <float.h>

#define COLOR_OUT 0xFFFFFFFF

void plot(uint32_t *pixels, int width, int height, double *data, uint32_t nb_elem, uint32_t color)
{
    double min = data[0], max = data[0];
    for(uint32_t k=0; k < nb_elem; k++){
        min = (min < data[k]) ? min : data[k];
        max = (max > data[k]) ? max : data[k];
    }

    // memset(pixels, COLOR_OUT, width * height * sizeof(uint32_t));

    // padding
    {
        double diff = (max - min);
        min -= 0.1 * diff;
        max += 0.1 * diff;
    }


    int i = 0;
    int j = (height - 1) * (1.0 - (data[0] - min) / (max - min));

    for(uint32_t k=1; k < nb_elem; k++){

        int new_i = k * (width - 1) / (nb_elem - 1);
        int new_j = (height - 1) * (1.0 - (data[k] - min) / (max - min));
        
        int jmin = (j < new_j) ? j : new_j;
        int jmax = (j > new_j) ? j : new_j;

        // if(jmax - jmin > new_i - i)
            for(int jj=jmin; jj<jmax; jj++){
                int ii = new_i + (i - new_i) * (new_j - jj)/(new_j - j);
                pixels[jj * width + ii] = color;
            }
        // else
            for(int ii=i; ii<new_i; ii++){
                int jj = new_j + (j - new_j) * (new_i - ii)/(new_i - i);
                pixels[jj * width + ii] = color;
            }
        
        i = new_i;
        j = new_j;
    }
}