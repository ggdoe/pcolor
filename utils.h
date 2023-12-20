#pragma once

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

void ring_map(double *x, double *y, int w_grid, int h_grid)
{
    for (int i = 0; i < h_grid * w_grid; i++) {
        double r = 0.5 * (1.0 + y[i]);
        double t = x[i] * 2 * M_PI;
        x[i] = r * cos(t);
        y[i] = r * sin(t);
    }
}

void pert_map(double *x, double *y, int w_grid, int h_grid, double amplitude)
{
    for (int i = 0; i < h_grid * w_grid; i++) {
        double r1 = amplitude * 2.0*((double)rand()/(double)RAND_MAX - 0.5);
        double r2 = amplitude * 2.0*((double)rand()/(double)RAND_MAX - 0.5);
        x[i] += r1;
        y[i] += r2;
    }
}

void colella_map(double *x, double *y, int w, int h)
{
    const double Lx = 1.0;
    const double Ly = 1.0;

    for (int i = 0; i < w*h; i++)
    {
        double sinsin = 0.07 * sin(2.0 * M_PI * x[i] / Lx) * sin(2.0 * M_PI * y[i] / Ly);
        x[i] += sinsin;
        y[i] += sinsin;
    }
}