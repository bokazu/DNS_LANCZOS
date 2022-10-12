#include "../all.h"

void vec_init(int n, double* vec)
{
    for (int i = 0; i < n; i++)
    {
        vec[i] = 0.0;
    }
}

void vec_init(int n, int* vec)
{
    for (int i = 0; i < n; i++)
    {
        vec[i] = 0;
    }
}

void sdz(int mat_dim, double* vec, double err)
{
    double a = 1. / cblas_dnrm2(mat_dim, vec, 1);
    cblas_dscal(mat_dim, a, vec, 1);

    // Check norm
    double norm = cblas_dnrm2(mat_dim, vec, 1);
    double eps = abs(1.0 - norm);
    if (eps > err)
    {
        printf("\x1b[31m NORM ERROR\x1b[39m\n");
    }
}
