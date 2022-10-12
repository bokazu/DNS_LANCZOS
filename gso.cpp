#include "all.h"

using namespace std;

void gso(int n, int k, double** u)
{
    double cef = 0.;
    for (int try_num = 0; try_num < k - 1; try_num++)
    {
        cef = -1.0 * cblas_ddot(n, u[try_num], 1, u[k + 1], 1);
        cblas_daxpy(n, cef, u[try_num], 1, u[k + 1], 1);
    }
}

int sdz(int n, double err,double* v)
{
    double a = 0.;
    int info_sdz = 0;//0 : ok , 1 : error occured
    a = 1. / cblas_dnrm2(n, v, 1);
    cblas_dscal(n, a, v, 1);

    // Check norm
    double norm = 0.;
    norm = cblas_dnrm2(n, v, 1);
    double eps = abs(1.0 - norm);
    if (eps > 1.0e-15)
    {
        info_sdz = 1;
        printf("\x1b[31m NORM ERROR\x1b[39m\n");
    }
    return info_sdz;
}

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