#ifndef __MYlIB_H_
#define __MYlIB_H_

#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>

#include "cblas.h"
#include "lapacke.h"


int dns_eigen(char c, char e, char info_mat, char info_ls,const int mat_dim,const int tri_mat_dim,std::string matfile_name, std::string outfile_name,std::string testfile_name,const int test=1, int precision=8);
double dns_lancozs_N(char info_ls, const int mat_dim, int tri_mat_dim, double *M, std::string outfile_name);
double dns_lancozs_V(char info_ls, const int mat_dim, int tri_mat_dim, double *M, std::string outfile_name, char c = 'A', std::string testfile_name = "test/test.txt");
void vec_init(int n, double* vec);
void vec_init(int n, int* vec);
void sdz(int mat_dim,double* vec, double err = 1.0e-15);

/*行列出力用の関数テンプレート*/
template <typename T>
void printvec(int vec_elements,T *v, int precision = 5)
{
    T vtmp;
    std::cout << "[";
    for (int col_num = 0; col_num < vec_elements; col_num++)
    {
        vtmp = v[col_num];
        std::cout << std::setw(7) << std::scientific << std::setprecision(precision) << std::left << vtmp;
        if (col_num < vec_elements - 1)
        {
            std::cout << "  ";
        }
    }
    std::cout << "]" << std::endl;
}


template <typename T>
void printmat(int mat_dim, T *A, int precision=5)
{
    T mtmp;
    for (int row_num = 0; row_num < mat_dim; row_num++)
    {
        std::cout << "[";
        for (int col_num = 0; col_num < mat_dim; col_num++)
        {
            mtmp = A[col_num + mat_dim * row_num];
            std::cout << std::scientific << std::setprecision(precision) << mtmp;
            if (col_num < mat_dim - 1)
            {
                std::cout << "  ";
            }
        }
        if (row_num < mat_dim - 1)
        {
            std::cout << "];" << std::endl;
        }
        else
        {
            std::cout << "]";
        }
    }
    std::cout << "]" << std::endl;
}

template <typename T>
void print_tri_diag_vec(int mat_dim, T *diag, T *sub_diag, int precision=5)
{
    T dtmp, subdtmp;
    for (int i = 0; i < mat_dim; i++)
    {
        for (int j = 0; j <= i - 2; j++)
        {
            std::cout << std::setw(7) << std::scientific << std::setprecision(precision) << std::left << 0.0 << "  ";
        }
        if (i == 0)
        {
            std::cout << std::setw(7) << std::scientific << std::setprecision(precision) << std::left << diag[i] << "  " << sub_diag[i] << "  ";
        }
        else if (i == mat_dim - 1)
        {
            std::cout << std::setw(7) << std::scientific << std::setprecision(precision) << std::left << sub_diag[i] << "  " << diag[i] << "  ";
        }
        else
        {
            std::cout << std::setw(7) << std::scientific << std::setprecision(precision) << std::left << sub_diag[i - 1] << "  " << diag[i] << "  "
                      << sub_diag[i] << "  ";
        }
        for (int k = i + 2; k < mat_dim; k++)
        {
            if (k > 0)
            {
                std::cout << std::setw(7) << std::scientific << std::setprecision(precision) << std::left << 0.0 << "  ";
            }
        }
        std::cout << std::endl;
    }
}
#endif
