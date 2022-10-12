#include "../all.h"

int dns_eigen(char c, char e, char info_mat, char info_ls,const int mat_dim,const int tri_mat_dim,std::string matfile_name, std::string outfile_name,std::string testfile_name,const int test, int precision){
    std::cout << "Start Matrix Diagonalization (DNS format)\n" ;
    std::cout << "=======================================================================\n"
              << "INPUT FILE : " << matfile_name << std::endl
              << "=======================================================================\n";
    
    std::cout << "Mode : ";
    if(c == 'A'){
        std::cout << "Diagonalization by lanczos method\n";
    }
    else if(c == 'T'){
        std::cout << "Test Mode\n";
    }
std::cout << "=======================================================================\n";
    /*=========================ファイルの存在確認=================================*/
    std::ifstream ifs_matfile(matfile_name);
    if(!ifs_matfile){
        std::cerr << "Could not open the file(line 31) -'" << matfile_name << "'" << std::endl;
    }
    /*=========================================================================*/
    
    double *M = new double[mat_dim * mat_dim]; //行列の要素を格納する
    vec_init(mat_dim*mat_dim,M);

    /*===========================配列に行列要素を代入する=========================*/
    int num = 0;
    double mat_val = 0.; 
    while(ifs_matfile >> mat_val){
        M[num] = mat_val;
        num++;
    }
    ifs_matfile.clear();
    ifs_matfile.seekg(std::ios::beg);

    if(info_mat == 'y'){
        printmat(mat_dim, M);
    }
    /*==========================================================================*/

    /*=============================lacnzos法====================================*/
    if(c == 'A'){
        if(e == 'N') dns_lancozs_N(info_ls, mat_dim, tri_mat_dim,M,outfile_name);
        if(e == 'V') dns_lancozs_V(info_ls, mat_dim, tri_mat_dim,M,outfile_name);
    }
    else if(c=='T'){//テストモード
        std::ofstream ofs_testfile(testfile_name);
        if(!ifs_matfile){
            std::cerr << "Could not open the file(line 31) -'" << testfile_name << "'" << std::endl;
        }   

        double *arr_test = new double[2];
        double *lw = new double[mat_dim];

        if(e == 'N'){
            ofs_testfile << "No." << "," <<"lapack" << "," <<"My ans" << "," << "error\n";
        
            LAPACKE_dsyev(LAPACK_COL_MAJOR, 'N', 'U', mat_dim, M, mat_dim,lw);
            arr_test[0] = lw[0];
            for(int t = 0;t < test; t++){
                arr_test[1] = dns_lancozs_N(info_ls, mat_dim, tri_mat_dim,M,outfile_name);
                ofs_testfile << arr_test[0] << "," << arr_test[1] << "," << abs(arr_test[0]- arr_test[1]) <<std::endl;      
            }
        }
        else if(e == 'V'){
            for(int t = 0;t < test; t++){
                dns_lancozs_V(info_ls, mat_dim, tri_mat_dim,M,outfile_name,c,testfile_name); 
            }
        }
        
        
        delete[] arr_test;
        delete[] lw;
    }
    delete[] M;
}
