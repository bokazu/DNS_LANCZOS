#include "../all.h"

double dns_lancozs_N(char info_ls, const int mat_dim, int tri_mat_dim, double *M, std::string outfile_name)
{
    int count = 0;
    double eps = 1.0;         //誤差を代入するようの変数
    double err = 1.0e-16;     //要求精度
    bool err_checker = true;  //誤差が要求精度の範囲内に収まっているかを確認するためのflag

    /*==========================初期ベクトルの用意=========================*/
    double **u = new double *[2];
    for(int k=0;k<2;k++){
        u[k] = new double[mat_dim];
    }

    std::random_device rand;
    std::mt19937 mt(rand());
    std::uniform_real_distribution<> rand1(0,1);
    for(int k=0;k<mat_dim;k++){
        u[0][k] = rand1(mt);
        u[1][k] = 0.0;
    }
    sdz(mat_dim,u[0]);
    /*=================================================================*/
    double *alpha = new double[tri_mat_dim];
    vec_init(tri_mat_dim,alpha);
    
    double *beta = new double[tri_mat_dim - 1];
    vec_init(tri_mat_dim - 1,beta);
    
    double *eval_even = new double[tri_mat_dim];
    vec_init(tri_mat_dim, eval_even);

    double *eval_odd = new double[tri_mat_dim];
    vec_init(tri_mat_dim, eval_odd);

    double *diag = new double[tri_mat_dim];
    vec_init(tri_mat_dim, diag);

    double *sub_diag = new double[tri_mat_dim - 1];
    vec_init(tri_mat_dim - 1, sub_diag);

    double *tri_diag_evec;//lapack_dstevに渡す用('N'の場合は参照されないので適当に用意しておく)
    // vec_init(tri_mat_dim*tri_mat_dim, tri_diag_evec);



    /*==========================lanczos法===============================*/
    for(int ls = 0;ls < tri_mat_dim; ls++){
        if(err_checker){
            count = ls;
            /*========メモリ削減のためのlanczosベクトルの更新===========*/
            if(ls > 0){
                if(ls % 2 == 0){
                    cblas_dscal(mat_dim, -beta[ls-1],u[1],1);
                }
                else{
                    cblas_dscal(mat_dim, -beta[ls-1],u[0],1);
                }
            }
            /*====================================================*/


            /*=============alpha, betaの計算とlancozsベクトルの更新==========*/
            if(ls % 2 == 0){
                cblas_dgemv(CblasColMajor, CblasNoTrans, mat_dim , mat_dim , 1.0 , M , mat_dim, u[0], 1 , 1.0 , u[1] , 1);
                alpha[ls] = cblas_ddot(mat_dim,u[1],1,u[0],1);
                cblas_daxpy(mat_dim,-alpha[ls], u[0], 1, u[1],1);
                if(ls != tri_mat_dim -1){
                    beta[ls] = cblas_dnrm2(mat_dim, u[1], 1);
                    cblas_dscal(mat_dim, 1./beta[ls],u[1],1);
                }
            }
            else{
                 cblas_dgemv(CblasColMajor, CblasNoTrans, mat_dim , mat_dim , 1.0 , M , mat_dim , u[1] , 1 , 1.0 , u[0] , 1);
                 alpha[ls] = cblas_ddot(mat_dim,u[1],1,u[0],1);
                cblas_daxpy(mat_dim, -alpha[ls], u[1],1,u[0],1);
                if(ls != tri_mat_dim -1){
                    beta[ls] = cblas_dnrm2(mat_dim, u[0],1);
                    cblas_dscal(mat_dim, 1./beta[ls], u[0],1);
                }
            }
            /*========================================================*/

            /*==================三重対角行列の数値対角化=================*/
            vec_init(tri_mat_dim, diag);
            vec_init(tri_mat_dim-1,sub_diag);

            if(ls % 2 == 0){
                cblas_dcopy(tri_mat_dim, alpha, 1 ,diag, 1);
                cblas_dcopy(tri_mat_dim - 1 , beta , 1 , sub_diag , 1);
               
                if(ls < tri_mat_dim - 1){
                    sub_diag[ls] = 0.;
                    int info = LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls+1, diag, sub_diag,tri_diag_evec,ls+1);
                    if (info != 0)
                    {
                        std::cout << "@ls = " << ls << " , LAPACKE_detev errored." << std::endl;
                    }
                }
                else{
                    int info = LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls+1, diag, sub_diag,tri_diag_evec,ls+1);
                    if (info != 0)
                    {
                        std::cout << "@ls = " << ls << " , LAPACKE_detev errored." << std::endl;
                    }
                }
                cblas_dcopy(tri_mat_dim, diag, 1, eval_even,1);
                if (info_ls == 'y')
                {
                    std::cout << "@ls = " << ls << " : eigen value = " << eval_even[0] << std::endl;
                }
                else if (info_ls == 's')
                {
                    std::cout << "@ls = " << ls << std::endl;
                }
            }
            else{
                
                cblas_dcopy(tri_mat_dim,alpha,1,diag,1);
                cblas_dcopy(tri_mat_dim-1,beta,1,sub_diag,1);
                
                if(ls < tri_mat_dim -1){
                    sub_diag[ls] = 0.;
                    int info = LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls+1, diag, sub_diag, tri_diag_evec, ls+1);
                    if (info != 0)
                    {
                        std::cout << "@ls = " << ls << " , LAPACKE_detev's error." << std::endl;
                    }
                }
                else{
                    int info = LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls+1, diag, sub_diag, tri_diag_evec,ls+1);
                    if (info != 0)
                    {
                        std::cout << "@ls = " << ls << " , LAPACKE_detev's error." << std::endl;
                    }
                }
                cblas_dcopy(tri_mat_dim, diag, 1, eval_odd, 1);
                if (info_ls == 'y')
                {
                    std::cout << "@ls = " << ls << " : eigen value = " << eval_odd[0] << std::endl;
                }
                else if (info_ls == 's')
                {
                    std::cout << "@ls = " << ls << std::endl;
                }
            }
            /*======================================================================*/

            /*============================収束状況の確認==============================*/
            if(ls > 0){
                eps = abs(eval_even[0] - eval_odd[0]);
                if(info_ls == 'y'){
                    std::cout << "eps = " << std::setprecision(16) << eps << std::endl;
                } 

                if(eps > err){
                    err_checker = true;
                }
                else{
                    err_checker = false;
                }
            }
            /*=====================================================================*/
        }
        else{
            std::cout << "break @" << count << std::endl;
            break;
        }
    }

    /*=========================基底状態の固有値============================*/
    double gstate_eval = 0.;
    if(count%2 == 0) gstate_eval = eval_even[0];
    else gstate_eval = eval_odd[0];

    std::cout << "eigen value of ground state = " << gstate_eval << std::endl;
    /*==========================配列リソースのリリース====================*/
    for(int i=0;i<2;i++){
        delete[] u[i];
    }
    delete[] u;
    delete[] alpha;
    delete[] beta;
    delete[] eval_even;
    delete[] eval_odd;
    delete[] diag;
    delete[] sub_diag;
    delete[] tri_diag_evec;

    return gstate_eval;
}

/*-----------------------------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------------------------------*/


double dns_lancozs_V(char info_ls, const int mat_dim, int tri_mat_dim, double *M, std::string outfile_name, char c, std::string testfile_name)
{
    int count = 0;
    double eps = 1.0;         //誤差を代入するようの変数
    double err = 1.0e-16;     //要求精度
    bool err_checker = true;  //誤差が要求精度の範囲内に収まっているかを確認するためのflag

    /*==========================初期ベクトルの用意=========================*/
    double **u = new double *[2];
    for(int k=0;k<2;k++){
        u[k] = new double[mat_dim];
    }

    std::random_device rand;
    std::mt19937 mt(rand());
    std::uniform_real_distribution<> rand1(0,1);
    for(int k=0;k<mat_dim;k++){
        u[0][k] = rand1(mt);
        u[1][k] = 0.0;
    }
    sdz(mat_dim,u[0]);

    double *evec = new double[mat_dim];
    cblas_dcopy(mat_dim,u[0],1,evec,1);
    /*=================================================================*/
    double *alpha = new double[tri_mat_dim];
    vec_init(tri_mat_dim,alpha);
    
    double *beta = new double[tri_mat_dim - 1];
    vec_init(tri_mat_dim - 1,beta);
    
    double *eval_even = new double[tri_mat_dim];
    vec_init(tri_mat_dim, eval_even);

    double *eval_odd = new double[tri_mat_dim];
    vec_init(tri_mat_dim, eval_odd);

    double *diag = new double[tri_mat_dim];
    vec_init(tri_mat_dim, diag);

    double *sub_diag = new double[tri_mat_dim - 1];
    vec_init(tri_mat_dim - 1, sub_diag);

    double *tri_diag_evec = new double[tri_mat_dim*tri_mat_dim];//lapack_dstevに渡す用('N'の場合は参照されないので適当に用意しておく)
    vec_init(tri_mat_dim*tri_mat_dim, tri_diag_evec);



    /*==========================lanczos法===============================*/
    for(int ls = 0;ls < tri_mat_dim; ls++){
        if(err_checker){
            count = ls;
            /*========メモリ削減のためのlanczosベクトルの更新===========*/
            if(ls > 0){
                if(ls % 2 == 0){
                    cblas_dscal(mat_dim, -beta[ls-1],u[1],1);
                }
                else{
                    cblas_dscal(mat_dim, -beta[ls-1],u[0],1);
                }
            }
            /*====================================================*/


            /*=============alpha, betaの計算とlancozsベクトルの更新==========*/
            if(ls % 2 == 0){
                cblas_dgemv(CblasColMajor, CblasNoTrans, mat_dim , mat_dim , 1.0 , M , mat_dim, u[0], 1 , 1.0 , u[1] , 1);
                alpha[ls] = cblas_ddot(mat_dim,u[1],1,u[0],1);
                cblas_daxpy(mat_dim,-alpha[ls], u[0], 1, u[1],1);
                if(ls != tri_mat_dim -1){
                    beta[ls] = cblas_dnrm2(mat_dim, u[1], 1);
                    cblas_dscal(mat_dim, 1./beta[ls],u[1],1);
                }
            }
            else{
                 cblas_dgemv(CblasColMajor, CblasNoTrans, mat_dim , mat_dim , 1.0 , M , mat_dim , u[1] , 1 , 1.0 , u[0] , 1);
                 alpha[ls] = cblas_ddot(mat_dim,u[1],1,u[0],1);
                cblas_daxpy(mat_dim, -alpha[ls], u[1],1,u[0],1);
                if(ls != tri_mat_dim -1){
                    beta[ls] = cblas_dnrm2(mat_dim, u[0],1);
                    cblas_dscal(mat_dim, 1./beta[ls], u[0],1);
                }
            }
            /*========================================================*/

            /*==================三重対角行列の数値対角化=================*/
            vec_init(tri_mat_dim, diag);
            vec_init(tri_mat_dim-1,sub_diag);

            if(ls % 2 == 0){
                cblas_dcopy(tri_mat_dim, alpha, 1 ,diag, 1);
                cblas_dcopy(tri_mat_dim - 1 , beta , 1 , sub_diag , 1);
               
                if(ls < tri_mat_dim - 1){
                    sub_diag[ls] = 0.;
                    int info = LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', ls+1, diag, sub_diag,tri_diag_evec,ls+1);
                    if (info != 0)
                    {
                        std::cout << "@ls = " << ls << " , LAPACKE_detev errored." << std::endl;
                    }
                }
                else{
                    int info = LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', ls+1, diag, sub_diag,tri_diag_evec,ls+1);
                    if (info != 0)
                    {
                        std::cout << "@ls = " << ls << " , LAPACKE_detev errored." << std::endl;
                    }
                }
                cblas_dcopy(tri_mat_dim, diag, 1, eval_even,1);
                if (info_ls == 'y')
                {
                    std::cout << "@ls = " << ls << " : eigen value = " << eval_even[0] << std::endl;
                }
                else if (info_ls == 's')
                {
                    std::cout << "@ls = " << ls << std::endl;
                }
            }
            else{
                
                cblas_dcopy(tri_mat_dim,alpha,1,diag,1);
                cblas_dcopy(tri_mat_dim-1,beta,1,sub_diag,1);
                
                if(ls < tri_mat_dim -1){
                    sub_diag[ls] = 0.;
                    int info = LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', ls+1, diag, sub_diag, tri_diag_evec, ls+1);
                    if (info != 0)
                    {
                        std::cout << "@ls = " << ls << " , LAPACKE_detev's error." << std::endl;
                    }
                }
                else{
                    int info = LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', ls+1, diag, sub_diag, tri_diag_evec,ls+1);
                    if (info != 0)
                    {
                        std::cout << "@ls = " << ls << " , LAPACKE_detev's error." << std::endl;
                    }
                }
                cblas_dcopy(tri_mat_dim, diag, 1, eval_odd, 1);
                if (info_ls == 'y')
                {
                    std::cout << "@ls = " << ls << " : eigen value = " << eval_odd[0] << std::endl;
                }
                else if (info_ls == 's')
                {
                    std::cout << "@ls = " << ls << std::endl;
                }
            }
            /*======================================================================*/

            /*============================収束状況の確認==============================*/
            if(ls > 0){
                eps = abs(eval_even[0] - eval_odd[0]);
                if(info_ls == 'y'){
                    std::cout << "eps = " << std::setprecision(16) << eps << std::endl;
                } 

                if(eps > err){
                    err_checker = true;
                }
                else{
                    err_checker = false;
                }
            }
            /*=====================================================================*/
        }
        else{
            std::cout << "break @" << count << std::endl;
            break;
        }
    }

    /*=========================基底状態の固有値============================*/
    double gstate_eval = 0.;
    if(count%2 == 0) gstate_eval = eval_even[0];
    else gstate_eval = eval_odd[0];
    std::cout << "eigen value of ground state = " << gstate_eval << std::endl;


    /*==============================基底状態の固有ベクトル===================================*/
    std::cout << "=======================================================================\n";
    std::cout << "start calculation of eigen vector\n";
    std::cout << "=======================================================================\n";
    vec_init(mat_dim, u[0]);
    vec_init(mat_dim, u[1]);
    cblas_dcopy(mat_dim, evec, 1, u[0], 1);
    for(int ls = 0;ls < count + 2; ls++){
        if (ls % 2 == 0)
        {
            if (ls == 0)
                cblas_dscal(mat_dim, tri_diag_evec[ls], evec, 1);
            else
                cblas_daxpy(mat_dim, tri_diag_evec[ls], u[0], 1, evec, 1);
        }
        else
        {
            cblas_daxpy(mat_dim, tri_diag_evec[ls], u[1], 1, evec, 1);
        }

        if (ls > 0)
        {
            if (ls % 2 == 0)
            {
                cblas_dscal(mat_dim, -beta[ls - 1], u[1], 1);
            }
            else
            {
                cblas_dscal(mat_dim, -beta[ls - 1], u[0], 1);
            }
        }

        if(ls % 2 == 0){
            cblas_dgemv(CblasColMajor, CblasNoTrans, mat_dim , mat_dim , 1.0 , M , mat_dim, u[0], 1 , 1.0 , u[1] , 1);
            cblas_daxpy(mat_dim,-alpha[ls], u[0], 1, u[1],1);
            if(ls != tri_mat_dim -1){
                cblas_dscal(mat_dim, 1./beta[ls],u[1],1);
            }
        }
        else{
            cblas_dgemv(CblasColMajor, CblasNoTrans, mat_dim , mat_dim , 1.0 , M , mat_dim , u[1] , 1 , 1.0 , u[0] , 1);
            cblas_daxpy(mat_dim, -alpha[ls], u[1],1,u[0],1);
            if(ls != tri_mat_dim -1){
                cblas_dscal(mat_dim, 1./beta[ls], u[0],1);
            }
        }
    }
    sdz(mat_dim,evec);
    /*==============================テストモードでの処理=====================================*/
    if(c == 'T'){
        std::ofstream ofs_testfile(testfile_name,std::ios::app);
        ofs_testfile << "eigen vector = [" << ",";
        for(int i=0;i<mat_dim;i++){
            ofs_testfile << evec[i] << ",";
        }
        ofs_testfile << "]\n";
    }
    /*==========================配列リソースのリリース====================*/
    for(int i=0;i<2;i++){
        delete[] u[i];
    }
    delete[] u;
    delete[] alpha;
    delete[] beta;
    delete[] eval_even;
    delete[] eval_odd;
    delete[] diag;
    delete[] sub_diag;
    delete[] tri_diag_evec;
    delete[] evec;
    return gstate_eval;
}
