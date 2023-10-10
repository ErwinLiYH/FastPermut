#include "permutation.hpp"

void cmkl_permut(unsigned int* seed_p, float* location_vecs, int* N_SNPs, int N_chros, float* annot_mats, int N_annots, int N_times, float* count_mats){
    //>>> timer------------------
    #if IF_TIMING==1
    struct timeval start0;
    gettimeofday(&start0, NULL);
    #endif
    //<<<------------------------

    // Generate the SNP range for each chromosome---------------------------------
    int SNP_ranges[N_chros+1];
    gen_snp_ranges(&(SNP_ranges[0]), N_SNPs, N_chros);
    //----------------------------------------------------------------------------

    // Create random values matrix for all chromosomes and all permutation steps--
    int RNs[N_chros][N_times];
    for (int i=0; i<N_chros; ++i){
        if (seed_p == NULL){
            unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
            gen_random_mkl(&(RNs[i][0]), N_times, 0, N_SNPs[i], seed);
        }else{
            gen_random_mkl(&(RNs[i][0]), N_times, 0, N_SNPs[i], *seed_p*i);
        }
    }

    //>>> store the random values in IF_STORE_RNS mode----
    #if IF_STORE_RNS==1
    std::ofstream outRNs("cmkl_permut_RNs", std::ios::binary);
    for (int i=0; i<N_chros; ++i){
        for (int j=0; j<N_times; ++j){
            outRNs.write(reinterpret_cast<char*>(&(RNs[i][j])), sizeof(int));
        }
    }
    outRNs.close();
    #endif
    //<<<-------------------------------------------------

    //----------------------------------------------------------------------------

    //>>> timer------------------
    #if IF_TIMING==1
    struct timeval start, end;
    gettimeofday(&start, NULL);
    #endif
    //<<<------------------------

    // do permutation-------------------------------------------------------------
    for (int i=0; i< N_chros; ++i){// chromosomes loop
        for (int j=0; j<N_times; ++j){ // permutations loop

            //>>> print process in TEST mode---------
            #if IF_TEST == 1
            if (j%PRINT_INTERVAL==0){
                std::cout << "process chromosome " << i << ", " << "permutation " << j << std::endl;
            }
            #endif
            //<<< -----------------------------------

            //>>> print origianl location_vecs and radom number for rotate in TEST_PRINT_ARRAY mode
            #if TEST_PRINT_ARRAY == 1
            std::cout << "random number: " << RNs[i][j] << std::endl;
            std::cout << "Location vectors before rotate:" << std::endl;
            print_location_vecs(location_vecs, N_chros, N_SNPs);
            #endif
            //<<< --------------------------------------------------------------------------

            // rotate location vector
            std::rotate(location_vecs+SNP_ranges[i], location_vecs+SNP_ranges[i]+RNs[i][j], location_vecs+SNP_ranges[i+1]);

            //>>> print rotaed location_vecs in TEST_PRINT_ARRAY mode-----
            #if TEST_PRINT_ARRAY == 1
            std::cout << "Location vectors after rotate:" << std::endl;
            print_location_vecs(location_vecs, N_chros, N_SNPs);
            #endif
            //<<<----------------------------------------------------------

            // calculate the matrix times vector by mkl ------------------------------------------
            #if ANNOT_MAT_ROW_OR_COLUMN_BASED == 1
            // when annotation matrix store in row-base layout
            cblas_sgemv(CblasRowMajor, CblasTrans, N_SNPs[i], N_annots, 1.0, annot_mats+SNP_ranges[i]*N_annots, N_annots, location_vecs+SNP_ranges[i], 1, 0.0, count_mats+i*N_times*N_annots+j*N_annots, 1);
            #else
            // when annotation matrix store in column-base layout
            cblas_sgemv(CblasColMajor, CblasTrans, N_SNPs[i], N_annots, 1.0, annot_mats+SNP_ranges[i], SNP_ranges[N_chros], location_vecs+SNP_ranges[i], 1, 0.0, count_mats+i*N_times*N_annots+j*N_annots, 1);
            #endif
            //------------------------------------------------------------------------------------

            //>>> print count matrixes in TEST_PRINT_ARRAY mode-----
            #if TEST_PRINT_ARRAY == 1
            std::cout << "count matrixes:" << i << std::endl;
            print_count_mats(count_mats, N_chros, N_times, N_annots);
            #endif
            //<<<---------------------------------------------------
        }
    }
    //----------------------------------------------------------------------------

    //>>> timer------------------
    #if IF_TIMING==1
    int n_threads = mkl_get_max_threads();
    gettimeofday(&end, NULL);
    double time_taken, time_taken0;
    time_taken = (end.tv_sec - start.tv_sec) * 1e6;
    time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;
    time_taken0 = (end.tv_sec - start0.tv_sec) * 1e6;
    time_taken0 = (time_taken0 + (end.tv_usec - start0.tv_usec)) * 1e-6;
    if (ANNOT_MAT_ROW_OR_COLUMN_BASED==0){
        std::cout << "Column-major variant" << std::endl;
    }else{
        std::cout << "Row-major variant" << std::endl;
    }
    std::cout << "MKL number of threads: " << n_threads << std::endl;
    std::cout << "cmkl_permut total time: " << std::fixed << std::setprecision(10) << time_taken0 << std::endl;
    std::cout << "cmkl_permut permutation time: " << std::fixed << std::setprecision(10) << time_taken << std::endl;
    std::cout.unsetf(std::ios_base::floatfield);
    std::cout << std::setprecision(6);
    #endif
    //<<<------------------------
}

void cmkl_permut_sparse(unsigned int* seed_p, float* location_vecs, int* N_SNPs, int N_chros, float* annot_mats, int N_annots, int N_times, float* count_mats){
    //>>> timer------------------
    #if IF_TIMING==1
    struct timeval start0;
    gettimeofday(&start0, NULL);
    #endif
    //<<<------------------------

    // used in mkl_sparse_s_mv to describe matrix, only SPARSE_MATRIX_TYPE_GENERAL needed
    // SPARSE_FILL_MODE_LOWER and SPARSE_DIAG_NON_UNIT will be ignored
    struct matrix_descr descr = {SPARSE_MATRIX_TYPE_GENERAL, SPARSE_FILL_MODE_LOWER, SPARSE_DIAG_NON_UNIT};

    // Generate the SNP range for each chromosome---------------------------------
    int SNP_ranges[N_chros+1];
    gen_snp_ranges(&(SNP_ranges[0]), N_SNPs, N_chros);
    //----------------------------------------------------------------------------

    // Create random values matrix for all chromosomes and all permutation steps--
    int RNs[N_chros][N_times];
    for (int i=0; i<N_chros; ++i){
        if (seed_p == NULL){
            unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
            gen_random_mkl(&(RNs[i][0]), N_times, 0, N_SNPs[i], seed);
        }else{
            gen_random_mkl(&(RNs[i][0]), N_times, 0, N_SNPs[i], *seed_p*i);
        }
    }

    //>>> store the random values in IF_STORE_RNS mode----
    #if IF_STORE_RNS==1
    std::ofstream outRNs("cmkl_permut_sparse_RNs", std::ios::binary);
    for (int i=0; i<N_chros; ++i){
        for (int j=0; j<N_times; ++j){
            outRNs.write(reinterpret_cast<char*>(&(RNs[i][j])), sizeof(int));
        }
    }
    outRNs.close();
    #endif
    //<<<-------------------------------------------------

    //----------------------------------------------------------------------------

    //Set up sparse matrixes---------------------------------------------------------
    float* values[N_chros];
    MKL_INT* indexes[N_chros];
    MKL_INT* pointerB[N_chros];
    MKL_INT* pointerE[N_chros];
    sparse_matrix_t smats[N_chros];
    for (int i=0; i<N_chros; ++i){
        #if ANNOT_MAT_ROW_OR_COLUMN_BASED == 1
        gen_sparse_mat(annot_mats+SNP_ranges[i]*N_annots, N_annots, N_SNPs[i], N_annots, &(values[i]), &(indexes[i]), &(pointerB[i]), &(pointerE[i]));
        mkl_sparse_s_create_csr(&(smats[i]), SPARSE_INDEX_BASE_ZERO, N_SNPs[i], N_annots, pointerB[i], pointerE[i], indexes[i], values[i]);
        #else
        gen_sparse_mat(annot_mats+SNP_ranges[i], N_annots, N_SNPs[i], SNP_ranges[N_chros], &(values[i]), &(indexes[i]), &(pointerB[i]), &(pointerE[i]));
        mkl_sparse_s_create_csc(&(smats[i]), SPARSE_INDEX_BASE_ZERO, N_SNPs[i], N_annots, pointerB[i], pointerE[i], indexes[i], values[i]);
        #endif
    }
    //----------------------------------------------------------------------------

    //>>> timer------------------
    #if IF_TIMING==1
    struct timeval start, end;
    gettimeofday(&start, NULL);
    #endif
    //<<<------------------------

    //begin permutation-----------------------------------------------------------
    for (int i=0; i<N_chros; ++i){
        for (int j=0; j<N_times; ++j){ // permutations loop

            //>>> print process in TEST mode---------
            #if IF_TEST == 1
            if (j%PRINT_INTERVAL==0){
                std::cout << "process chromosome " << i << ", " << "permutation " << j << std::endl;
            }
            #endif
            //<<<------------------------------------

            //>>> print origianl location_vecs and radom number for rotate in TEST_PRINT_ARRAY mode
            #if TEST_PRINT_ARRAY == 1
            std::cout << "random number:" << RNs[i][j] << std::endl;
            std::cout << "Location vectors before rotate:" << std::endl;
            print_location_vecs(location_vecs, N_chros, N_SNPs);
            #endif
            //<<<--------------------------------------------------------------------

            // rotate location vector
            std::rotate(location_vecs+SNP_ranges[i], location_vecs+SNP_ranges[i]+RNs[i][j], location_vecs+SNP_ranges[i+1]);

            //>>> print rotaed location_vecs in TEST_PRINT_ARRAY mode-----
            #if TEST_PRINT_ARRAY == 1
            std::cout << "Location vectors after rotate:" << std::endl;
            print_location_vecs(location_vecs, N_chros, N_SNPs);
            #endif
            //<<<---------------------------------------------------------

            // calculate the matrix times vector by mkl sparse
            int status = mkl_sparse_s_mv(SPARSE_OPERATION_TRANSPOSE, 1.0, smats[i], descr, location_vecs+SNP_ranges[i], 0.0 , count_mats+i*N_times*N_annots+j*N_annots);

            //>>> print count matrixes in TEST_PRINT_ARRAY mode-----
            #if TEST_PRINT_ARRAY == 1
            std::cout << "count matrixes:" << i << std::endl;
            print_count_mats(count_mats, N_chros, N_times, N_annots);
            #endif
            //<<<---------------------------------------------------
        }
    }
    //----------------------------------------------------------------------------

    //>>> timer------------------
    #if IF_TIMING==1
    int n_threads = mkl_get_max_threads();
    gettimeofday(&end, NULL);
    double time_taken, time_taken0;
    time_taken = (end.tv_sec - start.tv_sec) * 1e6;
    time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;
    time_taken0 = (end.tv_sec - start0.tv_sec) * 1e6;
    time_taken0 = (time_taken0 + (end.tv_usec - start0.tv_usec)) * 1e-6;
    if (ANNOT_MAT_ROW_OR_COLUMN_BASED==0){
        std::cout << "Column-major variant" << std::endl;
    }else{
        std::cout << "Row-major variant" << std::endl;
    }
    std::cout << "MKL number of threads: " << n_threads << std::endl;
    std::cout << "cmkl_permut_sparse total time: "<< std::fixed << std::setprecision(10) << time_taken0 << std::endl;
    std::cout << "cmkl_permut_sparse permutation time: "<< std::fixed << std::setprecision(10) << time_taken << std::endl;
    std::cout.unsetf(std::ios_base::floatfield);
    std::cout << std::setprecision(6);
    #endif
    //<<<------------------------
}