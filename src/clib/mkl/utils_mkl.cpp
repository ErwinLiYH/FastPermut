#include "permutation.hpp"

void gen_random_mkl(int* r, int N, int min, int max, unsigned int seed){
    int original_num_threads = mkl_get_max_threads();
    int num_cores = omp_get_num_procs();
    mkl_set_num_threads(num_cores);
    VSLStreamStatePtr stream;
    vslNewStream(&stream, VSL_BRNG_MT19937, seed);
    viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, N, r, min, max);
    mkl_set_num_threads(original_num_threads);
    #if TEST_PRINT_ARRAY==1
    std::cout << "Random numbers: " << std::endl;
    for (int i=0; i<N; ++i){
        std::cout<< *(r+i) << ", ";
    }
    std::cout << std::endl;
    #endif
}

void gen_sparse_mat(float* annot_mats, MKL_INT NCols, MKL_INT NRows, int lda,float** values, MKL_INT** indexes, MKL_INT** pointerB, MKL_INT** pointerE){
    // calculate number of non-zero values(nnz)----------------------
    int nnz=0;
    #if ANNOT_MAT_ROW_OR_COLUMN_BASED==1
    for (int i=0; i<NRows; ++i){
        for (int j=0; j<NCols; ++j){
            if (*(annot_mats+i*lda+j)!=0){
                nnz+=1;
            }
        }
    }
    #else
    for (int j=0; j<NCols; ++j){
        for (int i=0; i<NRows; ++i){
            if (*(annot_mats+j*lda+i)!=0){
                nnz+=1;
            }
        }
    }
    #endif
    //---------------------------------------------------------------

    // malloc memory for outputs-------------------------------------
    #if ANNOT_MAT_ROW_OR_COLUMN_BASED==1
    *pointerB = new MKL_INT[NRows];
    *pointerE = new MKL_INT[NRows];
    #else
    *pointerB = new MKL_INT[NCols];
    *pointerE = new MKL_INT[NCols];
    #endif
    *values = new float[nnz];
    *indexes = new MKL_INT[nnz];
    //---------------------------------------------------------------

    // fill outpus---------------------------------------------------
    nnz=0;// reuse nnz as non-zero value id
    #if ANNOT_MAT_ROW_OR_COLUMN_BASED==1
    for (int i=0; i<NRows; ++i){
        (*pointerB)[i] = nnz;
        for (int j=0; j<NCols; ++j){
            if (*(annot_mats+i*lda+j)!=0){
                (*values)[nnz] = *(annot_mats+i*lda+j);
                (*indexes)[nnz] = j;
                nnz+=1;
            }
        }
        (*pointerE)[i] = nnz;
    }
    #else
    for (int j=0; j<NCols; ++j){
        (*pointerB)[j] = nnz;
        for (int i=0; i<NRows; ++i){
            if (*(annot_mats+j*lda+i)!=0){
                (*values)[nnz] = *(annot_mats+j*lda+i);
                (*indexes)[nnz] = i;
                nnz+=1;
            }
        }
        (*pointerE)[j] = nnz;
    }
    #endif
    //---------------------------------------------------------------

    // test corectness by print--------------------------------------
    #if TEST_PRINT_ARRAY==1
    std::cout<< "nnz "<<nnz<<std::endl;
    std::cout << "Values:" << std::endl;
    for (int i=0; i<nnz; ++i){
        std::cout << (*values)[i] << ", ";
    }
    std::cout << std::endl;
    std::cout << "Indexes:" << std::endl;
    for (int i=0; i<nnz; ++i){
        std::cout << (*indexes)[i] << ", ";
    }
    std::cout << std::endl;
    std::cout << "pointerB:" << std::endl;
    #if ANNOT_MAT_ROW_OR_COLUMN_BASED==1
    for (int i=0; i<NRows; ++i){
    #else
    for (int i=0; i<NCols; ++i){
    #endif
        std::cout << (*pointerB)[i] << ", ";
    }
    std::cout<<std::endl;
    std::cout << "pointerE:" << std::endl;
    #if ANNOT_MAT_ROW_OR_COLUMN_BASED==1
    for (int i=0; i<NRows; ++i){
    #else
    for (int i=0; i<NCols; ++i){
    #endif
        std::cout << (*pointerE)[i] << ", ";
    }
    std::cout<<std::endl;
    #endif
    //---------------------------------------------------------------
}