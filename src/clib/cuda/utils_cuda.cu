#include "permutation.cuh"

__global__ void rotate_to_Kernel(float* vector, float* rotated_vector, int N, float* K){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N) {
        // Perform the rotation operation
        int newIdx = (idx + int(floor((*K)*N))) % N;
        // printf("thread %d, new id %d, K %f, N %d, original value %f\n",idx,newIdx,*K,N,*(vector+idx));
        // Write to a temporary variable to avoid overwriting the original vector
        *(rotated_vector+newIdx) = *(vector+idx);
    }
}

int init_cuda(){
    cudaError_t status = cudaFree(0);
    return status;
}

void gen_random_cuda(float* r, int N, unsigned long long seed){
    curandGenerator_t gen;
    curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_MT19937);
    curandSetPseudoRandomGeneratorSeed(gen, seed);
    curandGenerateUniform(gen, r, N);
    #if TEST_PRINT_ARRAY==1
    float* h_r = new float[N];
    cudaMemcpy(h_r, r, N*sizeof(float), cudaMemcpyDeviceToHost);
    print_1D_F_array(h_r, N);
    delete [] h_r;
    #endif
}

void gen_sparse_mat(float* annot_mats, int NCols, int NRows, int lda, float** values, int** indexes, int** pointers, int* size_nnz, int* size_p){
    // calculate number of non-zero values(nnz)----------------------
    int nnz=0;
    #if ANNOT_MAT_ROW_OR_COLUMN_BASED == 1
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
    #if ANNOT_MAT_ROW_OR_COLUMN_BASED == 1
    *pointers = new int[NRows+1];
    *size_p = NRows+1;
    #else
    *pointers = new int[NCols+1];
    *size_p = NCols+1;
    #endif
    *values = new float[nnz];
    *indexes = new int[nnz];
    *size_nnz = nnz;
    //---------------------------------------------------------------

    // fill outputs---------------------------------------------------
    nnz=0;// reuse nnz as non-zero value id
    #if ANNOT_MAT_ROW_OR_COLUMN_BASED == 1
    for (int i=0; i<NRows; ++i){
        (*pointers)[i] = nnz;
        for (int j=0; j<NCols; ++j){
            if (*(annot_mats+i*lda+j)!=0){
                (*values)[nnz] = *(annot_mats+i*lda+j);
                (*indexes)[nnz] = j;
                nnz+=1;
            }
        }
    }
    (*pointers)[NRows] = nnz;
    #else
    for (int j=0; j<NCols; ++j){
        (*pointers)[j] = nnz;
        for (int i=0; i<NRows; ++i){
            if (*(annot_mats+j*lda+i)!=0){
                (*values)[nnz] = *(annot_mats+j*lda+i);
                (*indexes)[nnz] = i;
                nnz+=1;
            }
        }
    }
    (*pointers)[NCols] = nnz;
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
    std::cout << "pointers:" << std::endl;
    #if ANNOT_MAT_ROW_OR_COLUMN_BASED == 1
    for (int i=0; i<NRows+1; ++i){
    #else
    for (int i=0; i<NCols+1; ++i){
    #endif
        std::cout << (*pointers)[i] << ", ";
    }
    std::cout<<std::endl;
    #endif
    //---------------------------------------------------------------
}