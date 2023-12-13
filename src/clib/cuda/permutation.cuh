#include <cmath>
#include <chrono>
#include <iomanip>
#include <cmath>
#include <sys/time.h>
#include <cublas_v2.h>
#include <curand.h>
#include <cusparse_v2.h>
#include "../common/utils.hpp"

__global__ void rotate_to_Kernel(float* vector, float* rotated_vector, int N, float* K);

void gen_random_cuda(float* r, int N, unsigned long long seed);

int init_cuda();

void cuda_permut(unsigned long long* seed_p, float* location_vecs, int* N_SNPs, int N_chros, float* annot_mats, int N_annots, int N_times, float* count_mats, int* RNsCol);

void gen_sparse_mat(float* annot_mats, int NCols, int NRows, int lda, float** values, int** indexes, int** pointers, int* size_nnz, int* size_p);

void cuda_permut_sparse(unsigned long long* seed_p, float* location_vecs, int* N_SNPs, int N_chros, float* annot_mats, int N_annots, int N_times, float* count_mats, int* RNsCol);

#if IF_TEST==1
#define CHECK_CUDA(func){                                                      \
    cudaError_t status = (func);                                               \
    if (status != cudaSuccess) {                                               \
        printf("CUDA API failed at line %d with error: %s (%d)\n",             \
               __LINE__, cudaGetErrorString(status), status);                  \
    }                                                                          \
}

#define CHECK_CUSPARSE(func){                                                  \
    cusparseStatus_t status = (func);                                          \
    if (status != CUSPARSE_STATUS_SUCCESS) {                                   \
        printf("CUSPARSE API failed at line %d with error: %s (%d)\n",         \
               __LINE__, cusparseGetErrorString(status), status);              \
    }                                                                          \
}

#define CHECK_CUBLAS(func){                                                    \
    cublasStatus_t status = (func);                                            \
    if (status != CUBLAS_STATUS_SUCCESS) {                                     \
        printf("CUBLAS API failed at line %d with error code: %d\n",           \
                __LINE__, status);                                             \
    }                                                                          \
}
#else
#define CHECK_CUDA(func) (func);
#define CHECK_CUSPARSE(func) (func);
#define CHECK_CUBLAS(func) (func);
#endif