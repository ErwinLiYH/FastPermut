#include <mkl.h>
#include <algorithm>
#include <random>
#include <chrono>
#include <sys/time.h>
#include <iomanip>
#include <omp.h>
#include <mpi.h>
#include <vector>
#include <cstdlib>
#include <sstream>
#include "../common/utils.hpp"

template <class T> using Vec = std::vector<T>;

void gen_sparse_mat(float* annot_mats, MKL_INT NCols, MKL_INT NRows, int lda,float** values, MKL_INT** indexes, MKL_INT** pointerB, MKL_INT** pointerE);

void gen_random_mkl(int* r, int N, int min, int max, unsigned int seed);

void cmkl_permut(unsigned int* seed_p, float* location_vecs, int* N_SNPs, int N_chros, float* annot_mats, int N_annots, int N_times, float* count_mats);

void cmkl_permut_omp(unsigned int* seed_p, float* location_vecs, int* N_SNPs, int N_chros, float* annot_mats, int N_annots, int N_times, float* count_mats);

void cmkl_permut_sparse(unsigned int* seed_p, float* location_vecs, int* N_SNPs, int N_chros, float* annot_mats, int N_annots, int N_times, float* count_mats);

void cmkl_permut_sparse_omp(unsigned int* seed_p, float* location_vecs, int* N_SNPs, int N_chros, float* annot_mats, int N_annots, int N_times, float* count_mats);

template <typename InputIterator, typename OutputIterator>
void rotate_to(InputIterator first, InputIterator middle, InputIterator last, OutputIterator destination) {
    // Copy from middle to last
    OutputIterator new_dest = std::copy(middle, last, destination);
    // Copy from first to middle
    std::copy(first, middle, new_dest);
}