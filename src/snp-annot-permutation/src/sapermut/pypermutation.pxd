cdef extern from "permutation.hpp":

    void cmkl_permut(unsigned int* seed_p, float* location_vecs, int* N_SNPs, int N_chros, float* annot_mats, int N_annots, int N_times, float* count_mats)

    void cmkl_permut_omp(unsigned int* seed_p, float* location_vecs, int* N_SNPs, int N_chros, float* annot_mats, int N_annots, int N_times, float* count_mats)

    void cmkl_permut_sparse(unsigned int* seed_p, float* location_vecs, int* N_SNPs, int N_chros, float* annot_mats, int N_annots, int N_times, float* count_mats)

    void cmkl_permut_sparse_omp(unsigned int* seed_p, float* location_vecs, int* N_SNPs, int N_chros, float* annot_mats, int N_annots, int N_times, float* count_mats)