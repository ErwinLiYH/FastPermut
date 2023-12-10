cdef extern from "permutation.cuh":

    void cuda_permut(unsigned long long* seed_p, float* location_vecs, int* N_SNPs, int N_chros, float* annot_mats, int N_annots, int N_times, float* count_mats)

    void cuda_permut_sparse(unsigned long long* seed_p, float* location_vecs, int* N_SNPs, int N_chros, float* annot_mats, int N_annots, int N_times, float* count_mats)

    int init_cuda()