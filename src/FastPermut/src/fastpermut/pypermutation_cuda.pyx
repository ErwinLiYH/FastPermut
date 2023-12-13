# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION

import numpy as np
import pandas as pd
from . cimport pypermutation_cuda
from . import timing_decorator

from libc.stdlib cimport malloc, free

def init_cuda():
    status = pypermutation_cuda.init_cuda()
    if status==0:
        print("CUDA initialized", flush=True)

@timing_decorator
def permutation_cuda(chromosome_array: np.ndarray, location_array: np.ndarray, annotation_array: np.ndarray, columns: list | None, permutation_times: int, to_sparse=True, seed=None, record_random_number=False) -> pd.DataFrame:
    _, N_SNPs_py = np.unique(chromosome_array, return_counts=True)
    cdef int N_times = permutation_times
    cdef int N_chros = len(N_SNPs_py)
    cdef int N_annots = annotation_array.shape[1]
    cdef int N_SNP = annotation_array.shape[0]
    cdef int[:]  N_SNPs = np.asarray(N_SNPs_py, dtype=np.intc)

    cdef float[:] location_vecs = location_array
    cdef float[:,:] annot_mats = annotation_array
    cdef float[:,:] count_mats = np.zeros((N_chros*N_times, N_annots), dtype=np.single)

    cdef int* RN_Column
    if record_random_number:
        RN_Column = <int*>malloc(N_chros * N_times * sizeof(int))
    else:
        RN_Column = NULL

    cdef unsigned long long seed_c
    cdef unsigned long long* seed_c_p = &seed_c
    if seed==None:
        seed_c_p = NULL
    else:
        seed_c = seed

    if to_sparse:
        pypermutation_cuda.cuda_permut_sparse(seed_c_p, &location_vecs[0], &N_SNPs[0], N_chros, &annot_mats[0, 0], N_annots, N_times, &count_mats[0, 0], RN_Column)
    else:
        pypermutation_cuda.cuda_permut(seed_c_p, &location_vecs[0], &N_SNPs[0], N_chros, &annot_mats[0, 0], N_annots, N_times, &count_mats[0, 0], RN_Column)

    result = pd.DataFrame(np.asarray(count_mats), columns=columns)
    cdef int[:] mv
    if record_random_number:
        # create memoryview from pointer
        mv = <int[:N_chros * N_times]>RN_Column
        # create numpy array from memoryview
        arr = np.asarray(mv, dtype=np.intc)
        # create a copy
        arr = arr.copy()

        # release memory of RN_Column
        free(RN_Column)

        result["RN"] = arr
    return result