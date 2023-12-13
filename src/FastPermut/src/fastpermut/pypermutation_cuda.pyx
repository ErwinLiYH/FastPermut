# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION

import numpy as np
import pandas as pd
from . cimport pypermutation_cuda
from . import timing_decorator


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

    if record_random_number:
        cdef int[:] RN_Column = np.zeros(N_chros*N_times, dtype=np.intc)
    else:
        cdef int* RN_Column = NULL

    cdef unsigned long long seed_c
    cdef unsigned long long* seed_c_p = &seed_c
    if seed==None:
        seed_c_p = NULL
    else:
        seed_c = seed

    if to_sparse:
        pypermutation_cuda.cuda_permut_sparse(seed_c_p, &location_vecs[0], &N_SNPs[0], N_chros, &annot_mats[0, 0], N_annots, N_times, &count_mats[0, 0], &RN_Column[0])
    else:
        pypermutation_cuda.cuda_permut(seed_c_p, &location_vecs[0], &N_SNPs[0], N_chros, &annot_mats[0, 0], N_annots, N_times, &count_mats[0, 0], &RN_Column[0])

    result = pd.DataFrame(np.asarray(count_mats), columns=columns)
    if record_random_number:
        result["RN"] = np.asarray(RN_Column)
    return result