# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION

import numpy as np
import pandas as pd
from . cimport pypermutation
import time


def timing_decorator(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()  # record the start time
        result = func(*args, **kwargs)  # call the original function
        end_time = time.time()  # record the end time
        print(f"Function {func.__name__} took {end_time - start_time} seconds to run.", flush=True)
        return result
    return wrapper

@timing_decorator
def permutation_mkl(chromosome_array: np.ndarray, location_array: np.ndarray, annotation_array: np.ndarray, columns: list | None, permutation_times: int, to_sparse=True, omp=True, seed=None) -> pd.DataFrame:
    _, N_SNPs_py = np.unique(chromosome_array, return_counts=True)
    cdef int N_times = permutation_times
    cdef int N_chros = len(N_SNPs_py)
    cdef int N_annots = annotation_array.shape[1]
    cdef int N_SNP = annotation_array.shape[0]
    cdef int[:] N_SNPs = np.asarray(N_SNPs_py, dtype=np.intc)

    cdef float[:] location_vecs = location_array
    cdef float[:,:] annot_mats = annotation_array
    cdef float[:,:] count_mats = np.zeros((N_chros*N_times, N_annots), dtype=np.single)

    cdef unsigned int seed_c
    cdef unsigned int* seed_c_p = &seed_c
    if seed==None:
        seed_c_p = NULL
    else:
        seed_c = seed

    if to_sparse:
        if omp:
            pypermutation.cmkl_permut_sparse_omp(seed_c_p, &location_vecs[0], &N_SNPs[0], N_chros, &annot_mats[0,0], N_annots, N_times, &count_mats[0,0])
        else:
            pypermutation.cmkl_permut_sparse(seed_c_p, &location_vecs[0], &N_SNPs[0], N_chros, &annot_mats[0,0], N_annots, N_times, &count_mats[0,0])
    else:
        if omp:
            pypermutation.cmkl_permut_omp(seed_c_p, &location_vecs[0], &N_SNPs[0], N_chros, &annot_mats[0,0], N_annots, N_times, &count_mats[0,0])
        else:
            pypermutation.cmkl_permut(seed_c_p, &location_vecs[0], &N_SNPs[0], N_chros, &annot_mats[0,0], N_annots, N_times, &count_mats[0,0])

    return pd.DataFrame(np.asarray(count_mats), columns=columns)

def permutation_mkl_mpi(chromosome_array: np.ndarray, location_array: np.ndarray, annotation_array: np.ndarray, columns: list | None, permutation_times: int, to_sparse=True, name="test", seed=None):
    import mpi4py
    mpi4py.rc.initialize = False
    mpi4py.rc.finalize = False

    from mpi4py import MPI
    start_time = time.time()
    required = MPI.THREAD_FUNNELED
    provided = MPI.Init_thread(required)
    if provided < required:
        print('Warning: The level of thread support is less than the required one')

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    remains = permutation_times%size

    py_N_times = int(permutation_times/size)

    if rank<remains:
        py_N_times+=1

    print("Total size: %d, rank: %d, permutation times: %d"%(size, rank, py_N_times) ,flush=True)

    _, N_SNPs_py = np.unique(chromosome_array, return_counts=True)
    cdef int N_times = py_N_times
    cdef int N_chros = len(N_SNPs_py)
    cdef int N_annots = annotation_array.shape[1]
    cdef int N_SNP = annotation_array.shape[0]
    cdef int[:] N_SNPs = np.asarray(N_SNPs_py, dtype=np.intc)

    cdef float[:] location_vecs = location_array
    cdef float[:,:] annot_mats = annotation_array
    cdef float[:,:] count_mats = np.zeros((N_chros*N_times, N_annots), dtype=np.single)

    cdef unsigned int seed_c
    cdef unsigned int* seed_c_p = &seed_c
    if seed==None:
        seed_c_p = NULL
    else:
        seed_c = seed

    if to_sparse:
        pypermutation.cmkl_permut_sparse_omp(seed_c_p, &location_vecs[0], &N_SNPs[0], N_chros, &annot_mats[0,0], N_annots, N_times, &count_mats[0,0])
    else:
        pypermutation.cmkl_permut_omp(seed_c_p, &location_vecs[0], &N_SNPs[0], N_chros, &annot_mats[0,0], N_annots, N_times, &count_mats[0,0])

    end_time = time.time()
    print(f"Function permutation_mkl_mpi rank {rank} took {end_time - start_time} seconds to run.", flush=True)

    res = pd.DataFrame(np.asarray(count_mats), columns=columns)
    res.to_csv("%s_R%d.csv"%(name, rank), index=False)

    MPI.Finalize()