from fastpermut import permutation
import fastpermut


chro, loca, annot = fastpermut.load_from_file("../test_data")
permutation.permutation_mkl_mpi(chro, loca, annot, ["first","second","third"], 10, to_sparse=False)