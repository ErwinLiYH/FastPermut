from fastpermut import permutation
import fastpermut


chro, loca, annot = fastpermut.load_from_file("../test_data")
csv = permutation.permutation_mkl(chro, loca, annot, ["first","second","third"], 10, to_sparse=False, omp=False)

csv.to_csv("result.csv", index=False)