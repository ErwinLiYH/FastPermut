from sapermut import permutation
import sapermut


chro, loca, annot = sapermut.load_from_file("../test_data")
csv = permutation.permutation_mkl(chro, loca, annot, ["first","second","third"], 10, to_sparse=False, omp=True)

csv.to_csv("result.csv", index=False)