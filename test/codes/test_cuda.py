from fastpermut import permutation_cuda
import fastpermut

chro, loca, annot = fastpermut.load_from_file("../test_data")
csv = permutation_cuda.permutation_cuda(chro, loca, annot, ["first","second","third"], 10, False)

csv.to_csv("result.csv", index=False)