#include "permutation.hpp"

int main(){
    float location_vecs[6]={1,2,10,11,12};
    int N_SNPs[2] = {2,3};
    int N_chros = 2;
    float annot_mats[10] = {123,234,345,456,567,1234,2345,3456,4567,5678};
    int N_annots = 2;
    int N_times = 10;
    float count_mats[20];

    cmkl_permut(location_vecs, N_SNPs, N_chros, annot_mats, N_annots, N_times, count_mats);

    std::cout<<count_mats[0]<<", "<<count_mats[1]<<std::endl;
}