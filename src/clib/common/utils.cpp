#include "utils.hpp"

void test_order(float* arr, int size){
    for (int i=0; i<size; ++i){
        std::cout << *(arr+i) << ", ";
    }
    std::cout<<std::endl;
}

void print_location_vecs(float* arr, int N_chros, int* N_SNPs){
    int SNP_ranges[N_chros+1];
    gen_snp_ranges(&(SNP_ranges[0]), N_SNPs, N_chros);
    for (int i=0; i<N_chros; ++i){
        std::cout << "location vector for chromosome " << i << std::endl;
        for (int j=0; j<N_SNPs[i]; ++j){
            std::cout << *(arr+SNP_ranges[i]+j) << ", ";
        }
        std::cout << std::endl;
    }
}

void print_annot_mats(float* arr, int N_chros, int* N_SNPs, int N_annots){
    int SNP_ranges[N_chros+1];
    gen_snp_ranges(&(SNP_ranges[0]), N_SNPs, N_chros);
    for (int i=0; i<N_chros; ++i){
        std::cout << "annotation matrix for chromosome " << i << std::endl;
        for (int j=0; j<N_SNPs[i]; ++j){
            for (int k=0; k<N_annots; ++k){
                #if ANNOT_MAT_ROW_OR_COLUMN_BASED == 1
                std::cout << *(arr+SNP_ranges[i]*N_annots+j*N_annots+k) << ", ";
                #else
                std::cout << *(arr+SNP_ranges[i]+j+SNP_ranges[N_annots+1]*k) << ", ";
                #endif
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}

void print_count_mats(float* arr, int N_chros, int N_times, int N_annots){
    for (int i=0; i<N_chros; ++i){
        std::cout << "count matrix for chromosome " << i << std::endl;
        for (int j=0; j<N_times; ++j){
            for (int k=0; k<N_annots; ++k){
                std::cout << *(arr+N_times*N_annots*i+j*N_annots+k) << ", ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}

void print_1D_F_array(float* array, int len){
    for (int i=0; i<len; ++i){
        std::cout<< *(array+i) << ", ";
    }
    std::cout << std::endl;
}

void gen_snp_ranges(int* SNP_ranges, int* N_SNPs, int N_chros){
    // Generate the SNP range for each chromosome
    // For example: N_SNPs = [4,6,10],
    // which means 4 SNP in chromosome 1, 6 in chromosome 2, 10 in chromosome 3
    // SNP_ranges will be: [0,4,10,20]
    // [0,4) is the range of chromosome 1
    // [4,10) is the range of chromosome 2
    // [10,20) is the range of chromosome 3
    int SNP_len = 0;
    for (int i=0; i<N_chros; ++i){
        SNP_ranges[i] = SNP_len;
        SNP_len+=N_SNPs[i];
    }
    SNP_ranges[N_chros] = SNP_len;
    #if IF_TEST == 1
    std::cout << "SNP range: " << std::endl;
    for (int i=0; i<N_chros+1; ++i){
        std::cout << SNP_ranges[i] << ", ";
    }
    std::cout << std::endl;
    #endif
}

void print_location_vecs_tofile(float* arr, int N_chros, int* N_SNPs, std::ofstream& outfile){
    int SNP_ranges[N_chros+1];
    gen_snp_ranges(&(SNP_ranges[0]), N_SNPs, N_chros);
    for (int i=0; i<N_chros; ++i){
        outfile << "location vector for chromosome " << i << std::endl;
        for (int j=0; j<N_SNPs[i]; ++j){
            outfile << *(arr+SNP_ranges[i]+j) << ", ";
        }
        outfile << std::endl;
    }
}

void print_annot_mats_tofile(float* arr, int N_chros, int* N_SNPs, int N_annots, std::ofstream& outfile){
    int SNP_ranges[N_chros+1];
    gen_snp_ranges(&(SNP_ranges[0]), N_SNPs, N_chros);
    for (int i=0; i<N_chros; ++i){
        outfile << "annotation matrix for chromosome " << i << std::endl;
        for (int j=0; j<N_SNPs[i]; ++j){
            for (int k=0; k<N_annots; ++k){
                #if ANNOT_MAT_ROW_OR_COLUMN_BASED == 1
                outfile << *(arr+SNP_ranges[i]*N_annots+j*N_annots+k) << ", ";
                #else
                outfile << *(arr+SNP_ranges[i]+j+SNP_ranges[N_annots+1]*k) << ", ";
                #endif
            }
            outfile << std::endl;
        }
        outfile << std::endl;
    }
}

void print_count_mats_tofile(float* arr, int N_chros, int N_times, int N_annots, std::ofstream& outfile){
    for (int i=0; i<N_chros; ++i){
        outfile << "count matrix for chromosome " << i << std::endl;
        for (int j=0; j<N_times; ++j){
            for (int k=0; k<N_annots; ++k){
                outfile << *(arr+N_times*N_annots*i+j*N_annots+k) << ", ";
            }
            outfile << std::endl;
        }
        outfile << std::endl;
    }
}