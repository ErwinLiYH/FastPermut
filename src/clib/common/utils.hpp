#include "config.hpp"
#include <iostream>
#include <fstream>

void test_order(float* arr, int size);

void print_location_vecs(float* arr, int N_chros, int* N_SNPs);

void print_annot_mats(float* arr, int N_chros, int* N_SNPs, int N_annots);

void print_count_mats(float* arr, int N_chros, int N_times, int N_annots);

void print_1D_F_array(float* array, int len);

void gen_snp_ranges(int* SNP_ranges, int* N_SNPs, int N_chros);

void print_location_vecs_tofile(float* arr, int N_chros, int* N_SNPs, std::ofstream& outfile);

void print_annot_mats_tofile(float* arr, int N_chros, int* N_SNPs, int N_annots, std::ofstream& outfile);

void print_count_mats_tofile(float* arr, int N_chros, int N_times, int N_annots, std::ofstream& outfile);