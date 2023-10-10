#include "permutation.cuh"

void cuda_permut(unsigned long long* seed_p, float* location_vecs, int* N_SNPs, int N_chros, float* annot_mats, int N_annots, int N_times, float* count_mats){
    //>>> timer------------------
    #if IF_TIMING==1
    struct timeval start0, end0;
    gettimeofday(&start0, NULL);
    #endif
    //<<<------------------------

    // Generate the SNP range for each chromosome---------------------------------
    int SNP_ranges[N_chros+1];
    gen_snp_ranges(&(SNP_ranges[0]), N_SNPs, N_chros);
    //----------------------------------------------------------------------------

    //Create uniform distributions for all chromosomes---------------------------
    float* d_RNs;
    CHECK_CUDA( cudaMalloc(&d_RNs, N_chros*N_times*sizeof(float)) )
    for (int i=0; i<N_chros; ++i){
        if (seed_p == NULL) {
            unsigned long long seed = std::chrono::system_clock::now().time_since_epoch().count();
            gen_random_cuda(d_RNs+N_times*i, N_times, seed);
        }else{
            gen_random_cuda(d_RNs+N_times*i, N_times, *seed_p*i);
        }
    }

    //>>> store the random values in IF_STORE_RNS mode----
    #if IF_STORE_RNS==1
    float* RNs = new float[N_chros*N_times];
    CHECK_CUDA( cudaMemcpy(RNs, d_RNs, N_chros*N_times*sizeof(float), cudaMemcpyDeviceToHost) )
    std::ofstream outfile("cuda_permut_RNs", std::ios::binary);
    for (int i=0; i<N_chros; ++i){
        for (int j=0; j<N_times; ++j){
            int RN = std::floor(*(RNs+N_times*i+j)*N_SNPs[i]);
            outfile.write(reinterpret_cast<char*>(&RN), sizeof(int));
        }
    }
    outfile.close();
    delete[] RNs;
    #endif
    //<<<-------------------------------------------------

    //----------------------------------------------------------------------------

    // copy memory to device------------------------------------------------------
    float* d_location_vecs;
    CHECK_CUDA( cudaMalloc((void**)&d_location_vecs, SNP_ranges[N_chros]*sizeof(float)) )
    CHECK_CUDA( cudaMemcpy(d_location_vecs, location_vecs, SNP_ranges[N_chros]*sizeof(float), cudaMemcpyHostToDevice) )
    float* d_annot_mats;
    CHECK_CUDA( cudaMalloc((void**)&d_annot_mats, N_annots*SNP_ranges[N_chros]*sizeof(float)) )
    CHECK_CUDA( cudaMemcpy(d_annot_mats, annot_mats, N_annots*SNP_ranges[N_chros]*sizeof(float), cudaMemcpyHostToDevice) )
    float* d_count_mats;
    CHECK_CUDA( cudaMalloc((void**)&d_count_mats, N_chros*N_times*N_annots*sizeof(float)) )
    float* d_location_vecs_rotated;
    CHECK_CUDA( cudaMalloc((void**)&d_location_vecs_rotated, SNP_ranges[N_chros]*sizeof(float)) )
    //----------------------------------------------------------------------------

    // Create a handle for CUBLAS
    cublasHandle_t handle;
    cublasCreate(&handle);

    // alpha and bate used in cublas
    const float a=1.0;
    const float b=0.0;

    //>>> timer------------------
    #if IF_TIMING==1
    struct timeval start, end;
    gettimeofday(&start, NULL);
    #endif
    //<<<------------------------

    // do permutation-------------------------------------------------------------
    for (int i=0; i< N_chros; ++i){// chromosomes loop
        for (int j=0; j<N_times; ++j){ // permutations loop

            //>>> print process in TEST mode---------
            #if IF_TEST == 1
            if (j%PRINT_INTERVAL==0){
                std::cout<< "process chromosome " << i << ", " << "permutation " << j << std::endl;
            }
            #endif
            //<<<------------------------------------

            //>>> print origianl location_vecs and radom number for rotate in TEST_PRINT_ARRAY mode
            #if TEST_PRINT_ARRAY == 1
            float random_number;
            CHECK_CUDA( cudaMemcpy(&random_number, d_RNs+N_times*i+j, sizeof(float), cudaMemcpyDeviceToHost) )
            std::cout << "random number: " << std::floor(random_number*N_SNPs[i]) << std::endl;
            std::cout << "Location vectors before rotate" << std::endl;
            print_location_vecs(location_vecs, N_chros, N_SNPs);
            #endif
            //<<<-------------------------------------------

            // rotate it in GPU
            rotate_to_Kernel<<<(N_SNPs[i]+255)/256, 256>>>(d_location_vecs+SNP_ranges[i], d_location_vecs_rotated+SNP_ranges[i], N_SNPs[i], d_RNs+N_times*i+j);

            //>>> print rotaed location_vecs in TEST_PRINT_ARRAY mode-----
            #if TEST_PRINT_ARRAY == 1
            // cudaDeviceSynchronize();
            float* rotated = new float[SNP_ranges[N_chros]];
            CHECK_CUDA( cudaMemcpy(rotated, d_location_vecs_rotated, SNP_ranges[N_chros]*sizeof(float), cudaMemcpyDeviceToHost) )
            // std::cout << "rotated copy: " << std::endl;
            // print_1D_F_array(rotated, N_SNPs[i]);
            print_location_vecs(rotated, N_chros, N_SNPs);
            delete[] rotated;
            #endif
            //<<<-----------------------------------------------------------

            // calculate the matrix times vector by cuBLAS------------------------------------------------------
            #if ANNOT_MAT_ROW_OR_COLUMN_BASED == 1
            // d_annot_matsin row-based memory layout
            cublasSgemv(handle, CUBLAS_OP_N, N_annots, N_SNPs[i], &a, d_annot_mats+SNP_ranges[i]*N_annots, N_annots, d_location_vecs_rotated+SNP_ranges[i], 1, &b, d_count_mats+i*N_times*N_annots+j*N_annots, 1);
            #else
            // d_annot_matsin column-based memory layout
            cublasSgemv(handle, CUBLAS_OP_T, N_SNPs[i], N_annots, &a, d_annot_mats+SNP_ranges[i], SNP_ranges[N_chros], d_location_vecs_rotated+SNP_ranges[i], 1, &b, d_count_mats+i*N_times*N_annots+j*N_annots, 1);
            #endif
            //---------------------------------------------------------------------------------------------------

            //>>> print count matrixes in TEST_PRINT_ARRAY mode-----
            #if TEST_PRINT_ARRAY == 1
            CHECK_CUDA( cudaMemcpy(count_mats, d_count_mats, N_chros*N_times*N_annots*sizeof(float), cudaMemcpyDeviceToHost) )
            std::cout << "count matrixes: " << i << std::endl;
            print_count_mats(count_mats, N_chros, N_times, N_annots);
            #endif
            //<<<---------------------------------------------------
        }
    }
    //----------------------------------------------------------------------------

    //>>> timer------------------
    #if IF_TIMING==1
    gettimeofday(&end, NULL);
    #endif
    //<<<------------------------

    //destroy cublas handel
    cublasDestroy(handle);

    //copy count matrix back
    cudaMemcpy(count_mats, d_count_mats, N_chros*N_times*N_annots*sizeof(float), cudaMemcpyDeviceToHost);

    //free cuda memory
    cudaFree(d_location_vecs);
    cudaFree(d_annot_mats);
    cudaFree(d_count_mats);
    cudaFree(d_location_vecs_rotated);

    //>>> timer------------------
    #if IF_TIMING==1
    gettimeofday(&end0, NULL);
    double time_taken, time_taken0;
    time_taken = (end.tv_sec - start.tv_sec) * 1e6;
    time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;
    time_taken0 = (end0.tv_sec - start0.tv_sec) * 1e6;
    time_taken0 = (time_taken0 + (end0.tv_usec - start0.tv_usec)) * 1e-6;
    if (ANNOT_MAT_ROW_OR_COLUMN_BASED==0){
        std::cout << "Column-major variant" << std::endl;
    }else{
        std::cout << "Row-major variant" << std::endl;
    }
    std::cout << "cuda_permut total time: "<< std::fixed << std::setprecision(10) << time_taken0 << std::endl;
    std::cout << "cuda_permut permutation time: " << std::fixed << std::setprecision(10) << time_taken << std::endl;
    std::cout.unsetf(std::ios_base::floatfield);
    std::cout << std::setprecision(6);
    #endif
    //<<<------------------------
}

void cuda_permut_sparse(unsigned long long* seed_p, float* location_vecs, int* N_SNPs, int N_chros, float* annot_mats, int N_annots, int N_times, float* count_mats){
    //>>> timer------------------
    #if IF_TIMING==1
    struct timeval start0, end0;
    gettimeofday(&start0, NULL);
    #endif
    //<<<------------------------

    // Generate the SNP range for each chromosome---------------------------------
    int SNP_ranges[N_chros+1];
    gen_snp_ranges(&(SNP_ranges[0]), N_SNPs, N_chros);
    //----------------------------------------------------------------------------

    // Create uniform distributions for all chromosomes---------------------------
    float* d_RNs;
    CHECK_CUDA( cudaMalloc(&d_RNs, N_chros*N_times*sizeof(float)) )
    for (int i=0; i<N_chros; ++i){
        if (seed_p == NULL) {
            unsigned long long seed = std::chrono::system_clock::now().time_since_epoch().count();
            gen_random_cuda(d_RNs+N_times*i, N_times, seed);
        }else{
            gen_random_cuda(d_RNs+N_times*i, N_times, *seed_p*i);
        }
    }

    //>>> store the random values in IF_STORE_RNS mode----
    #if IF_STORE_RNS==1
    float* RNs = new float[N_chros*N_times];
    CHECK_CUDA( cudaMemcpy(RNs, d_RNs, N_chros*N_times*sizeof(float), cudaMemcpyDeviceToHost) )
    std::ofstream outfile("cuda_permut_sparse_RNs", std::ios::binary);
    for (int i=0; i<N_chros; ++i){
        for (int j=0; j<N_times; ++j){
            int RN = std::floor(*(RNs+N_times*i+j)*N_SNPs[i]);
            outfile.write(reinterpret_cast<char*>(&RN), sizeof(int));
        }
    }
    outfile.close();
    delete[] RNs;
    #endif
    //<<<-------------------------------------------------

    //----------------------------------------------------------------------------

    // copy memory to device------------------------------------------------------
    float* d_location_vecs;
    CHECK_CUDA( cudaMalloc((void**)&d_location_vecs, SNP_ranges[N_chros]*sizeof(float)) )
    CHECK_CUDA( cudaMemcpy(d_location_vecs, location_vecs, SNP_ranges[N_chros]*sizeof(float), cudaMemcpyHostToDevice) )
    float* d_count_mats;
    CHECK_CUDA( cudaMalloc((void**)&d_count_mats, N_chros*N_times*N_annots*sizeof(float)) )
    float* d_location_vecs_rotated;
    CHECK_CUDA( cudaMalloc((void**)&d_location_vecs_rotated, SNP_ranges[N_chros]*sizeof(float)) )
    //----------------------------------------------------------------------------

    //Set up sparse matrixes and dense vector----------------------------------------
    size_t bufferSize = 0;
    void* dBuffer;

    cusparseDnVecDescr_t dense_location_vec[N_chros];
    for (int i=0; i<N_chros; ++i){
        CHECK_CUSPARSE( cusparseCreateDnVec(&(dense_location_vec[i]), N_SNPs[i], d_location_vecs_rotated+SNP_ranges[i], CUDA_R_32F) )
    }

    int size_nnz[N_chros];
    int size_p[N_chros];

    float* values[N_chros];
    float* values_d[N_chros];

    int* indexes[N_chros];
    int* indexes_d[N_chros];

    int* pointers[N_chros];
    int* pointers_d[N_chros];

    cusparseSpMatDescr_t smats[N_chros];
    for (int i=0; i<N_chros; ++i){
        #if ANNOT_MAT_ROW_OR_COLUMN_BASED == 1
        gen_sparse_mat(annot_mats+SNP_ranges[i]*N_annots, N_annots, N_SNPs[i], N_annots, &(values[i]), &(indexes[i]), &(pointers[i]), &(size_nnz[i]), &(size_p[i]));
        #else
        gen_sparse_mat(annot_mats+SNP_ranges[i], N_annots, N_SNPs[i], SNP_ranges[N_chros], &(values[i]), &(indexes[i]), &(pointers[i]), &(size_nnz[i]), &(size_p[i]));
        #endif

        CHECK_CUDA( cudaMalloc((void**)&(values_d[i]), size_nnz[i]*sizeof(float)) )
        CHECK_CUDA( cudaMemcpy(values_d[i], values[i], size_nnz[i]*sizeof(float), cudaMemcpyHostToDevice) )

        CHECK_CUDA( cudaMalloc((void**)&(indexes_d[i]), size_nnz[i]*sizeof(int)) )
        CHECK_CUDA( cudaMemcpy(indexes_d[i], indexes[i], size_nnz[i]*sizeof(int), cudaMemcpyHostToDevice) )

        CHECK_CUDA( cudaMalloc((void**)&(pointers_d[i]), size_p[i]*sizeof(int)) )
        CHECK_CUDA( cudaMemcpy(pointers_d[i], pointers[i], size_p[i]*sizeof(int), cudaMemcpyHostToDevice) )

        #if ANNOT_MAT_ROW_OR_COLUMN_BASED == 1
        CHECK_CUSPARSE (cusparseCreateCsr(&(smats[i]), N_SNPs[i], N_annots, size_nnz[i], pointers_d[i], indexes_d[i], values_d[i], CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F) )
        #else
        CHECK_CUSPARSE( cusparseCreateCsc(&(smats[i]), N_SNPs[i], N_annots, size_nnz[i], pointers_d[i], indexes_d[i], values_d[i], CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F) )
        #endif
    }
    //----------------------------------------------------------------------------

    // Create a handle for CUBLAS
    cusparseHandle_t handle;
    CHECK_CUSPARSE( cusparseCreate(&handle) )

    // alpha and bate used in cublas
    const float a=1.0;
    const float b=0.0;

    //>>> timer------------------
    #if IF_TIMING==1
    struct timeval start, end;
    gettimeofday(&start, NULL);
    #endif
    //<<<------------------------

    // do permutation-------------------------------------------------------------
    for (int i=0; i< N_chros; ++i){// chromosomes loop
        for (int j=0; j<N_times; ++j){ // permutations loop

            //>>> print process in TEST mode---------
            #if IF_TEST == 1
            if (j%PRINT_INTERVAL==0){
                std::cout<< "process chromosome " << i << ", " << "permutation " << j << std::endl;
            }
            #endif
            //<<<------------------------------------

            //>>> print origianl location_vecs and radom number for rotate in TEST_PRINT_ARRAY mode
            #if TEST_PRINT_ARRAY == 1
            float random_number;
            CHECK_CUDA( cudaMemcpy(&random_number, d_RNs+N_times*i+j, sizeof(float), cudaMemcpyDeviceToHost) )
            std::cout << "random number: " << std::floor(random_number*N_SNPs[i]) << "/" << random_number << std::endl;
            std::cout << "Location vectors before rotate" << std::endl;
            CHECK_CUDA( cudaMemcpy(location_vecs, d_location_vecs, SNP_ranges[N_chros]*sizeof(float), cudaMemcpyDeviceToHost) )
            print_location_vecs(location_vecs, N_chros, N_SNPs);
            #endif
            //<<<-------------------------------------------

            // rotate it in GPU
            rotate_to_Kernel<<<(N_SNPs[i]+255)/256, 256>>>(d_location_vecs+SNP_ranges[i], d_location_vecs_rotated+SNP_ranges[i], N_SNPs[i], d_RNs+N_times*i+j);

            //>>> test the status of kernel function--------------
            #if IF_TEST
            cudaDeviceSynchronize();
            cudaError_t err = cudaGetLastError();
            if (err != cudaSuccess) {
                printf("Error: %s\n", cudaGetErrorString(err));
            }
            cudaDeviceSynchronize();
            #endif
            //<<<-------------------------------------------------

            //>>> print rotaed location_vecs in TEST_PRINT_ARRAY mode-----
            #if TEST_PRINT_ARRAY == 1
            float* rotated = new float[SNP_ranges[N_chros]];
            CHECK_CUDA( cudaMemcpy(rotated, d_location_vecs_rotated, SNP_ranges[N_chros]*sizeof(float), cudaMemcpyDeviceToHost) )
            print_location_vecs(rotated, N_chros, N_SNPs);
            delete[] rotated;
            #endif
            //<<<-----------------------------------------------------------

            // calculate the matrix times vector by cuBLAS------------------------------------------------------
            cusparseDnVecDescr_t count_vec;
            CHECK_CUSPARSE( cusparseCreateDnVec(&count_vec, N_annots, d_count_mats+i*N_times*N_annots+j*N_annots, CUDA_R_32F) )
            if (j==0){
                CHECK_CUSPARSE( cusparseSpMV_bufferSize(handle, CUSPARSE_OPERATION_TRANSPOSE, &a, smats[i], dense_location_vec[i], &b, count_vec, CUDA_R_32F, CUSPARSE_SPMV_ALG_DEFAULT, &bufferSize) )
                #if IF_TEST
                std::cout << "bufferSize: " << bufferSize << std::endl;
                #endif
                if (bufferSize!=0){
                    CHECK_CUDA( cudaMalloc(&dBuffer, bufferSize) )
                }
            }
            CHECK_CUSPARSE( cusparseSpMV(handle, CUSPARSE_OPERATION_TRANSPOSE, &a, smats[i], dense_location_vec[i], &b, count_vec, CUDA_R_32F, CUSPARSE_SPMV_ALG_DEFAULT, dBuffer) )
            CHECK_CUSPARSE( cusparseDestroyDnVec(count_vec) )
            //---------------------------------------------------------------------------------------------------

            //>>> print count matrixes in TEST_PRINT_ARRAY mode-----
            #if TEST_PRINT_ARRAY == 1
            CHECK_CUDA( cudaMemcpy(count_mats, d_count_mats, N_chros*N_times*N_annots*sizeof(float), cudaMemcpyDeviceToHost) )
            std::cout << "count matrixes: " << i << std::endl;
            print_count_mats(count_mats, N_chros, N_times, N_annots);
            #endif
            //<<<---------------------------------------------------
        }
    }
    //----------------------------------------------------------------------------

    //copy count matrix back
    cudaMemcpy(count_mats, d_count_mats, N_chros*N_times*N_annots*sizeof(float), cudaMemcpyDeviceToHost);

    //>>> timer------------------
    #if IF_TIMING==1
    gettimeofday(&end, NULL);
    #endif
    //<<<------------------------

    //Destroy sparse matrixes and dense vectors
    for (int i=0; i<N_chros; ++i){
        cusparseDestroySpMat(smats[i]);
        cusparseDestroyDnVec(dense_location_vec[i]);
        delete[] values[i];
        cudaFree(values_d[i]);
        delete[] indexes[i];
        cudaFree(indexes_d[i]);
        delete[] pointers[i];
        cudaFree(pointers_d[i]);
    }

    //free cuda memory
    cudaFree(d_location_vecs);
    cudaFree(d_count_mats);
    cudaFree(d_location_vecs_rotated);
    cudaFree(dBuffer);

    //Destroy handle
    cusparseDestroy(handle);

    //>>> timer------------------
    #if IF_TIMING==1
    gettimeofday(&end0, NULL);
    double time_taken, time_taken0;
    time_taken = (end.tv_sec - start.tv_sec) * 1e6;
    time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;
    time_taken0 = (end0.tv_sec - start0.tv_sec) * 1e6;
    time_taken0 = (time_taken0 + (end0.tv_usec - start0.tv_usec)) * 1e-6;
    if (ANNOT_MAT_ROW_OR_COLUMN_BASED==0){
        std::cout << "Column-major variant" << std::endl;
    }else{
        std::cout << "Row-major variant" << std::endl;
    }
    std::cout << "cuda_permut_sparse total time: "<< std::fixed << std::setprecision(10) << time_taken0 << std::endl;
    std::cout << "cuda_permut_sparse permutation time: " << std::fixed << std::setprecision(10) << time_taken << std::endl;
    std::cout.unsetf(std::ios_base::floatfield);
    std::cout << std::setprecision(6);
    #endif
    //<<<------------------------
}