# HPC permutation test

> libpermut.so & libpermut_cuda.so

## Install to Cirrus

1. Install MKL versions

   ```bash
   module load intel-20.4/cc
   module load intel-20.4/cmkl
   module load openmpi
   cd clib/mkl
   make
   make install
   ```

2. Install CUDA version

   ```bash
   module load nvidia/nvhpc-nompi/22.11
   cd clib/cuda
   make
   make install
   ```