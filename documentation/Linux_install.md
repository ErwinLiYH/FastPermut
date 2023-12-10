# Install to General GNU/Linux

In this instruction, we use Debian 12 x86_64 as example.

1. Set up system environment
   1. MKL
      This project uses Intel MKL as CPU version's calculation library.
      Intel compiler `icpc` and MKL library `libmkl` are needed.
      So, only recommend using CPU version on the Intel platform.
      Please refer [Intel official instruction for Intel compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html) to install the Intel compiler and MKL library.

   2. OpenMPI

      Except for MKL, CPU version needs MPI for installing Python package `mpi4py`.

      ```bash
      apt update
      apt install libopenmpi-dev
      ```

   3. CUDA
      GPU version is based on Nvidia CUDA.
      You need a Nvidia GPU, then install the Nvidia GPU driver and CUDA tool kit.
      Nvidia compiler `nvcc` and library `cuBLAS`, `cuSPARSE`, `cuRAND` are needed.
         1. Install driver, [More information](https://wiki.debian.org/NvidiaGraphicsDrivers)
         ```bash
         # Add "contrib", "non-free" and "non-free-firmware" components to /etc/apt/sources.list
         deb http://deb.debian.org/debian/ bookworm main contrib non-free non-free-firmware

         # Update the list of available packages, then we can install the nvidia-driver package, plus the necessary firmware
         apt update
         apt install nvidia-driver firmware-misc-nonfree
         ```
         2. Install CUDA [More information](https://developer.nvidia.com/hpc-sdk-downloads)
         ```
         wget https://developer.download.nvidia.com/hpc-sdk/23.11/nvhpc_2023_2311_Linux_x86_64_cuda_12.3.tar.gz
         tar xpzf nvhpc_2023_2311_Linux_x86_64_cuda_12.3.tar.gz
         nvhpc_2023_2311_Linux_x86_64_cuda_12.3/install
         ```
   4. Install C/C++ compiler
      Used to compile Cython code when build and install Python package
      ```bash
      sudo apt-get install gcc g++
      ```

2. Set up Python environment

   Install conda:
   ```bash
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   sh Miniconda3-latest-Linux-x86_64.sh
   ```
   Restart the terminal, then set up a virtual environment:
   ```
   conda create -n fastpermut python=3.9
   conda activate fastpermut # activate virtual environment
   ```

   [More information](https://docs.conda.io/projects/miniconda/en/latest/) about conda.

3. Set C library and header files path

   ```bash
   export PERMUT_PATH=/path/to/library/parent # parent directory to library path
   ```
	If not set, default value is `~`

	For example: `~`, the library will be installed to `~/hpcPermut`

4. Install C library

   ```bash
   # set environment
   export PERMUT_PATH=/path/to/library/parent # set C library path(if not set)

   # Install MKL version
   cd src/clib/mkl
   make
   make install

   # Install CUDA version
   cd src/clib/cuda
   make
   make install
   ```

5. Install python package

   ```bash
   # set environment
   conda activate fastpermut # activate virtual environment(if not set)
   export PERMUT_PATH=/path/to/library/parent # set C library path(if not set)
   export CC=gcc
   export CXX=g++
   cd src/snp-annot-permutation

   # install both MKL(CPU) and CUDA(GPU) versions
   export INSTALL_SELECTOR=all # optional, if not set INSTALL_SELECTOR, it will install both by default
   pip install .

   # install only MKL(CPU) version
   export INSTALL_SELECTOR=cpu
   pip install .

   # install only CUDA(GPU) version
   export INSTALL_SELECTOR=gpu
   pip install .
   ```
