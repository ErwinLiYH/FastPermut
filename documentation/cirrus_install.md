# Install to Cirrus

1. Setup Python

   ```bash
   # set up python
   module load python/3.9.13 # load python
   python -m venv --system-site-packages /path/to/myvenv # create virtual environment
   extend-venv-activate /path/to/myvenv

   source /path/to/myvenv/bin/activate # activate virture environment
   ```

   After setup the Python environment, you only need "source" it when login to Cirrus again.

   More information about Python env setup on Cirrus: [Cirrus Python documentation](https://cirrus.readthedocs.io/en/main/user-guide/python.html#installing-your-own-python-packages-with-pip)

2. Set C library and header files path

   ```bash
   export PERMUT_PATH=/path/to/library/parent # parent directory to library path
   ```
	If not set, the default value is `~`

	Recommend to use `/work/user_name/` in Cirrus. For example: `/work/m22oc/m22oc/s2303098`, the library will be installed to `/work/m22oc/m22oc/s2303098/hpcPermut`
4. Install C library

   ```bash
   # set environment
   export PERMUT_PATH=/path/to/library/parent # set C library path(if not set)

   # Install MKL version
   module load openmpi # load openmpi (if not in python env)
   # in Cirrus, python/3.9.13 already included openmpi
   module load intel-20.4/cc # load intel-compiler
   module load intel-20.4/cmkl # load MKL
   cd src/clib/mkl
   make
   make install

   # Install CUDA version
   module load nvidia/nvhpc-nompi/22.11 # load CUDA
   cd src/clib/cuda
   make
   make install
   ```

5. Install python package

   ```bash
   # set environment
   source /path/to/myvenv/bin/activate # activate virture environment(if not set)
   export PERMUT_PATH=/path/to/library/parent # set C library path(if not set)
   export CC=gcc
   export CXX=g++
   cd src/snp-annot-permutation

   # install both MKL(CPU) and CUDA(GPU) versions
   module load intel-20.4/cc # load intel-compiler(if not set)
   module load intel-20.4/cmkl # load MKL(if not set)
   module load nvidia/nvhpc-nompi/22.11 # load CUDA(if not set)
   export INSTALL_SELECTOR=all # optional, if not set INSTALL_SELECTOR, it will install both by default
   pip install .

   # install only MKL(CPU) version
   module load intel-20.4/cc # load intel-compiler(if not set)
   module load intel-20.4/cmkl # load MKL(if not set)
   export INSTALL_SELECTOR=cpu
   pip install .

   # install only CUDA(GPU) version
   module load nvidia/nvhpc-nompi/22.11 # load CUDA(if not set)
   export INSTALL_SELECTOR=gpu
   pip install .
   ```
