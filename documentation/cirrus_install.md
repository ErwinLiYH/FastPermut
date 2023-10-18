# Install to Cirrus

1. Setup python

   ```bash
   # set up python
   module load python/3.9.13 # load python
   python -m venv --system-site-packages /path/to/myvenv # create virture environment
   extend-venv-activate /path/to/myvenv
   
   source /path/to/myvenv/bin/activate # activate virture environment
   ```

   After setup python environment, you only need "source" it when login to Cirrus again.

   More information about Python env setup on Cirrus: [Cirrus Python documentation](https://cirrus.readthedocs.io/en/main/user-guide/python.html#installing-your-own-python-packages-with-pip)

2. Set C library and header files path

   ```bash
   export PERMUT_PATH=/path/to/library
   ```
	If not set this, default install path is `~/hpcPermut`
3. Install C library

   ```bash
   # Install MKL version
   module load openmpi # load openmpi (if not in python env)
   module load intel-20.4/cc # load intel-compiler
   module load intel-20.4/cmkl # load MKL
   export PERMUT_PATH=/path/to/library
   cd src/clib/mkl
   make
   make install
   
   # Install CUDA version
   module load nvidia/nvhpc-nompi/22.11 # load CUDA
   cd src/clib/cuda
   make
   make install
   ```

4. Install python package

   ```bash
   source /path/to/myvenv/bin/activate
   module load intel-20.4/cc
   module load intel-20.4/cmkl
   module load nvidia/nvhpc-nompi/22.11
   export PERMUT_PATH=/path/to/library
   export CC=gcc
   export CXX=g++
   cd src/snp-annot-permutation
   pip install .
   ```

   
