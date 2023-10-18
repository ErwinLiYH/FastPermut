# Install to Cirrus

1. Setup python

   ```bash
   # set up python
   module load python/3.9.13 # load python
   python -m venv --system-site-packages /path/to//myvenv # create virture environment
   extend-venv-activate /path/to//myvenv
   source /path/to//myvenv # activate virture environment
   ```

   After setup python environment, you only need "source" it when login to Cirrus again.

   More information about Python env setup on Cirrus: [Cirrus Python documentation](https://cirrus.readthedocs.io/en/main/user-guide/python.html#installing-your-own-python-packages-with-pip)

2. Install C library

   ```bash
   # Install MKL version
   module load openmpi # load openmpi (if not in python env)
   module load intel-20.4/cc # load intel-compiler
   module load intel-20.4/cmkl # load MKL
   cd clib/mkl
   make
   make install
   
   # Install CUDA version
   module load nvidia/nvhpc-nompi/22.11 # load CUDA
   cd clib/cuda
   make
   make install
   ```

3. Install python package

   ```bash
   module load intel-20.4/cc
   module load intel-20.4/cmkl
   module load nvidia/nvhpc-nompi/22.11
   export CC=gcc
   export CXX=g++
   cd snp-annot-permutation
   pip install .
   ```

   
