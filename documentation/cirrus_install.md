# Install to Cirrus

[Cirrus](https://www.cirrus.ac.uk/) is an HPC system located at the University of Edinburgh designed to solve computational, simulation, modelling, and data science challenges.

The following instructions describe how to setup `FastPermut` on this system.

1. Setup Python

   ```bash
   # set up python
   module load python/3.9.13 # load python
   python -m venv --system-site-packages /path/to/myvenv # create virtual environment
   extend-venv-activate /path/to/myvenv

   source /path/to/myvenv/bin/activate # activate virtual environment
   ```

   After the Python environment has been setup, you only need "source" the virtual environment when you login to Cirrus again.

   More information about Python env setup on Cirrus: [Cirrus Python documentation](https://docs.cirrus.ac.uk/user-guide/python/).

2. Set the C library and header files path

   ```bash
   export PERMUT_PATH=/path/to/library/parent # parent directory to library path
   ```
	If it is not set, the default value will be  `~` (your home environment).

	We recommend you use `/work/user_name/` in Cirrus. For example: `/work/m22oc/m22oc/s2303098`, the library will be installed to `/work/m22oc/m22oc/s2303098/hpcPermut`.
4. Install the C library

   ```bash
   # set environment
   export PERMUT_PATH=/path/to/library/parent # set C library path(if not set)

   # Install both MKL and CUDA versions
   module load openmpi # load openmpi (if not in python env)
   # in Cirrus, python/3.9.13 already included openmpi
   module load intel-20.4/cc # load intel-compiler
   module load intel-20.4/cmkl # load MKL
   module load nvidia/nvhpc-nompi/22.11 # load CUDA
   cd src/clib
   make
   make install

   # Or install only MKL version
   module load openmpi # load openmpi (if not in python env)
   # in Cirrus, python/3.9.13 already included openmpi
   module load intel-20.4/cc # load intel-compiler
   module load intel-20.4/cmkl # load MKL
   cd src/clib
   make cpu
   make install_cpu

   # Or install only CUDA version
   module load nvidia/nvhpc-nompi/22.11 # load CUDA
   cd src/clib
   make gpu
   make install_gpu
   ```

5. Install the Python package

   The version of Python package must be same as the version of C library installed. For exmple, if only CPU version of C library installed, only CPU version of Python package can be installed.

   ```bash
   # set environment
   source /path/to/myvenv/bin/activate # activate virture environment(if not set)
   export PERMUT_PATH=/path/to/library/parent # set C library path(if not set)
   cd src/FastPermut/

   # Install both MKL(CPU) and CUDA(GPU) versions
   module load intel-20.4/cc # load intel-compiler(if not set)
   module load intel-20.4/cmkl # load MKL(if not set)
   module load nvidia/nvhpc-nompi/22.11 # load CUDA(if not set)
   export INSTALL_SELECTOR=all # optional, if not set INSTALL_SELECTOR, it will install both by default

   # Or install only MKL(CPU) version
   module load intel-20.4/cc # load intel-compiler(if not set)
   module load intel-20.4/cmkl # load MKL(if not set)
   export INSTALL_SELECTOR=cpu

   # Or install only CUDA(GPU) version
   module load nvidia/nvhpc-nompi/22.11 # load CUDA(if not set)
   export INSTALL_SELECTOR=gpu

   # set compilor for Cython
   export CC=gcc
   export CXX=g++

   # Install
   pip install .
   ```
