# SNP-Annotation permutation test

> python package snp-annot-permutation

## Install to Cirrus

1. Install the hpcPermut C library
2. Change the `include_dir` and `library_dir` if the `INSTALL_PATH` changed when installing hpcPermut. Default `INSTALL_PATH` is `~/hpcPermut`, default `include_dir` and `library_dir` are `~/hpcPermut/include` and `~/hpcPermut/lib`

3. Make a custom python environment and activate (see [Cirrus Python documentation](https://cirrus.readthedocs.io/en/main/user-guide/python.html#installing-your-own-python-packages-with-pip))

4. Install package

   ```bash
   module load intel-20.4/cc
   module load intel-20.4/cmkl
   module load nvidia/nvhpc-nompi/22.11
   export CC=gcc
   export CXX=g++
   cd snp-annot-permutation
   pip install .
   ```

## Example codes

see [test code](./test)
