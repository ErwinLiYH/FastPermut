#!/bin/bash

# This script is configured for rapid turnaround in the short queue
# This only works for jobs *less than 20 minutes on 2 nodes or fewer*.
# For more nodes set "--qos=standard" and delete the reservation line.

#SBATCH --account=your_account_name
#SBATCH --job-name=mkl_test
#SBATCH --time=00:01:00
#SBATCH --nodes=1
# #SBATCH --ntasks=4
#SBATCH --qos=standard

# You can probably leave these options mostly as they are

#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --exclusive
#SBATCH --tasks-per-node=36
#SBATCH --cpus-per-task=1
#SBATCH --partition=standard

# used to get the runtime of sparse version code against the number of permutation steps

source py/env/path/bin/activate
module load intel-20.4/cc
module load intel-20.4/cmkl
export LD_PRELOAD=/mnt/lustre/indy2lfs/sw/intel/compilers_and_libraries_2020.4.304/linux/mkl/lib/intel64_lin/libmkl_def.so:/mnt/lustre/indy2lfs/sw/intel/compilers_and_libraries_2020.4.304/linux/mkl/lib/intel64_lin/libmkl_avx2.so:/mnt/lustre/indy2lfs/sw/intel/compilers_and_libraries_2020.4.304/linux/mkl/lib/intel64_lin/libmkl_core.so:/mnt/lustre/indy2lfs/sw/intel/compilers_and_libraries_2020.4.304/linux/mkl/lib/intel64_lin/libmkl_intel_lp64.so:/mnt/lustre/indy2lfs/sw/intel/compilers_and_libraries_2020.4.304/linux/mkl/lib/intel64_lin/libmkl_intel_thread.so:/mnt/lustre/indy2lfs/sw/intel/compilers_and_libraries_2020.4.304/linux/compiler/lib/intel64_lin/libiomp5.so

python ../../codes/test_mkl.py