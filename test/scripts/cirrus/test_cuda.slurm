#!/bin/bash

# This script is configured for rapid turnaround in the short queue
# This only works for jobs *less than 20 minutes on 2 nodes or fewer*.
# For more nodes set "--qos=standard" and delete the reservation line.

#SBATCH --account=your_account_name
#SBATCH --job-name=test_cuda
#SBATCH --time=00:01:00
#SBATCH --nodes=1
# #SBATCH --ntasks=4
#SBATCH --qos=short
#SBATCH --gres=gpu:1

# You can probably leave these options mostly as they are

#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --exclusive
#SBATCH --partition=gpu

# used to get the runtime of sparse version code against the number of threads

source py/env/path/bin/activate

python ../../codes/test_cuda.py