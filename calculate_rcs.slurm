#!/bin/bash


#SBATCH -D ./
#SBATCH -J matlab_pct_batch_job
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --nodes=1
#SBATCH --tasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --time=12:30:00
 
module load slurm_setup
module switch spack/23.1.0
 
# IMPORTANT
# Default settings of Intel MPI module may disrupt
# functionality of Parallel-Computing-Toolbox!
# Do one of the following solutions: 
# (1) Unload Intel MPI module:
module rm intel-mpi
module rm mpi.intel
# (2) If Intel MPI module is mandatory, uncomment next 2 lines
module load intel-mpi/2018.4.274
# export KMP_AFFINITY=granularity=thread,none
 
module load matlab/R2023b_Update4-generic # EDIT HERE (see supported releases)
 
 
# Example: matrix-matrix multiplication C = A*B
#          with A of size NROWA x NCOLA and
#          B of size NROWB x NCOLB
dir=./Bench_100/*/*
interval_start=1
interval_end=100
frequency=94e9
 
# Run MATLAB
# => Using option -r don't add file extension .m to the function call!
# => MATLAB commandline arguments are case-sensitive!
matlab -nodisplay -singleCompThread \
       -r "calculate_rcs('$dir' ,[$interval_start $interval_end], $frequency);"