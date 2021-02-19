#!/bin/bash -l
#SBATCH -p short
#SBATCH -t 00:49:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=4G
#SBATCH -o SW_1d_MyEq_gamma05.out
 
module load matlab/r2016b
srun matlab_multithread -nosplash -r "SW_1D_parfor_myequation_gamma05(24); exit(0)"
