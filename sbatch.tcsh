#!/bin/tcsh
#SBATCH --ntasks=96
#SBATCH --cpus-per-task=2

setenv OMP_NUM_THREADS 2

mpiexec ./body3 part_500.txt particles_out 200 0.5
