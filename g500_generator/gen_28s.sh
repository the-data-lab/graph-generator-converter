#!/bin/sh
# one hour timelimit:
#SBATCH --time 01:00:00

# default queue, 512 processors (two nodes worth)
#SBATCH -p defq -n 4 
# -N number of machines used in this app
# -n total number of processes used in this app
module load openmpi

mpirun ./g500_2d_tuple/generator_test_mpi 28 16 1 1 

