#!/bin/bash
#SBATCH -t 1:00:00     # walltime,
#SBATCH -N 1           # number of cluster nodes
#SBATCH -n 4     # number of total tasks

module load intel/19.1.3 openmpi/3.1.6

mpirun -np 1 ./pfdd_openmpi -in inBLi.fcc-pfdd > out-CuGB447-iso-edge-Sxz5
