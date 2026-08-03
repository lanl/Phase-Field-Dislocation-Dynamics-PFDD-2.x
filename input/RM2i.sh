#!/bin/bash
#SBATCH -t 3:00:00     # walltime,
#SBATCH -N 1           # number of cluster nodes
#SBATCH -n 8     # number of total tasks

module load intel/19.1.3 openmpi/3.1.6

mpirun -np 8 ./pfdd_openmpi -in inM2i.fcc-pfdd > out-NonCuGB472SlabY-iso-mix-Sxz0077-a202-core158ScaledMu
