#!/bin/bash
#SBATCH --job-name=AMM7-S
#SBATCH --partition=compute
#SBATCH --time=00:30:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks-per-socket=32
#SBATCH --nodes=2


module load NEMO/prg-env

echo $PWD

export OMP_NUM_THREADS=1
export I_MPI_SHM=off
unset I_MPI_PMI_LIBRARY

mpiexec.hydra -print-rank-map -ppn 1 -np 16 ./xios_server.exe : -np 96 ./nemo

exit
