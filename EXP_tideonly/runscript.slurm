#!/bin/bash
#SBATCH --job-name=amm7-surge
#SBATCH --time=00:05:00
#SBATCH --account=n01-CLASS
#SBATCH --partition=standard
#SBATCH --qos=standard

#SBATCH --nodes=2
#SBATCH --ntasks-per-core=1

# load modules
#. /work/n01/n01/jdha/modules/gnu-mpich
. /work/n01/n01/annkat/SE-NEMO_UPD/SE-NEMO/scripts/env/gnu-mpich

export OMP_NUM_THREADS=1
echo `date` : Launch Job

#SRUN_CMD=`./slurm_setup -S 1 -s 1 -m 1 -C 95 -g 999 -N 128 -t 00:05:00 -a n01-CLASS -j amm7-surge`  # set nodes=1
SRUN_CMD=`./slurm_setup -S 4 -s 16 -m 2 -C 96 -g 2 -N 128 -t 00:05:00 -a n01-CLASS -j amm7-surge -v`  # set nodes=2

chmod u+x ./myscript_wrapper.sh

echo $SRUN_CMD
eval $SRUN_CMD
