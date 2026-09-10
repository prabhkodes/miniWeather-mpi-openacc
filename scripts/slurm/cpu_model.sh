#!/bin/bash
#SBATCH --job-name=2dEulerEquations
#SBATCH --account=ICT25_MHPC
#SBATCH --partition=dcgp_usr_prod     # DCGP (CPU)
#SBATCH --qos=normal                  # Queue
#SBATCH --time=00:10:00               # 10 minutes
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4           # 4 MPI ranks
#SBATCH --cpus-per-task=8             # 8 threads per rank
#SBATCH --output=2dEulerEquations_%j.out
#SBATCH --error=2dEulerEquations_%j.err

module purge

# Compiler
module load gcc/12.2.0

# MPI
module load openmpi/4.1.6--gcc--12.2.0-cuda-12.2

# NetCDF Fortran compile with MPI + GCC
module load netcdf-fortran/4.6.1--openmpi--4.1.6--gcc--12.2.0-spack0.22

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close

srun ./serial/model
