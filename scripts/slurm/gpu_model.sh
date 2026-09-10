#!/bin/bash
#SBATCH --job-name=2dEulerEquations
#SBATCH --account=ICT25_MHPC_0
#SBATCH --partition=boost_usr_prod    # Booster (GPU)
#SBATCH --qos=boost_qos_dbg           # Cola de debug (máx ~30 min)
#SBATCH --time=00:10:00               # 10 minutos para probar
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4           # 4 MPI ranks
#SBATCH --cpus-per-task=8             # 8 hilos CPU por rank
#SBATCH --gres=gpu:1                  # 1 GPU
#SBATCH --output=jacobi_%j.out
#SBATCH --error=jacobi_%j.err

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
