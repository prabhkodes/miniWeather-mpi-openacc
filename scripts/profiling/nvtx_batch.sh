#!/bin/bash
#SBATCH --job-name=NVTX_Batch
#SBATCH --account=ICT25_MHPC_0
#SBATCH --partition=boost_usr_prod    # Booster (GPU)
#SBATCH --qos=boost_qos_dbg           # Cola de debug (máx ~30 min)
#SBATCH --time=00:10:00               # 10 minutos para probar
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4           # 1 MPI rank
#SBATCH --cpus-per-task=1             # 1 hilo CPU por rank
#SBATCH --gres=gpu:4                  # 4 GPUs
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

module purge
module load nvhpc/24.5
module load hpcx-mpi/2.19
module load netcdf-fortran/4.6.1--hpcx-mpi--2.19--nvhpc--24.5
module load binutils/2.42

echo "--- Building Model ---"
make clean || exit 1
make -j 4 USE_OPENACC=1 || exit 1

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close

echo "--- Running Simulation ---"
echo "Nodes = ${SLURM_NNODES}"
echo "Tasks / Node = ${SLURM_NTASKS_PER_NODE}"
echo "CPUs / Task = ${SLURM_CPUS_PER_TASK} (OMP Threads)"

srun nsys profile --nic-metrics=true \
    --trace=cuda,nvtx,mpi \
    -o "logs/%q{SLURM_JOB_ID}N%q{SLURM_JOB_NUM_NODES}%q{SLURM_PROCID}" ./model 100 1000 0

echo "--- Ending Simulation ---"
