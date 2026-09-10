#!/usr/bin/env bash
set -euo pipefail

# Sweep launcher for CPU (batch.sh) and GPU (gpu_batch.sh).
# Configures node lists by type and submits sbatch with minimal overrides.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CPU_SCRIPT="${ROOT_DIR}/batch.sh"
GPU_SCRIPT="${ROOT_DIR}/gpu_batch.sh"

# Node lists to test (space-separated). Can be overridden by environment variables.
CPU_NODE_LIST=(${CPU_NODE_LIST:-"1"})
GPU_NODE_LIST=(${GPU_NODE_LIST:-"1"})

# Parameters adjustable via environment variables
CPU_NTASKS_PER_NODE=${CPU_NTASKS_PER_NODE:-2}
GPU_NTASKS_PER_NODE=${GPU_NTASKS_PER_NODE:-1}
GPUS_PER_NODE=${GPUS_PER_NODE:-4}
CPU_JOB_PREFIX=${CPU_JOB_PREFIX:-cpu}
GPU_JOB_PREFIX=${GPU_JOB_PREFIX:-gpu}

# Submit a job to SLURM with specified parameters
# Arguments: script_path nodes job_name [additional_sbatch_args...]
submit_job() {
  local script="$1"
  local nodes="$2"
  local job_name="$3"
  shift 3

  echo "Submitting ${job_name} with ${nodes} node(s)"
  sbatch --nodes="${nodes}" --job-name="${job_name}" "$@" "${script}"
}

# CPU/OpenMP sweep (batch.sh)
for n in "${CPU_NODE_LIST[@]}"; do
  submit_job "${CPU_SCRIPT}" "${n}" "${CPU_JOB_PREFIX}_n${n}" \
    --ntasks-per-node="${CPU_NTASKS_PER_NODE}"
done

# GPU/OpenACC sweep (gpu_batch.sh)
for n in "${GPU_NODE_LIST[@]}"; do
  submit_job "${GPU_SCRIPT}" "${n}" "${GPU_JOB_PREFIX}_n${n}" \
    --ntasks-per-node="${GPU_NTASKS_PER_NODE}" \
    --gres="gpu:${GPUS_PER_NODE}"
done
