#!/usr/bin/env bash
set -euo pipefail

# Collects key metrics from logs and dumps them into a text file.
# Usage: ./collect_results.sh [log_directory] [output.txt]

# Default log directory and output file
LOG_DIR="${1:-logs}"
OUT_FILE="${2:-sweep_summary.txt}"

# Check if log directory exists
if [[ ! -d "${LOG_DIR}" ]]; then
  echo "Log directory does not exist: ${LOG_DIR}" >&2
  exit 1
fi

# Enable nullglob to handle case when no files match
shopt -s nullglob
LOGS=("${LOG_DIR}"/*.out "${LOG_DIR}"/*.log)
if [[ ${#LOGS[@]} -eq 0 ]]; then
  echo "No .out/.log files found in ${LOG_DIR}" >&2
  exit 1
fi

# Extract first matching field from log file, return "NA" if not found
# Arguments: pattern field_number file_path
first_field_or_na() {
  local pattern="$1"
  local field="$2"
  local file="$3"
  local val
  val=$(grep -m1 "${pattern}" "${file}" | awk -v f="${field}" '{print $f}')
  [[ -z "${val}" ]] && val="NA"
  echo "${val}"
}

# Generate summary report with metrics from all log files
{
  echo "Run Summary (generated $(date))"
  echo -e "file\tnodes\ttasks_per_node\tgpus\tdelta_mass\tdelta_energy\tcpu_time_s"
  for f in "${LOGS[@]}"; do
    # Extract configuration and performance metrics from each log file
    nodes=$(first_field_or_na "^Nodes" 3 "${f}")
    tasks=$(first_field_or_na "^Tasks / Node" 5 "${f}")
    gpus=$(first_field_or_na "Number of GPUs available" 5 "${f}")
    mass=$(first_field_or_na "Fractional Delta Mass" 4 "${f}")
    energy=$(first_field_or_na "Fractional Delta Energy" 4 "${f}")
    time_s=$(first_field_or_na "USED CPU TIME" 4 "${f}")
    echo -e "$(basename "${f}")\t${nodes}\t${tasks}\t${gpus}\t${mass}\t${energy}\t${time_s}"
  done
} > "${OUT_FILE}"

echo "Summary written to ${OUT_FILE}"
