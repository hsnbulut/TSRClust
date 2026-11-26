#!/usr/bin/env bash
# PREPARED ONLY: the caller must first transfer and run the hash/config acceptance check.
set -euo pipefail
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
source "$script_dir/configuration.sh"
cd "$TSR_REVISION_ROOT"
python3 scripts/check_transfer.py --verify-only
mkdir -p "$TSR_OUTPUT_ROOT/logs"
args=(--job-name=tsr_deps_accept --account="$TSR_ACCOUNT" --qos="$TSR_QOS" --partition="$TSR_PARTITION" --nodes=1 --ntasks=1 --cpus-per-task="$TSR_INSTALL_CPU_COUNT" --mem="$TSR_INSTALL_MEMORY" --time="$TSR_INSTALL_WALLTIME" --chdir="$TSR_REVISION_ROOT" --output="$TSR_OUTPUT_ROOT/logs/install_accept_%j.out" --error="$TSR_OUTPUT_ROOT/logs/install_accept_%j.err" --export="ALL,TSR_ACCEPTANCE_VERIFY_ONLY=${TSR_ACCEPTANCE_VERIFY_ONLY:-0}")
if [[ ${1:-} == --print-only ]];then printf '%q ' sbatch "${args[@]}" slurm/install_and_verify.sbatch;printf '\n';exit 0;fi
sbatch "${args[@]}" slurm/install_and_verify.sbatch
