#!/usr/bin/env bash
set -euo pipefail
mode=${1:?pilot/tier1/tier2};retry_file=${2:-}
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
source "$script_dir/configuration.sh"
: "${TSR_REVISION_ROOT:?}" "${TSR_SOURCE_ROOT:?}" "${TSR_OUTPUT_ROOT:?}" "${TSR_PARTITION:?}" "${TSR_CPU_COUNT:?}" "${TSR_MEMORY:?}" "${TSR_WALLTIME:?}"
[[ ${TSR_ENV_VERIFIED:-0} == 1 ]] || { echo 'Inspect preflight and allocation first; TSR_ENV_VERIFIED is unset' >&2; exit 2; }
[[ $mode == pilot || ${TSR_FULL_APPROVED:-} == REVIEWED_AND_APPROVED ]] || { echo 'Full jobs require implementation/pilot review approval' >&2; exit 2; }
concurrency=${TSR_ARRAY_CONCURRENCY:-2}
[[ $concurrency =~ ^[0-9]+$ && $concurrency -ge 1 && $concurrency -le 20 ]] || { echo 'Concurrency must be 1..20' >&2;exit 2; }
[[ $TSR_CPU_COUNT =~ ^[0-9]+$ && $TSR_CPU_COUNT -ge 1 ]] || exit 2
cd "$TSR_REVISION_ROOT"
python3 scripts/check_transfer.py --verify-only
export TSR_MANIFEST="$TSR_REVISION_ROOT/config/tasks_${mode}.csv"
python3 scripts/dependency_acceptance_receipt.py
if [[ $mode == pilot ]];then
 : "${TSR_PILOT_RUN_LABEL:?Set a NEW run label after dependency acceptance (e.g. acceptance_JOBID)}"
 [[ $TSR_PILOT_RUN_LABEL =~ ^[A-Za-z0-9_-]+$ ]] || exit 2
 export TSR_RUN_OUTPUT="$TSR_OUTPUT_ROOT/pilot_${TSR_PILOT_RUN_LABEL}"
else
 export TSR_RUN_OUTPUT="$TSR_OUTPUT_ROOT/$mode"
fi
export TSR_TASK_LIST=${retry_file:-}
mkdir -p "$TSR_RUN_OUTPUT" "$TSR_OUTPUT_ROOT/logs"
# Counting CSV records does not load/fill models on a login node.
if [[ -n $retry_file ]];then count=$(wc -l < "$retry_file");else count=$(( $(wc -l < "$TSR_MANIFEST") - 1 ));fi
[[ $count -gt 0 ]] || { echo 'No tasks to submit';exit 0; }
groups=$(( (count + TSR_CPU_COUNT - 1) / TSR_CPU_COUNT ))
args=(--job-name="tsr_${mode}" --account="$TSR_ACCOUNT" --qos="$TSR_QOS" --partition="$TSR_PARTITION" --ntasks="$TSR_CPU_COUNT" --cpus-per-task=1 --nodes=1 --mem="$TSR_MEMORY" --time="$TSR_WALLTIME" --array="1-${groups}%${concurrency}" --chdir="$TSR_REVISION_ROOT" --output="$TSR_OUTPUT_ROOT/logs/${mode}_%A_%a.out" --error="$TSR_OUTPUT_ROOT/logs/${mode}_%A_%a.err" --export=ALL)
sbatch "${args[@]}" slurm/array.sbatch
