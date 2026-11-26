#!/usr/bin/env bash
set -euo pipefail
: "${TSR_OUTPUT_ROOT:?}"
job=${1:?SLURM job ID}
[[ $job =~ ^[0-9]+(_[0-9]+)?$ ]] || { echo 'Invalid SLURM job ID' >&2;exit 2; }
mkdir -p "$TSR_OUTPUT_ROOT/diagnostics"
sacct -j "$job" --format=JobID,JobName,State,ExitCode,Elapsed,TotalCPU,MaxRSS,ReqMem,AllocCPUS --units=K -P > "$TSR_OUTPUT_ROOT/diagnostics/sacct_${job}.csv"
shopt -s nullglob
files=("$TSR_OUTPUT_ROOT/logs/"*_${job}.out "$TSR_OUTPUT_ROOT/logs/"*_${job}.err "$TSR_OUTPUT_ROOT/logs/"*_${job}_*.out "$TSR_OUTPUT_ROOT/logs/"*_${job}_*.err)
findings="$TSR_OUTPUT_ROOT/diagnostics/log_findings_${job}.txt"
if [[ ${#files[@]} == 0 ]];then
 printf 'No matching log files for job %s\n' "$job" > "$findings"
elif command -v rg >/dev/null 2>&1;then
 rg -n -i 'error|warning|failed|out.of.memory|time.limit|segfault' "${files[@]}" > "$findings" || true
else
 grep -n -E -i 'error|warning|failed|out.of.memory|time.limit|segfault' "${files[@]}" > "$findings" || true
fi
