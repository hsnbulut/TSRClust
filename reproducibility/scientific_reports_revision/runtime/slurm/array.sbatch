#!/usr/bin/env bash
#SBATCH --account=hbulut
#SBATCH --qos=normal
#SBATCH --partition=debug
set -euo pipefail
source "${TSR_REVISION_ROOT:?Source config/truba.env first}/slurm/configuration.sh"
source "$TSR_REVISION_ROOT/slurm/runtime_environment.sh"
if [[ -n ${TSR_R_LIB:-} ]];then export R_LIBS_USER="$TSR_R_LIB";fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 RCPP_PARALLEL_NUM_THREADS=1
if [[ -n ${TSR_TASK_LIST:-} ]];then mapfile -t ids < "$TSR_TASK_LIST";else
 count=$(( $(wc -l < "$TSR_MANIFEST") - 1 ));mapfile -t ids < <(seq 1 "$count")
fi
first=$(( (SLURM_ARRAY_TASK_ID-1)*TSR_CPU_COUNT ));pids=()
for ((w=0;w<TSR_CPU_COUNT;w++));do
 idx=$((first+w));[[ $idx -lt ${#ids[@]} ]] || break
 srun --exclusive --export=ALL --job-name="task_${ids[$idx]}" -N1 -n1 -c1 python3 scripts/run_task_measured.py "$TSR_MANIFEST" "${ids[$idx]}" "$TSR_RUN_OUTPUT" &
 pids+=("$!")
done
status=0;for pid in "${pids[@]}";do wait "$pid" || status=1;done
exit "$status"
