#!/usr/bin/env bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
count=$(( $(wc -l < config/tasks_pilot.csv) - 1 ))
for ((i=1;i<=count;i++));do bash scripts/tsr_rscript.sh scripts/run_task.R config/tasks_pilot.csv "$i" results/pilot;done
bash scripts/tsr_rscript.sh scripts/check_complete.R config/tasks_pilot.csv results/pilot
bash scripts/tsr_rscript.sh scripts/aggregate.R config/tasks_pilot.csv results/pilot
