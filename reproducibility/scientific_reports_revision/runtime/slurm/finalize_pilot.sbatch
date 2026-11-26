#!/usr/bin/env bash
#SBATCH --account=hbulut
#SBATCH --qos=normal
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
set -euo pipefail
source "${TSR_REVISION_ROOT:?}/slurm/configuration.sh"
source "$TSR_REVISION_ROOT/slurm/runtime_environment.sh"
: "${TSR_RUN_OUTPUT:?Set the exact new pilot output directory}" "${TSR_PILOT_JOB_ID:?Set the submitted pilot array job ID}"
export TSR_ARTIFACT_DIR="$TSR_RUN_OUTPUT/artifacts"
launcher="$TSR_REVISION_ROOT/scripts/tsr_rscript.sh"
bash slurm/collect_diagnostics.sh "$TSR_PILOT_JOB_ID"
"$launcher" scripts/check_complete.R config/tasks_pilot.csv "$TSR_RUN_OUTPUT"
"$launcher" scripts/aggregate.R config/tasks_pilot.csv "$TSR_RUN_OUTPUT"
"$launcher" scripts/pilot_acceptance.R config/tasks_pilot.csv "$TSR_RUN_OUTPUT"
python3 scripts/verify_pilot_resources.py "$TSR_RUN_OUTPUT" config/tasks_pilot.csv
printf 'PILOT_ACCEPTANCE_OK output=%s\n' "$TSR_RUN_OUTPUT"
