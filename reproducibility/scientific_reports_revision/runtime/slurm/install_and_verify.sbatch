#!/usr/bin/env bash
#SBATCH --job-name=tsr_deps_accept
#SBATCH --account=hbulut
#SBATCH --qos=normal
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
set -euo pipefail
source "${TSR_REVISION_ROOT:?Source config/truba.env first}/slurm/configuration.sh"
cd "$TSR_REVISION_ROOT"
export TSR_ACCEPTANCE_RUN_DIR="$TSR_OUTPUT_ROOT/acceptance/job_${SLURM_JOB_ID}"
mkdir -p "$TSR_ACCEPTANCE_RUN_DIR"
acceptance_stage=initialization
acceptance_command='initialize acceptance'
finish_acceptance() {
 rc=$?
 trap - EXIT
 if [[ $rc != 0 ]];then
  python3 "$TSR_REVISION_ROOT/scripts/dependency_acceptance_receipt.py" \
   --record-failure "$TSR_ACCEPTANCE_RUN_DIR" --stage "$acceptance_stage" \
   --exit-code "$rc" --command "$acceptance_command" || echo 'FAILURE_RECEIPT_WRITE_FAILED' >&2
  printf 'DEPENDENCY_ACCEPTANCE_FAILED job=%s stage=%s exit=%s\n' "$SLURM_JOB_ID" "$acceptance_stage" "$rc" >&2
 fi
 exit "$rc"
}
trap finish_acceptance EXIT
run_acceptance_stage() {
 acceptance_stage=$1;shift
 printf -v acceptance_command '%q ' "$@"
 "$@"
}
run_acceptance_stage begin python3 scripts/dependency_acceptance_receipt.py --begin "$TSR_ACCEPTANCE_RUN_DIR"
acceptance_stage=runtime_environment
acceptance_command='source slurm/runtime_environment.sh'
source "$TSR_REVISION_ROOT/slurm/runtime_environment.sh"
[[ $(uname -s) == Linux ]] || { echo 'This job requires Linux' >&2;exit 2; }
launcher="$TSR_REVISION_ROOT/scripts/tsr_rscript.sh"
# Always verify first. A completed stamp alone is not evidence of a working library.
acceptance_stage=dependency_precheck
acceptance_command='scripts/tsr_rscript.sh scripts/verify_dependencies.R'
if run_acceptance_stage verify_dependencies "$launcher" scripts/verify_dependencies.R "$TSR_R_LIB" > "$TSR_ACCEPTANCE_RUN_DIR/dependency_precheck.log" 2>&1;then
 echo 'DEPENDENCIES_ALREADY_VALID: installation skipped'
 cat "$TSR_ACCEPTANCE_RUN_DIR/dependency_precheck.log"
else
 cat "$TSR_ACCEPTANCE_RUN_DIR/dependency_precheck.log" >&2
 if [[ ${TSR_ACCEPTANCE_VERIFY_ONLY:-0} == 1 || -d $TSR_R_LIB ]];then
  echo 'DEPENDENCY_PRECHECK_FAILED: existing packages preserved; no automatic reinstall' >&2
  exit 2
 fi
 # Only a genuinely absent library is eligible for a first installation.
 run_acceptance_stage first_install "$launcher" --install scripts/install_dependencies.R "$TSR_R_LIB"
 run_acceptance_stage verify_dependencies "$launcher" scripts/verify_dependencies.R "$TSR_R_LIB"
fi
run_acceptance_stage library_child_process "$launcher" tests/library_child_process.R
run_acceptance_stage startup_acceptance_contract python3 tests/test_startup_acceptance_regression.py
run_acceptance_stage startup_isolation python3 tests/test_r_startup_isolation.py
run_acceptance_stage minimal_competitors "$launcher" tests/minimal_competitors.R
run_acceptance_stage pre_truba_checks "$launcher" tests/pre_truba_checks.R
run_acceptance_stage run_tests "$launcher" tests/run_tests.R
run_acceptance_stage workflow_checks "$launcher" tests/workflow_checks.R
run_acceptance_stage scientific_preservation python3 tests/test_scientific_preservation.py
run_acceptance_stage success_receipt python3 scripts/dependency_acceptance_receipt.py --write "$TSR_ACCEPTANCE_RUN_DIR"
printf 'DEPENDENCY_ACCEPTANCE_OK job=%s output=%s\n' "$SLURM_JOB_ID" "$TSR_ACCEPTANCE_RUN_DIR"
# No simulation array is submitted by this job.
