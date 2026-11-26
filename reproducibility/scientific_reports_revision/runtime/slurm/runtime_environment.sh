#!/usr/bin/env bash
# Must only be sourced after SLURM has allocated a compute job.
[[ -n ${SLURM_JOB_ID:-} ]] || { echo 'Refusing runtime/package execution outside a SLURM compute job' >&2; return 2; }
[[ ${SLURM_JOB_ACCOUNT:-hbulut} == hbulut ]] || { echo 'Wrong SLURM account' >&2; return 2; }
[[ ${SLURM_JOB_QOS:-normal} == normal ]] || { echo 'Wrong SLURM QoS' >&2; return 2; }
node_list=${SLURM_JOB_NODELIST:-${SLURM_NODELIST:-}}
[[ -n $node_list ]] || { echo 'SLURM allocation node list is missing' >&2;return 2; }
current_node=$(hostname -s)
node_found=0
while IFS= read -r allocated_node;do
 [[ ${allocated_node%%.*} == "$current_node" ]] && node_found=1
done < <(scontrol show hostnames "$node_list")
[[ $node_found == 1 ]] || { echo 'Current host is not an allocated compute node' >&2;return 2; }
cd "$TSR_REVISION_ROOT"
python3 scripts/check_transfer.py --verify-only
# Remove inherited Mac/user overrides without changing any user configuration file.
unset R_HOME R_ARCH R_LIBS R_LIBS_USER R_LIBS_SITE R_MAKEVARS_USER R_MAKEVARS_SITE
module purge
module load "$TSR_GCC_MODULE"
module load "$TSR_R_MODULE"
export R_MAKEVARS_USER=/dev/null
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 RCPP_PARALLEL_NUM_THREADS=1 MAKEFLAGS=-j1
source "$TSR_REVISION_ROOT/slurm/r_startup_environment.sh"
export SLURM_EXPORT_ENV=ALL
bash scripts/tsr_rscript.sh --probe
module list 2>&1
