#!/usr/bin/env bash
# The single production R launcher. No module loading, authentication or submission.
set -euo pipefail
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
export TSR_REVISION_ROOT=$(cd "$script_dir/.." && pwd -P)
cd "$TSR_REVISION_ROOT"
: "${TSR_R_LIB:?Source config/truba.env, or set the explicit local validation library}"
source slurm/r_startup_environment.sh
export TSR_RSCRIPT="$TSR_REVISION_ROOT/scripts/tsr_rscript.sh"
install_mode=0
probe_mode=0
if [[ ${1:-} == --probe ]];then probe_mode=1;shift;fi
if [[ ${1:-} == --install ]];then
  install_mode=1;shift
  [[ ${1:-} == scripts/install_dependencies.R ]] || { echo 'Install mode is limited to the dependency installer' >&2;exit 2; }
fi
run_r() {
 if [[ $(uname -s) == Darwin ]];then
  export R_HOME=${TSR_LOCAL_R_HOME:?Set TSR_LOCAL_R_HOME to the R 4.5.2 installation}
  if [[ ${1:-} == -e ]];then
   "$R_HOME/bin/exec/R" --vanilla --slave "$@"
  else
   local target=${1:?R script};shift
   "$R_HOME/bin/exec/R" --vanilla --slave --file="$target" --args "$@"
  fi
 else
  Rscript --vanilla "$@"
 fi
}
if [[ $probe_mode == 1 ]];then
 run_r -e 'stopifnot(as.character(getRversion())=="4.5.2");cat(R.version.string,"\nR_HOME=",R.home(),"\nR_LIBS_USER=",Sys.getenv("R_LIBS_USER"),"\n",sep="")'
 exit 0
fi
if [[ $install_mode == 1 ]];then
 # Target may be absent on first installation. Installer itself enforces compute host.
 unset R_LIBS R_LIBS_USER
 run_r -e 'stopifnot(as.character(getRversion())=="4.5.2")'
else
 run_r -e 'source("R/library_environment.R"); initialize_revision_library(required=TRUE); assert_revision_environment();'
fi
run_r "$@"
