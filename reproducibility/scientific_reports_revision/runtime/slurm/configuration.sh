#!/usr/bin/env bash
# Shared configuration; no installation, transfer, module loading or submission here.
tsr_slurm_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
source "$tsr_slurm_dir/../config/truba.env"
[[ $TSR_ACCOUNT == hbulut && $TSR_QOS == normal && $TSR_PARTITION == debug ]] || { echo 'Verified hbulut/normal/debug configuration mismatch' >&2; return 2; }
[[ $TSR_GCC_MODULE == comp/gcc/14.1.0 && $TSR_R_MODULE == apps/R/4.5.2-gcc-14.1.0 ]] || { echo 'Verified module names mismatch' >&2; return 2; }
unset tsr_slurm_dir
