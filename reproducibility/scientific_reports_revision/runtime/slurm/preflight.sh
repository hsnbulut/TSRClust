#!/usr/bin/env bash
set -euo pipefail
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
source "$script_dir/configuration.sh"
: "${TSR_REVISION_ROOT:?Set revision root}" "${TSR_PARTITION:?Set verified candidate partition}"
mkdir -p "$TSR_REVISION_ROOT/results/environment"
{
 date -u
 sinfo -o '%P %a %l %D %c %m'
 scontrol show partition "$TSR_PARTITION"
 sacctmgr show assoc user="$(id -un)" format=User,Account,Partition,QOS || true
 module avail R 2>&1 || true
 module list 2>&1 || true
} > "$TSR_REVISION_ROOT/results/environment/truba_environment.txt"
cat "$TSR_REVISION_ROOT/results/environment/truba_environment.txt"
# Read-only scheduler inspection only: perform package/R tests in a compute job.
