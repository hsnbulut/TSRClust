#!/usr/bin/env bash
set -euo pipefail
bash "${TSR_REVISION_ROOT:?}/slurm/submit_common.sh" pilot "${1:-}"
