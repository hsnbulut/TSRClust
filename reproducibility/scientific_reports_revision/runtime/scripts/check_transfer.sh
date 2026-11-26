#!/usr/bin/env bash
set -euo pipefail
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
source "$script_dir/../config/truba.env"
cd "$TSR_REVISION_ROOT"
python3 scripts/check_transfer.py "$@"
