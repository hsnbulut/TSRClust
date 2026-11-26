#!/usr/bin/env bash
# Compatibility alias: all ordinary local execution uses the common launcher.
set -euo pipefail
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
exec "$script_dir/tsr_rscript.sh" "$@"
