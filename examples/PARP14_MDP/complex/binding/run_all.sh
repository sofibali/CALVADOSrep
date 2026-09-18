#!/bin/bash
# Run every prepared binding set, one after another.
#   bash run_all.sh [PER_GPU]
set -e
CWD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
for s in "$CWD"/*/run.sh; do
    echo "=== $(basename $(dirname $s)) ==="
    bash "$s" "${1:-3}"
done
