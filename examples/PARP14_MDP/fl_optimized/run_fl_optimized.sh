#!/bin/bash
# Launch all 25 fl_optimized replicates in parallel (or use GNU parallel)
set -e

CWD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_EXE="/home/sbali/miniconda3/envs/CALVADOS/bin/python"
PARALLEL=${1:-4}  # number of parallel jobs (default 4)

cd "$CWD"

# Collect replicate dirs
DIRS=$(ls -d seed-*_sample-* | sort)

# Run with xargs -P for parallelism
echo "$DIRS" | xargs -I {} -P $PARALLEL bash -c '
    cd "{}" || exit 1
    if [ -f "*.dcd" ] 2>/dev/null; then
        echo "{}: already has trajectory, skipping"
        exit 0
    fi
    echo "  Starting {}..."
    '"$PYTHON_EXE"' run.py > run.log 2>&1 || echo "  {}: FAILED"
    echo "  {}: done"
'

echo "All 25 replicates finished"
