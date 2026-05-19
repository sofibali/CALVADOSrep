#!/bin/bash
# Launch ALL prepared PARP14 fragment simulations in parallel.
# 240 cores → can run ~50-60 sims simultaneously (each uses 4 threads).
#
# Usage:
#   bash run_all_fragments.sh              # default: 60 parallel jobs
#   bash run_all_fragments.sh 80           # 80 parallel jobs
#   bash run_all_fragments.sh 60 fast      # only run fragments <=300 res
#   bash run_all_fragments.sh 60 skip-done # skip already-completed
#
# Each fragment has 25 replicates (5 seeds x 5 samples). Total prepared:
# 48 fragments × 25 reps = 1200 simulations.

set -uo pipefail

CWD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_EXE="/home/sbali/miniconda3/envs/CALVADOS/bin/python"
LOG_DIR="$CWD/_run_logs"
mkdir -p "$LOG_DIR"

PARALLEL=${1:-60}
MODE=${2:-all}

echo "============================================================"
echo "Launching PARP14 fragment simulations"
echo "  Parallel jobs:  $PARALLEL"
echo "  Mode:           $MODE  (all | fast | skip-done)"
echo "  Log dir:        $LOG_DIR"
echo "============================================================"

# Build list of all sim dirs (one per replicate, across all fragments)
SIM_DIRS=()
for frag_dir in "$CWD"/*/; do
    [ -d "$frag_dir" ] || continue
    frag_name="$(basename "$frag_dir")"
    [ -f "$frag_dir/metadata.json" ] || continue

    # Fast mode: skip fragments > 300 residues
    if [ "$MODE" = "fast" ]; then
        nres=$(grep -oP '"n_residues":\s*\K[0-9]+' "$frag_dir/metadata.json")
        if [ -n "$nres" ] && [ "$nres" -gt 300 ]; then
            continue
        fi
    fi

    for rep_dir in "$frag_dir"seed-*_sample-*/; do
        [ -d "$rep_dir" ] || continue
        rep_dir="${rep_dir%/}"

        # Skip if --skip-done and DCD already exists
        if [ "$MODE" = "skip-done" ]; then
            if ls "$rep_dir"/*.dcd 2>/dev/null | grep -q '.'; then
                continue
            fi
        fi

        SIM_DIRS+=("$rep_dir")
    done
done

TOTAL=${#SIM_DIRS[@]}
echo "  Sim dirs to run: $TOTAL"

if [ "$TOTAL" -eq 0 ]; then
    echo "Nothing to run."
    exit 0
fi

run_one() {
    local sim_dir="$1"
    local frag_name=$(basename "$(dirname "$sim_dir")")
    local rep_name=$(basename "$sim_dir")
    local log_file="$LOG_DIR/${frag_name}__${rep_name}.log"

    cd "$sim_dir" || return 1
    "$PYTHON_EXE" run.py > "$log_file" 2>&1
    local rc=$?
    if [ $rc -eq 0 ]; then
        echo "  [OK]   $frag_name/$rep_name"
    else
        echo "  [FAIL] $frag_name/$rep_name  (see $log_file)"
    fi
    return $rc
}

export -f run_one
export PYTHON_EXE LOG_DIR

# Use xargs for parallel execution
printf '%s\n' "${SIM_DIRS[@]}" | xargs -I{} -P "$PARALLEL" bash -c 'run_one "$@"' _ {}

echo ""
echo "============================================================"
echo "All fragment runs complete"
echo "============================================================"
