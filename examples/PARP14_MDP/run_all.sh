#!/bin/bash
# Run all PARP14 CALVADOS simulations
# 4 sets x 25 replicates x 5 ns = 500 ns total (125 ns per set)
#
# Usage:
#   bash run_all.sh              # run all sequentially
#   bash run_all.sh parallel     # run all in parallel (background)
#   bash run_all.sh fl           # run only full-length
#   bash run_all.sh md           # run only macrodomains
#   bash run_all.sh core         # run only core construct
#   bash run_all.sh mka          # run only MD+KH+ART construct

set -e
cd "$(dirname "$0")"

MODE="${1:-sequential}"
PIDS=()

run_sim() {
    local dir="$1"
    if [ ! -f "$dir/config.yaml" ]; then
        echo "SKIP: $dir (no config.yaml)"
        return
    fi
    if [ -f "$dir/restart.chk" ]; then
        echo "RUN (restart): $dir"
    else
        echo "RUN (fresh):   $dir"
    fi

    if [ "$MODE" = "parallel" ]; then
        (cd "$dir" && python run.py) &
        PIDS+=($!)
    else
        (cd "$dir" && python run.py)
    fi
}

run_set() {
    local prefix="$1"
    local label="$2"
    echo "=================================================="
    echo "$label (25 replicates x 5 ns = 125 ns)"
    echo "=================================================="
    for seed in 1 2 3 4 5; do
        for sample in 0 1 2 3 4; do
            run_sim "${prefix}_seed-${seed}_sample-${sample}"
        done
    done
}

# Determine which sets to run
case "$MODE" in
    fl)   run_set "fl" "Full-Length PARP14 (1801 res)" ;;
    md)   run_set "md" "Macrodomains only (586 res)" ;;
    core) run_set "core" "Core construct (1051 res)" ;;
    mka)  run_set "mka" "MD+KHb+ART construct (930 res)" ;;
    *)
        run_set "fl" "Full-Length PARP14 (1801 res)"
        run_set "md" "Macrodomains only (586 res)"
        run_set "core" "Core construct (1051 res)"
        run_set "mka" "MD+KHb+ART construct (930 res)"
        ;;
esac

# Wait for parallel jobs
if [ "$MODE" = "parallel" ] && [ ${#PIDS[@]} -gt 0 ]; then
    echo ""
    echo "Waiting for ${#PIDS[@]} parallel jobs..."
    for pid in "${PIDS[@]}"; do
        wait "$pid"
    done
fi

echo ""
echo "=================================================="
echo "SIMULATIONS COMPLETE"
echo "=================================================="
echo "Protocol: 25 replicates x 5 ns = 125 ns per set"
echo "Equilibration: discard first 0.5 ns of each replicate"
echo "Effective: 25 x 4.5 ns = 112.5 ns compiled per set"
