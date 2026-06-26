#!/bin/bash
# Extend PARP14 CALVADOS simulations from checkpoint.
# Round 2: md/core/mka 10 ns -> 20 ns (+10 ns). FL already at 20 ns (skipped).
#
# Usage:
#   bash run_extend.sh              # all sets, sequential
#   bash run_extend.sh parallel     # all sets, parallel
#   bash run_extend.sh fl           # FL only
#   bash run_extend.sh md           # md only
#   bash run_extend.sh core         # core only
#   bash run_extend.sh mka          # mka only

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
    if [ ! -f "$dir/restart.chk" ]; then
        echo "SKIP: $dir (no checkpoint)"
        return
    fi
    echo "EXTEND: $dir"

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
    echo "Extending: $label"
    echo "=================================================="
    for seed in 1 2 3 4 5; do
        for sample in 0 1 2 3 4; do
            run_sim "${prefix}/seed-${seed}_sample-${sample}"
        done
    done
}

case "$MODE" in
    fl)       run_set "fl" "Full-Length (-> 50 ns)" ;;
    md)       run_set "md" "Macrodomains (-> 25 ns)" ;;
    core)     run_set "core" "Core (-> 25 ns)" ;;
    mka)      run_set "mka" "MKA (-> 25 ns)" ;;
    parallel)
        run_set "md" "MD (-> 20 ns)"
        run_set "core" "CORE (-> 20 ns)"
        run_set "mka" "MKA (-> 20 ns)"
        ;;
    *)
        run_set "md" "MD (-> 20 ns)"
        run_set "core" "CORE (-> 20 ns)"
        run_set "mka" "MKA (-> 20 ns)"
        ;;
esac

if [ "$MODE" = "parallel" ] && [ ${#PIDS[@]} -gt 0 ]; then
    echo ""
    echo "Waiting for ${#PIDS[@]} parallel jobs..."
    for pid in "${PIDS[@]}"; do
        wait "$pid"
    done
fi

echo ""
echo "=================================================="
echo "EXTENSION COMPLETE"
echo "=================================================="
