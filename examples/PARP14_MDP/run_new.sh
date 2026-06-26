#!/bin/bash
# Run new PARP14 constructs: norrm (1474 res) and noart (1275 res)
# 25 replicates x 20 ns each
#
# Usage:
#   bash run_new.sh              # both sets, sequential
#   bash run_new.sh parallel     # both sets, parallel
#   bash run_new.sh norrm        # No-RRM only
#   bash run_new.sh noart        # No-ART only

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
    echo "$label (25 replicates x 20 ns = 500 ns)"
    echo "=================================================="
    for seed in 1 2 3 4 5; do
        for sample in 0 1 2 3 4; do
            run_sim "${prefix}/seed-${seed}_sample-${sample}"
        done
    done
}

case "$MODE" in
    norrm)    run_set "norrm" "No-RRM (1474 res)" ;;
    noart)    run_set "noart" "No-ART (1275 res)" ;;
    parallel)
        run_set "norrm" "No-RRM (1474 res)"
        run_set "noart" "No-ART (1275 res)"
        ;;
    *)
        run_set "norrm" "No-RRM (1474 res)"
        run_set "noart" "No-ART (1275 res)"
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
echo "SIMULATIONS COMPLETE"
echo "=================================================="
echo "Protocol: 25 replicates x 20 ns = 500 ns per set"
echo "Equilibration: discard first 0.5 ns of each replicate"
