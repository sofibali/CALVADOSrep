#!/bin/bash
# Run MD3-ART (595 res) construct: 25 replicates x 20 ns
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

echo "=================================================="
echo "MD3-ART (595 res, 25 replicates x 20 ns)"
echo "=================================================="
for seed in 1 2 3 4 5; do
    for sample in 0 1 2 3 4; do
        run_sim "md3art/seed-${seed}_sample-${sample}"
    done
done

if [ "$MODE" = "parallel" ] && [ ${#PIDS[@]} -gt 0 ]; then
    echo ""
    echo "Waiting for ${#PIDS[@]} parallel jobs..."
    for pid in "${PIDS[@]}"; do
        wait "$pid"
    done
fi

echo ""
echo "=================================================="
echo "MD3-ART COMPLETE"
echo "=================================================="
