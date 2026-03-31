#!/bin/bash
# Run PARP14 CALVADOS simulations
# 7 constructs x 25 replicates x 20 ns = 3500 ns total
#
# Usage:
#   bash run_all.sh              # run all sequentially
#   bash run_all.sh parallel     # run all in parallel (background)
#   bash run_all.sh fl           # run only full-length
#   bash run_all.sh md3art       # run only MD3-ART
#   bash run_all.sh md           # run only macrodomains

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
    echo "$label (25 replicates x 20 ns)"
    echo "=================================================="
    for seed in 1 2 3 4 5; do
        for sample in 0 1 2 3 4; do
            run_sim "${prefix}/seed-${seed}_sample-${sample}"
        done
    done
}

# Determine which sets to run
case "$MODE" in
    fl)     run_set "fl"     "Full Length (1-1801)" ;;
    norrm)  run_set "norrm"  "KH1-ART (315-1801)" ;;
    noart)  run_set "noart"  "KH1-WWE (315-1602)" ;;
    core)   run_set "core"   "KH7-ART (738-1801)" ;;
    mka)    run_set "mka"    "MD1-ART (790-1801)" ;;
    md)     run_set "md"     "MD1-MD3 (790-1388)" ;;
    md3art) run_set "md3art" "MD3-ART (1207-1801)" ;;
    parallel)
        run_set "fl"     "Full Length (1-1801)"
        run_set "norrm"  "KH1-ART (315-1801)"
        run_set "noart"  "KH1-WWE (315-1602)"
        run_set "core"   "KH7-ART (738-1801)"
        run_set "mka"    "MD1-ART (790-1801)"
        run_set "md"     "MD1-MD3 (790-1388)"
        run_set "md3art" "MD3-ART (1207-1801)"
        ;;
    *)
        run_set "fl"     "Full Length (1-1801)"
        run_set "norrm"  "KH1-ART (315-1801)"
        run_set "noart"  "KH1-WWE (315-1602)"
        run_set "core"   "KH7-ART (738-1801)"
        run_set "mka"    "MD1-ART (790-1801)"
        run_set "md"     "MD1-MD3 (790-1388)"
        run_set "md3art" "MD3-ART (1207-1801)"
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
