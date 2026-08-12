#!/bin/bash
# Extend fl_optimized from checkpoint (domain + custom restraints retained).
# Each run appends to the DCD until this round's target. Arg: "parallel" or "serial".
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
MODE="${1:-serial}"
PIDS=()
PY=/home/sbali/miniconda3/envs/CALVADOS/bin/python
for d in "$HERE"/seed-*_sample-*/; do
  [ -f "$d/restart.chk" ] || { echo "SKIP $d (no checkpoint)"; continue; }
  echo "EXTEND $d"
  if [ "$MODE" = parallel ]; then ( cd "$d" && $PY run.py ) & PIDS+=($!); else ( cd "$d" && $PY run.py ); fi
done
[ "$MODE" = parallel ] && wait
echo "fl_optimized extension done."
