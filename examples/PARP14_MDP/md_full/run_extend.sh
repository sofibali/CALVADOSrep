#!/bin/bash
# Extend md_full from checkpoint (+EXTRA_NS set by setup_extend_md_full.py).
# Each run appends to the DCD. Arg: "parallel" (all at once) or "serial" (default).
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
echo "md_full extension done."
