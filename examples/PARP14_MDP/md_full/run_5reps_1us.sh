#!/bin/bash
# Extend md_full seed-{1..5}_sample-0 from checkpoint to 1us total.
# Each run appends to the DCD. Arg: "parallel" (all at once) or "serial" (default).
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
MODE="${1:-serial}"
PY=/home/sbali/miniconda3/envs/CALVADOS/bin/python
d="$HERE/seed-1_sample-0/"
[ -f "$d/restart.chk" ] || { echo "SKIP $d (no checkpoint)"; exit 1; }
echo "EXTEND $d"
if [ "$MODE" = parallel ]; then ( cd "$d" && $PY run.py ) & else ( cd "$d" && $PY run.py ); fi
d="$HERE/seed-2_sample-0/"
[ -f "$d/restart.chk" ] || { echo "SKIP $d (no checkpoint)"; exit 1; }
echo "EXTEND $d"
if [ "$MODE" = parallel ]; then ( cd "$d" && $PY run.py ) & else ( cd "$d" && $PY run.py ); fi
d="$HERE/seed-3_sample-0/"
[ -f "$d/restart.chk" ] || { echo "SKIP $d (no checkpoint)"; exit 1; }
echo "EXTEND $d"
if [ "$MODE" = parallel ]; then ( cd "$d" && $PY run.py ) & else ( cd "$d" && $PY run.py ); fi
d="$HERE/seed-4_sample-0/"
[ -f "$d/restart.chk" ] || { echo "SKIP $d (no checkpoint)"; exit 1; }
echo "EXTEND $d"
if [ "$MODE" = parallel ]; then ( cd "$d" && $PY run.py ) & else ( cd "$d" && $PY run.py ); fi
d="$HERE/seed-5_sample-0/"
[ -f "$d/restart.chk" ] || { echo "SKIP $d (no checkpoint)"; exit 1; }
echo "EXTEND $d"
if [ "$MODE" = parallel ]; then ( cd "$d" && $PY run.py ) & else ( cd "$d" && $PY run.py ); fi
[ "$MODE" = parallel ] && wait
echo "md_full 5-rep 1us extension done."
