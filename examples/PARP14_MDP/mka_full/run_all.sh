#!/bin/bash
# Launch all prepared mka_full replicates. Arg: "parallel" or "serial" (default serial).
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
MODE="${1:-serial}"
PY=/home/sbali/miniconda3/envs/CALVADOS/bin/python
for d in "$HERE"/seed-*_sample-*/; do
  [ -f "$d/run.py" ] || continue
  if [ "$MODE" = parallel ]; then ( cd "$d" && $PY run.py ) & else ( cd "$d" && $PY run.py ); fi
done
[ "$MODE" = parallel ] && wait
echo "mka_full done."
