#!/bin/bash
# Launch the 5 fl_go 2-us states. Arg: "parallel" or "serial" (default).
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
MODE="${1:-serial}"
PY=/home/sbali/miniconda3/envs/CALVADOS/bin/python
for d in "$HERE"/state-*/; do
  [ -f "$d/run.py" ] || continue
  echo "RUN $d"
  if [ "$MODE" = parallel ]; then ( cd "$d" && $PY run.py ) & else ( cd "$d" && $PY run.py ); fi
done
[ "$MODE" = parallel ] && wait
echo "fl_go done."
