#!/bin/bash
# Launch all prepared md_full replicates. Arg: "parallel" or "serial" (default serial).
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
MODE="${1:-serial}"
for d in "$HERE"/seed-*_sample-*/; do
  [ -f "$d/run.py" ] || continue
  if [ "$MODE" = parallel ]; then ( cd "$d" && python run.py & ); else ( cd "$d" && python run.py ); fi
done
wait
