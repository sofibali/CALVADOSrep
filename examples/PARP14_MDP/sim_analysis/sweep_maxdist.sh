#!/usr/bin/env bash
# Empirical sensitivity of SAA to the ray probe length (MAX_DIST).
set -u
PY=/home/sbali/miniconda3/envs/CALVADOS/bin/python
export OMP_NUM_THREADS=4
for D in 1.0 2.0 3.0 4.0 5.0 6.0 8.0 10.0; do
  echo "#### max_dist = $D nm"
  $PY -u analyze_accessibility.py --set fl --target-frames 800 \
      --max-dist "$D" --tag "_md${D}" 2>&1 | grep -E "SAA=|Analyzed"
done
