#!/usr/bin/env bash
# DSS crosslink feasibility (Euclidean + Xwalk-style SASD) for every PARP14 set
# with trajectories, feeding the construct-discrimination analysis.
#
# Defaults are derived, not tuned: --min-seq-sep 9 = ceil(3.0/0.38)+1, and
# --target-frames 2000 picks the stride per trajectory (verified to reproduce
# the every-frame result on fl: 135/150 either way, ~24x cheaper).
set -u
PY=/home/sbali/miniconda3/envs/CALVADOS/bin/python
TOP=${TOP:-500}   # inter-domain pairs get the budget first (see run_sasd_pass)
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}

SETS="fl fl_optimized norrm noart core md_full mka_full mka_wwe_full \
core_full_go core_wwe_full_go kh1_art_full kh1_wwe_full fl_wwe_full_go \
md2_art_full md2_wwe_full md3_art_full md3_wwe_full md1_md2 md2_md3"

echo "=== crosslink sweep: $(date) | sasd-top=$TOP ==="
for s in $SETS; do
  echo ""; echo "######## $s ########"
  start=$(date +%s)
  $PY -u analyze_lys_contacts.py --set "$s" --sasd --sasd-top "$TOP" 2>&1 \
    | grep -E "SASD|replicates completed|Saved: data|ERROR|Traceback|no trajectory"
  echo "---- $s took $(( $(date +%s) - start ))s ----"
done
echo ""; echo "=== done: $(date) ==="
