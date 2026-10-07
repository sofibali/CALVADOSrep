#!/bin/bash
# 4 Go conditions x 5 reps, one replicate per GPU across 0-3.
set -u
cd "$(dirname "${BASH_SOURCE[0]}")"
export PYTHON_EXE=/home/sbali/miniconda3/envs/calvados/bin/python
GPUS=(0 1 2 3); NGPU=4
for tag in p9_dtx3l_got_r15_k50 p9_dtx3l_got_r15_k100 p9_dtx3l_got_r12_k50 p9_dtx3l_got_r12_k100; do
  echo "=== $tag START $(date -Is) ==="
  ( cd binding_gotight/$tag
    i=0
    for d in $(ls -d rep-* | sort -V); do echo "${GPUS[$((i % NGPU))]} $d"; i=$((i+1)); done | \
    xargs -P "$NGPU" -L1 bash -c '
      gpu="$0"; d="$1"; cd "$d" || exit 1
      compgen -G "*.dcd" > /dev/null && exit 0
      CUDA_VISIBLE_DEVICES="$gpu" "$PYTHON_EXE" run.py > run.log 2>&1 || echo "  $d: FAILED"
      echo "  $d: done"'
  ) 2>&1 | sed "s/^/[$tag] /"
  echo "=== $tag DONE $(date -Is) ($(find binding_gotight/$tag -name '*.dcd'|wc -l)/5) ==="
done
echo "=== ALL GO-TIGHT COMPLETE $(date -Is) ==="
