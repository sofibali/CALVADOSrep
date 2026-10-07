#!/bin/bash
# binding_xl: 21 reproducible crosslinks, IMProv-tiered thresholds, 5 replicates.
set -u
cd "$(dirname "${BASH_SOURCE[0]}")"
export PYTHON_EXE=/home/sbali/miniconda3/envs/calvados/bin/python
GPUS=(0 1 2 3); NGPU=4
cd binding_xl/p9_dtx3l_xl30_k20
i=0
for d in $(ls -d rep-* | sort -V); do echo "${GPUS[$((i % NGPU))]} $d"; i=$((i+1)); done | \
xargs -P "$NGPU" -L1 bash -c '
  gpu="$0"; d="$1"; cd "$d" || exit 1
  compgen -G "*.dcd" > /dev/null && { echo "  $d: has trajectory, skipping"; exit 0; }
  echo "  starting $d on GPU $gpu"
  CUDA_VISIBLE_DEVICES="$gpu" "$PYTHON_EXE" run.py > run.log 2>&1 || echo "  $d: FAILED"
  echo "  $d: done"'
echo "=== binding_xl COMPLETE $(date -Is) ==="
