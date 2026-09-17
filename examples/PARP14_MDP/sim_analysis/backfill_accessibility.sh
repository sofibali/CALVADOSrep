#!/usr/bin/env bash
# Backfill solid-angle accessibility (SAA / cone / shell-density) for every
# PARP14 simulation set that has trajectories on disk.
#
# Run sets SEQUENTIALLY, one invocation each: analyze_accessibility.py merges
# into data/accessibility_stats.npz, so this is incremental and crash-safe
# (a set already done stays done), but parallel invocations would race on the
# read-modify-write of that file.
#
# Sampling: --target-frames picks the stride PER TRAJECTORY so every replicate
# contributes a comparable number of frames. This matters because the sets span
# two very different populations:
#   fl / fl_optimized / noart / core / norrm  -> 3-4.2k frames, 25 replicates
#                                                from distinct AF3 seed x sample
#                                                starts (structural diversity)
#   *_full / *_go extensions                  -> 100k frames, 5 replicates
#                                                (conformational time)
# A single fixed stride cannot serve both: --stride 50 gives the 2 us runs a
# sensible 1 frame/ns but leaves a 3514-frame run with only ~69 frames.
# Domain reorientation decorrelates on tens of ns, so 2000 frames/replicate is
# still oversampled -- analyzing every frame costs ~50x more for no extra
# independent sampling (verified: replicate spread +/-0.001 SAA on the long runs).
set -u
PY=/home/sbali/miniconda3/envs/CALVADOS/bin/python
TARGET_FRAMES=${TARGET_FRAMES:-2000}
# Cap BLAS threads: the ray-casting matmuls are small, and unbounded threading
# spent ~50 core-minutes of sys time per 47 s of wall clock in contention.
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
export OPENBLAS_NUM_THREADS=$OMP_NUM_THREADS
export MKL_NUM_THREADS=$OMP_NUM_THREADS

SETS="fl fl_optimized md_full md2_md3 md3_wwe_full md3_art_full md2_wwe_full \
md2_art_full mka_wwe_full mka_full core_wwe_full_go core_full_go kh1_wwe_full \
kh1_art_full fl_wwe_full_go noart core norrm"

echo "=== accessibility backfill: $(date) | target-frames=$TARGET_FRAMES | threads=$OMP_NUM_THREADS ==="
for s in $SETS; do
  echo ""
  echo "############ $s ############"
  start=$(date +%s)
  $PY -u analyze_accessibility.py --set "$s" --target-frames "$TARGET_FRAMES" 2>&1 \
      | grep -E "SAA=|Analyzed|replicates done|Merging|Saved: accessibility_stats|ERROR|Traceback|WARN"
  echo "---- $s took $(( $(date +%s) - start ))s ----"
done
echo ""
echo "=== backfill complete: $(date) ==="
