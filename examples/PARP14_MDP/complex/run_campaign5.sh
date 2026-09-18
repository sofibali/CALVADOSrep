#!/bin/bash
# Binding campaign, 5-replicate screen.
#
# WHY 5 AND NOT 20 -- measured on the first three completed p9_dtx3l replicates
# (500 ns each, analysed at 0.5 ns):
#   * 31-37 association/dissociation events PER REPLICATE. The 40 nm box put the
#     system in the frequent-encounter regime it was designed for, so a single
#     replicate is already self-averaging -- this is not a rare-event problem.
#   * bound fraction 0.069 +/- 0.013 (SD across reps)
#       n= 5 -> 95% CI +/-0.012
#       n=20 -> 95% CI +/-0.006
#     The extra precision from 20 reps changes no conclusion.
#   * longest single contact in 1500 ns was 6 ns; ZERO episodes over 10 ns.
#     Nothing stays bound, so there is no "still bound at the end" replicate to
#     extend. Extension is triggered on episode lifetime instead -- see
#     followup_tier.sh.
#
# 12 sets x 5 reps = 60 replicates, ~27 GPU-h, ~7 h on 4 GPUs (was ~27 h).
# Sets already holding >=NREPS trajectories are skipped and marked complete.
set -u
cd "$(dirname "${BASH_SOURCE[0]}")"
export PYTHON_EXE=/home/sbali/miniconda3/envs/calvados/bin/python
GPUS=(0 1 2 3)
NREPS=${NREPS:-5}
NGPU=${#GPUS[@]}

SETS="p9_dtx3l p9_dtx3l_docked dtx3l_homo dtx3l_homo_docked p9_homo \
p14_dtx3l p14_dtx3l_docked p14_p9 p14_p9_docked ternary ternary_docked p14_homo"

mkdir -p logs

# Don't double-book GPUs still finishing replicates from a previous launcher.
# Match on the replicate's working directory rather than on a command-line
# pattern: `pgrep -f "python run.py"` also matches any monitoring shell whose
# own command line happens to contain that string, which deadlocks this wait.
sims_running() {
    for pid in $(pgrep -x python 2>/dev/null); do
        case "$(readlink /proc/$pid/cwd 2>/dev/null)" in
            */complex/binding/*) return 0 ;;
        esac
    done
    return 1
}
while sims_running; do sleep 60; done
echo "=== GPUs clear, starting 5-replicate screen $(date -Is) ==="

for s in $SETS; do
    sys=${s%_docked}
    have=$(find "binding/$s" -name '*.dcd' 2>/dev/null | wc -l)
    if [ "$have" -ge "$NREPS" ]; then
        echo "=== $s : already has $have trajectories, skipping ==="
        touch "logs/${s}.complete"; continue
    fi
    echo "=== $s : START $(date -Is) ==="
    ( cd "binding/$s"
      i=0
      for d in $(ls -d rep-* | sort -V | head -n "$NREPS"); do
          echo "${GPUS[$((i % NGPU))]} $d"; i=$((i + 1))
      done | xargs -P "$NGPU" -L1 bash -c '
          gpu="$0"; d="$1"
          cd "$d" || exit 1
          if compgen -G "*.dcd" > /dev/null; then echo "  $d: has trajectory, skipping"; exit 0; fi
          echo "  starting $d on GPU $gpu"
          CUDA_VISIBLE_DEVICES="$gpu" "$PYTHON_EXE" run.py > run.log 2>&1 || echo "  $d: FAILED"
          echo "  $d: done"
      '
    ) 2>&1 | sed "s/^/[$s] /"
    echo "=== $s : DONE $(date -Is)  ($(find binding/$s -name '*.dcd' | wc -l) trajectories) ==="
    touch "logs/${s}.complete"
done
echo "=== CAMPAIGN COMPLETE $(date -Is) ==="
