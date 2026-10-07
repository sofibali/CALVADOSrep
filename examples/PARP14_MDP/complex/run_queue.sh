#!/bin/bash
# Follow-up queue for the binding campaign.
#
# 1. ternary_docked rep-6..rep-20 (15 more replicates)
#    The only set with a long-lived state: one replicate held a continuous
#    PARP14-DTX3L contact for 278 ns, against a next-longest of 13 ns anywhere
#    else. With n=1 we cannot tell a rare-but-real bound state from a
#    docked-start kinetic trap. 20 replicates total gives the recurrence rate.
#    This is a REPLICATE problem, not a frame-rate problem.
#
# 2. Crosslink-restrained runs + matched control (20 replicates)
#    REPLACES the box-size titration, which cannot do the job: at 25.9 uM a
#    nanomolar complex is >98% bound, so bound fraction SATURATES and cannot
#    distinguish HyRes's implied 1 pM from the experimental ~10 nM. A box where
#    f is sensitive to a nM Kd would need to be ~550-1200 nm. Out of reach.
#
#    Instead: restrain to the three published BS3 crosslinks (Ashok 2022,
#    doi 10.1042/BCJ20210722, FL PARP9 + FL DTX3L) at r = 2.0 nm (BS3 reach):
#      PARP9 K557-DTX3L K401 / K632-K401 / K557-K363
#      xl_h5 / xl_h20 / xl_h100   harmonic, k = 5/20/100 kJ/mol/nm^2
#      xl_go15                    Go, 15 kJ/mol well depth (reversible)
#    Control = the existing unrestrained binding/p9_dtx3l (7 reps, done).
#    Binding is an INPUT here, so bound fraction is no longer a result; the
#    question is whether the interface concentrates onto the crosslinks and what
#    else appears as a consequence -- which needs the difference map, because
#    Go-KH restraints already manufactured a spurious contact once here.
#
# ~15 GPU-hours total, ~4 h on 4 GPUs. One replicate per GPU (measured:
# packing 2-4 per GPU costs ~55% of aggregate throughput on this box).
set -u
cd "$(dirname "${BASH_SOURCE[0]}")"
export PYTHON_EXE=/home/sbali/miniconda3/envs/calvados/bin/python   # lowercase = the CUDA env
GPUS=(0 1 2 3)
NGPU=${#GPUS[@]}
mkdir -p logs

# "<tree>/<set>:<first>:<last>"
QUEUE=(
  "binding/ternary_docked:6:20"
  "binding_go/p9_dtx3l_xl_h20:1:5"
  "binding_go/p9_dtx3l_xl_go15:1:5"
  "binding_go/p9_dtx3l_xl_h5:1:5"
  "binding_go/p9_dtx3l_xl_h100:1:5"
)

sims_running() {
    for pid in $(pgrep -x python 2>/dev/null); do
        case "$(readlink /proc/$pid/cwd 2>/dev/null)" in
            */complex/binding/*|*/complex/binding_go/*) return 0 ;;
        esac
    done
    return 1
}
while sims_running; do sleep 60; done
echo "=== QUEUE START $(date -Is) ==="

for entry in "${QUEUE[@]}"; do
    path="${entry%%:*}"; rest="${entry#*:}"; first="${rest%%:*}"; last="${rest##*:}"
    tag=$(basename "$path")
    [ -d "$path" ] || { echo "=== $tag : MISSING, skipped ==="; continue; }
    echo "=== $tag : rep-$first..rep-$last  START $(date -Is) ==="
    ( cd "$path"
      i=0
      for n in $(seq "$first" "$last"); do
          [ -d "rep-$n" ] || continue
          echo "${GPUS[$((i % NGPU))]} rep-$n"; i=$((i + 1))
      done | xargs -P "$NGPU" -L1 bash -c '
          gpu="$0"; d="$1"
          cd "$d" || exit 1
          if compgen -G "*.dcd" > /dev/null; then echo "  $d: has trajectory, skipping"; exit 0; fi
          echo "  starting $d on GPU $gpu"
          CUDA_VISIBLE_DEVICES="$gpu" "$PYTHON_EXE" run.py > run.log 2>&1 || echo "  $d: FAILED"
          echo "  $d: done"
      '
    ) 2>&1 | sed "s/^/[$tag] /"
    echo "=== $tag : DONE $(date -Is)  ($(find "$path" -name '*.dcd' | wc -l) trajectories) ==="
    touch "logs/${tag}.complete"
done
echo "=== QUEUE COMPLETE $(date -Is) ==="
