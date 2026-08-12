#!/bin/bash
# One-shot pipeline: featurize+TICA (stage 1) -> select states (stage 2) ->
# full-atom backmap + clash removal (stage 3). Edit the config block or override
# any variable on the command line, e.g.:
#
#   SET=fl_optimized FEATURES=ca K=4 REFERENCE=input/parp14.pdb \
#       FDOMAINS=input/domains.yaml bash tica_pipeline/run_pipeline.sh
#
# Defaults reproduce the md_full rigid-body-pose, 4 spread states result.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
cd "$ROOT"                       # run from examples/PARP14_MDP (so set/dcd paths resolve)
PY="${PY:-python}"              # or set PY=/path/to/envs/calvados-tica/bin/python

# ----------------------------- config -----------------------------
SET="${SET:-md_full}"
FEATURES="${FEATURES:-pose}"          # pose | pose_iface | ca | interface_ca | com | ...
CA_STRIDE="${CA_STRIDE:-25}"
TICA_LAG="${TICA_LAG:-200}"           # frames; 200 = 2 ns
K="${K:-4}"
REP_MODE="${REP_MODE:-spread}"        # spread | pcca
REFERENCE="${REFERENCE:-md_full/input/ref_allatom.pdb}"   # all-atom ref w/ sidechains
FDOMAINS="${FDOMAINS:-md_full/input/domains.yaml}"
# ------------------------------------------------------------------

TAG="${SET}_${FEATURES}"
NPZ="data/tica_${TAG}.npz"
STATES="states/${TAG}"

echo "==== Stage 1: featurize + TICA  (set=$SET feature=$FEATURES lag=${TICA_LAG}fr) ===="
$PY tica_pipeline/01_featurize_tica.py --set "$SET" --features "$FEATURES" \
    --ca-stride "$CA_STRIDE" --tica-lag "$TICA_LAG" --out "$NPZ"

echo "==== Stage 2: select $K states  (rep-mode=$REP_MODE) ===="
$PY tica_pipeline/02_select_states.py --tica-npz "$NPZ" --k "$K" \
    --rep-mode "$REP_MODE" --msm-lag "$TICA_LAG" --out-dir "$STATES"

echo "==== Stage 3: full-atom backmap + clash removal ===="
$PY tica_pipeline/03_backmap_minimize.py --states-dir "$STATES" \
    --reference "$REFERENCE" --fdomains "$FDOMAINS"

echo "==== DONE: full-atom states -> ${STATES}_allatom/ ===="
