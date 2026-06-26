#!/bin/bash
# PARP14 batch analysis launcher — runs all 4 steps sequentially with logging.
# Started by the 6h scheduled task.

set -u  # treat undefined vars as errors; DO NOT use -e so failures continue
cd /home/sbali/CALVADOS/examples/PARP14_MDP

PY=/home/sbali/miniconda3/envs/CALVADOS/bin/python
export OPENBLAS_NUM_THREADS=8
export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8
export NUMEXPR_MAX_THREADS=8

LOG_DIR=_batch_logs
mkdir -p "$LOG_DIR"
TS=$(date +%F_%H-%M-%S)
SUM="$LOG_DIR/summary_${TS}.log"
echo "=== PARP14 batch started $(date) ===" | tee "$SUM"

run_step() {
    local name="$1"; shift
    local log="$LOG_DIR/${TS}_${name}.log"
    echo "" | tee -a "$SUM"
    echo "── [$name] $(date '+%H:%M:%S') — running: $*" | tee -a "$SUM"
    "$@" > "$log" 2>&1
    local rc=$?
    if [ $rc -eq 0 ]; then
        echo "   [$name] OK ($(date '+%H:%M:%S'))" | tee -a "$SUM"
    else
        echo "   [$name] FAILED rc=$rc — see $log" | tee -a "$SUM"
    fi
    return 0  # always continue
}

# ─── STEP 1 ──────────────────────────────────────────────────────
echo "" | tee -a "$SUM"
echo "════ STEP 1 — trajectory analyses ════" | tee -a "$SUM"
run_step "01_energy"      $PY analyze_all.py --include-fragments --energy --workers 32
run_step "02_dmap_cmap"   $PY analyze_all.py --include-fragments --dmap --cmap --workers 32
run_step "03_fnc_conf"    $PY analyze_all.py --include-fragments --fnc --conf-prop --workers 32
run_step "04_access"      $PY analyze_accessibility.py --include-fragments --nrays 200

# ─── STEP 2 ──────────────────────────────────────────────────────
echo "" | tee -a "$SUM"
echo "════ STEP 2 — clustering 8 named sets ════" | tee -a "$SUM"

NAMED_SETS=(fl fl_optimized md core mka norrm noart md3art)
declare -a OK_SETS

for SET in "${NAMED_SETS[@]}"; do
    if [ ! -d "$SET" ]; then
        echo "   [cluster_$SET] SKIP — no directory" | tee -a "$SUM"
        continue
    fi
    # Has any trajectory?
    if ! ls "$SET"/seed-*_sample-*/parp14*.dcd > /dev/null 2>&1; then
        echo "   [cluster_$SET] SKIP — no .dcd files" | tee -a "$SUM"
        continue
    fi
    run_step "cluster_$SET" $PY cluster_states.py --set "$SET" \
        --features ca --ca-stride 25 \
        --reduce tica --tica-lag 80 \
        --n-microstates 200 --pcca \
        --msm-lag 80 --k 5 --ck-test \
        --no-silhouette --workers 16
    if [ $? -eq 0 ]; then
        OK_SETS+=("$SET")
    fi
done

# ─── STEP 3 ──────────────────────────────────────────────────────
echo "" | tee -a "$SUM"
echo "════ STEP 3 — clustering fragments ════" | tee -a "$SUM"

for FRAG_DIR in fragments/*/; do
    [ -d "$FRAG_DIR" ] || continue
    NAME=$(basename "$FRAG_DIR")
    [ -f "$FRAG_DIR/metadata.json" ] || continue
    # Trajectory check
    if ! ls "$FRAG_DIR"seed-*_sample-*/parp14*.dcd > /dev/null 2>&1; then
        echo "   [cluster_frag_$NAME] SKIP — no .dcd files" | tee -a "$SUM"
        continue
    fi
    run_step "cluster_frag_$NAME" $PY cluster_states.py --set "$NAME" \
        --features ca --ca-stride 25 \
        --reduce tica --tica-lag 80 \
        --n-microstates 200 --pcca \
        --msm-lag 80 --k 5 --ck-test \
        --no-silhouette --workers 16
    if [ $? -eq 0 ]; then
        OK_SETS+=("$NAME")
    fi
done

# ─── STEP 4 ──────────────────────────────────────────────────────
echo "" | tee -a "$SUM"
echo "════ STEP 4 — PyMOL sessions for ${#OK_SETS[@]} sets ════" | tee -a "$SUM"

for SET in "${OK_SETS[@]}"; do
    run_step "pse_$SET" $PY make_cluster_pse.py --auto --set "$SET" --out-pml-only
done

echo "" | tee -a "$SUM"
echo "=== PARP14 batch finished $(date) ===" | tee -a "$SUM"
echo "Summary log: $SUM"
