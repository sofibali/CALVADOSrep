#!/bin/bash
# Screen multiple featurization strategies side-by-side for one CALVADOS
# simulation set. Produces diagnostic plots only (ACF + 2D density of
# PCA/TICA with K-means overlay + ITS validation), no full analysis.
#
# Usage:
#   bash sweep_cluster_features.sh <set_name> [K] [tica_lag]
#   bash sweep_cluster_features.sh fl_optimized
#   bash sweep_cluster_features.sh md 5 50
#
# Each featurization writes to a distinct suffix so plots don't overwrite.

set -uo pipefail

SET=${1:?Usage: $0 SET_NAME [K] [tica_lag]}
K=${2:-5}
TICA_LAG=${3:-50}

PY=/home/sbali/miniconda3/envs/CALVADOS/bin/python
SCRIPT="$(dirname "$0")/cluster_states.py"

# Override BLAS threads to avoid segfault on many-core machines
export OPENBLAS_NUM_THREADS=8
export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8

ITS_LAGS="10 25 50 100 200 500"

echo "============================================================"
echo "Featurization screen for set: $SET  (k=$K, tica_lag=$TICA_LAG)"
echo "ITS lags: $ITS_LAGS frames"
echo "============================================================"

# Featurization configs to try
CONFIGS=(
    "com|--features com"
    "com_rg|--features com --add-rg"
    "ca_stride10|--features ca --ca-stride 10"
    "ca_stride25|--features ca --ca-stride 25"
)

for entry in "${CONFIGS[@]}"; do
    name="${entry%%|*}"
    flags="${entry##*|}"

    echo
    echo "─── $name ──────────────────────────────────────"
    $PY $SCRIPT --set "$SET" --k "$K" --reduce tica \
        --tica-lag "$TICA_LAG" --its-lags $ITS_LAGS \
        --screen $flags
    rc=$?
    if [ $rc -ne 0 ]; then
        echo "  $name: failed (rc=$rc)"
    fi
done

echo
echo "============================================================"
echo "Screen complete. Compare plots in figures/:"
echo "  cluster_acf_${SET}*.png      <- pick features with slow decay"
echo "  cluster_density_tica_${SET}*.png  <- well-separated clusters"
echo "  cluster_density_pca_${SET}*.png"
echo "  cluster_its_${SET}*.png      <- plateau = converged kinetic ts"
echo "============================================================"
echo "Then run full analysis with chosen featurization:"
echo "  $PY $SCRIPT --set $SET --k $K --reduce tica --tica-lag $TICA_LAG \\"
echo "      --features <com|ca> [--ca-stride N] [--add-rg]"
