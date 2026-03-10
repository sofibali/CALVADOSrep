#!/bin/bash
#
# run_composition_analysis.sh - Run full domain composition analysis pipeline
#
# This script runs all analysis needed to understand how domain composition
# affects MD1, MD2, MD3, and ART active sites.
#
# Prerequisites:
#   - AlphaFold3 structures in ../alphafold_outputs/
#   - SASA analysis in ../sasa_per_residue.csv
#   - pLDDT analysis in ./pLDDT_per_residue.csv
#   - fpocket installed (conda install -c conda-forge fpocket)
#
# Usage:
#   bash run_composition_analysis.sh
#

set -e  # Exit on error

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASE_DIR="$(dirname "$SCRIPT_DIR")"

echo "=============================================="
echo "PARP14 Domain Composition Analysis Pipeline"
echo "=============================================="
echo ""
echo "Base directory: $BASE_DIR"
echo "Script directory: $SCRIPT_DIR"
echo ""

# Create output directories
echo "Creating output directories..."
mkdir -p "$BASE_DIR/results/per_structure"
mkdir -p "$BASE_DIR/results/per_residue"
mkdir -p "$BASE_DIR/results/pocket_analysis"
mkdir -p "$BASE_DIR/figures/domain_composition_effects"
mkdir -p "$BASE_DIR/figures/active_site_accessibility"
mkdir -p "$BASE_DIR/figures/interdomain_contacts"
mkdir -p "$BASE_DIR/logs"

# Timestamp for logs
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="$BASE_DIR/logs/composition_analysis_${TIMESTAMP}.log"

echo "Log file: $LOG_FILE"
echo ""

# Check for required input files
echo "Checking input files..."
if [ ! -d "$BASE_DIR/alphafold_outputs" ]; then
    echo "ERROR: alphafold_outputs directory not found!"
    exit 1
fi

if [ ! -f "$BASE_DIR/sasa_per_residue.csv" ]; then
    echo "WARNING: sasa_per_residue.csv not found. Some analyses will be skipped."
fi

if [ ! -f "$SCRIPT_DIR/pLDDT_per_residue.csv" ]; then
    echo "WARNING: pLDDT_per_residue.csv not found. Some analyses will be skipped."
fi

# Check for fpocket
if ! command -v fpocket &> /dev/null; then
    echo "WARNING: fpocket not found. Install with: conda install -c conda-forge fpocket"
    echo "         Pocket analysis will be skipped."
    SKIP_FPOCKET=true
else
    SKIP_FPOCKET=false
fi

echo ""
echo "=============================================="
echo "Step 1: Inter-Domain Distance Analysis"
echo "=============================================="
echo ""

python "$SCRIPT_DIR/09_interdomain_dynamics.py" \
    --input "$BASE_DIR/alphafold_outputs" \
    --output "$BASE_DIR/results/per_structure/interdomain_distances.csv" \
    --stats-output "$BASE_DIR/results/per_structure/interdomain_statistics.csv" \
    --summary-output "$BASE_DIR/results/per_structure/key_domain_summary.csv" \
    --workers 8 \
    2>&1 | tee -a "$LOG_FILE"

echo ""
echo "=============================================="
echo "Step 2: Pocket Detection (fpocket)"
echo "=============================================="
echo ""

if [ "$SKIP_FPOCKET" = true ]; then
    echo "Skipping fpocket analysis (not installed)"
else
    python "$SCRIPT_DIR/08_fpocket_analysis.py" \
        --input "$BASE_DIR/alphafold_outputs" \
        --output "$BASE_DIR/results/pocket_analysis" \
        --csv-output "$BASE_DIR/results/per_structure/fpocket_results.csv" \
        --workers 8 \
        --seed "seed-1_sample-0" \
        2>&1 | tee -a "$LOG_FILE"
fi

echo ""
echo "=============================================="
echo "Step 3: Active Site Analysis"
echo "=============================================="
echo ""

python "$SCRIPT_DIR/10_active_site_analysis.py" \
    --sasa "$BASE_DIR/sasa_per_residue.csv" \
    --plddt "$SCRIPT_DIR/pLDDT_per_residue.csv" \
    --fpocket "$BASE_DIR/results/per_structure/fpocket_results.csv" \
    --output "$BASE_DIR/results/per_structure/active_site_metrics.csv" \
    --summary-output "$BASE_DIR/results/per_structure/active_site_summary.csv" \
    2>&1 | tee -a "$LOG_FILE"

echo ""
echo "=============================================="
echo "Step 4: Generate Publication Figures"
echo "=============================================="
echo ""

python "$SCRIPT_DIR/11_composition_effect_figures.py" \
    --input "$BASE_DIR/results/per_structure" \
    --output "$BASE_DIR/figures/domain_composition_effects" \
    2>&1 | tee -a "$LOG_FILE"

echo ""
echo "=============================================="
echo "Analysis Complete!"
echo "=============================================="
echo ""
echo "Results saved to:"
echo "  - $BASE_DIR/results/per_structure/"
echo "  - $BASE_DIR/figures/domain_composition_effects/"
echo ""
echo "Key output files:"
echo "  - interdomain_distances.csv      Inter-domain distance metrics"
echo "  - key_domain_summary.csv         MD1/MD2/MD3/ART summary by composition"
echo "  - fpocket_results.csv            Pocket detection results"
echo "  - active_site_metrics.csv        Combined active site analysis"
echo "  - active_site_summary.csv        Summary by composition category"
echo ""
echo "Figures generated:"
echo "  - fig1_composition_effect_heatmap.pdf"
echo "  - fig2_accessibility_boxplots.pdf"
echo "  - fig3_interdomain_distances.pdf"
echo "  - fig4_pocket_volume_distribution.pdf"
echo "  - fig5_seed_variability.pdf"
echo "  - fig6_composition_summary.pdf"
echo "  - fig7_contact_network.pdf"
echo ""
echo "Log file: $LOG_FILE"
