#!/bin/bash
#
# Run the complete PARP14 AlphaFold3 analysis pipeline
#
# Usage: bash run_full_pipeline.sh [background]
#   background - run with nohup in background

# Check if running in background mode
if [[ "$1" == "background" ]]; then
    TIMESTAMP=$(date +%Y%m%d_%H%M%S)
    LOGFILE="pipeline_${TIMESTAMP}.log"
    echo "Starting pipeline in background mode..."
    echo "Log file: $LOGFILE"
    echo "Monitor progress: tail -f $LOGFILE"
    echo "Check if running: ps aux | grep run_full_pipeline"
    nohup bash "$0" _run > "$LOGFILE" 2>&1 &
    echo "Pipeline started with PID: $!"
    exit 0
fi

# Internal flag for actual execution
if [[ "$1" != "_run" ]]; then
    echo "=========================================="
    echo "PARP14 AlphaFold3 Analysis Pipeline"
    echo "=========================================="
    echo ""
    echo "⚠️  Running in foreground mode"
    echo "   To run in background (survives disconnect):"
    echo "   bash run_full_pipeline.sh background"
    echo ""
    read -p "Continue in foreground? [y/N] " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

set -e  # Exit on error
set -o pipefail  # Capture errors in pipes

echo "=========================================="
echo "PARP14 AlphaFold3 Analysis Pipeline"
echo "Started: $(date)"
echo "=========================================="
echo ""

# Configuration
FASTA_FILE="../PARP14.fasta"
DOMAINS_FILE="../Domain_boundries.csv"
MAX_DOMAINS=4
NUM_SEEDS=5
AF3_MODELS="/path/to/alphafold3/models"
AF3_DATABASES="/path/to/alphafold3/databases"

# Step 1: Generate FASTA files
echo "Step 1: Generating FASTA files..."
python 01_generate_fasta.py \
    --fasta "$FASTA_FILE" \
    --domains "$DOMAINS_FILE" \
    --output fasta_files/ \
    --max-size $MAX_DOMAINS

echo ""

# Step 2: Generate AlphaFold3 inputs
echo "Step 2: Generating AlphaFold3 input files..."
python 02_generate_af3_inputs.py \
    --fasta-dir fasta_files/ \
    --output alphafold_inputs/ \
    --seed $NUM_SEEDS

echo ""

# Step 3: Run AlphaFold3 (comment out if running manually)
echo "Step 3: Running AlphaFold3 predictions..."
echo "⚠️  This step takes many hours - consider running manually"
# bash 03_run_alphafold.sh \
#     alphafold_inputs/ \
#     alphafold_outputs/ \
#     "$AF3_MODELS" \
#     "$AF3_DATABASES"

echo ""
echo "⊙ Skipping AlphaFold3 execution (uncomment in script to run)"
echo "  Run manually: bash 03_run_alphafold.sh alphafold_inputs/ alphafold_outputs/ ..."
echo ""

# Step 4: Structure inventory
echo "Step 4: Generating structure inventory..."
python 09_structure_inventory.py \
    --input ../alphafold_outputs/ \
    --output structure_inventory/

echo ""

# Step 5: Distance analysis
echo "Step 5: Running distance analysis..."
python 04_distance_analysis.py \
    --input ../alphafold_outputs/ \
    --domains "$DOMAINS_FILE" \
    --output distance_analysis/

echo ""

# Step 6: pLDDT analysis
echo "Step 6: Running pLDDT analysis..."
python 05_pLDDT_analysis.py \
    --input ../alphafold_outputs/ \
    --output-residue pLDDT_per_residue.csv

echo ""

# Step 7: SASA analysis
echo "Step 7: Running SASA analysis..."
python 06_sasa_analysis.py \
    --input ../alphafold_outputs/ \
    --output sasa_per_residue.csv

echo ""

# Summary
echo "=========================================="
echo "Pipeline Complete!"
echo "=========================================="
echo ""
echo "Output files:"
echo "  - structure_inventory/structure_summary.csv"
echo "  - structure_inventory/*.pdf (visualizations)"
echo "  - distance_analysis/domain_distances.csv"
echo "  - distance_analysis/contact_maps/*.pdf"
echo "  - distance_analysis/distance_matrices/*.npz"
echo "  - pLDDT_per_residue.csv"
echo "  - sasa_per_residue.csv"
echo ""
echo "Next steps:"
echo "  1. Open visualization notebook:"
echo "     jupyter notebook 08_visualization.ipynb"
echo "  2. Launch interactive dashboard:"
echo "     streamlit run ../interactive_dashboard.py"
echo ""
echo "=========================================="
