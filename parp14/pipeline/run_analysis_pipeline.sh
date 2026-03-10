#!/bin/bash
#
# Run only the analysis portion of PARP14 pipeline (skip input generation and AlphaFold3)
# This script analyzes already-completed AlphaFold3 structures
#
# Usage: bash run_analysis_pipeline.sh [background]
#   background - run with nohup in background (survives SSH disconnect)

# Check if running in background mode
if [[ "$1" == "background" ]]; then
    TIMESTAMP=$(date +%Y%m%d_%H%M%S)
    LOGFILE="analysis_pipeline_${TIMESTAMP}.log"
    echo "Starting analysis pipeline in background mode..."
    echo "Log file: $LOGFILE"
    echo "Monitor progress: tail -f $LOGFILE"
    echo "Check if running: ps aux | grep run_analysis_pipeline"
    nohup bash "$0" _run > "$LOGFILE" 2>&1 &
    PID=$!
    echo "Pipeline started with PID: $PID"
    echo ""
    echo "Useful commands:"
    echo "  tail -f $LOGFILE              # Watch progress"
    echo "  tail -100 $LOGFILE            # See last 100 lines"
    echo "  ps -p $PID -o pid,etime,cmd   # Check if still running"
    echo "  kill $PID                     # Stop the pipeline"
    exit 0
fi

# Internal flag for actual execution
if [[ "$1" != "_run" ]]; then
    echo "=========================================="
    echo "PARP14 Analysis Pipeline"
    echo "=========================================="
    echo ""
    echo "⚠️  Running in foreground mode"
    echo "   To run in background (survives disconnect):"
    echo "   bash run_analysis_pipeline.sh background"
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
echo "PARP14 Analysis Pipeline"
echo "Started: $(date)"
echo "=========================================="
echo ""

# Configuration
DOMAINS_FILE="../Domain_boundries.csv"
INPUT_DIR="../alphafold_outputs"
OUTPUT_BASE="."

# Verify input directory exists
if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: Input directory not found: $INPUT_DIR"
    exit 1
fi

# Count structures
STRUCTURE_COUNT=$(find "$INPUT_DIR" -maxdepth 1 -type d | wc -l)
echo "Found $((STRUCTURE_COUNT - 1)) structure directories in $INPUT_DIR"
echo ""

# Step 1: Structure inventory
echo "=========================================="
echo "Step 1: Structure Inventory"
echo "=========================================="
python 09_structure_inventory.py \
    --input "$INPUT_DIR" \
    --output structure_inventory/

if [[ -f structure_inventory/structure_summary.csv ]]; then
    echo ""
    echo "Structure summary:"
    head -20 structure_inventory/structure_summary.csv | column -t -s,
fi

echo ""

# Step 2: Distance analysis
echo "=========================================="
echo "Step 2: Distance Analysis"
echo "=========================================="
python 04_distance_analysis.py \
    --input "$INPUT_DIR" \
    --domains "$DOMAINS_FILE" \
    --output distance_analysis/

if [[ -f distance_analysis/domain_distances.csv ]]; then
    NUM_MEASUREMENTS=$(tail -n +2 distance_analysis/domain_distances.csv | wc -l)
    echo ""
    echo "Generated $NUM_MEASUREMENTS distance measurements"
fi

echo ""

# Step 3: Confidence scores
echo "=========================================="
echo "Step 3: Confidence Score Analysis"
echo "=========================================="
python 04b_extract_confidence_scores.py \
    --input "$INPUT_DIR" \
    --output confidence_scores.csv

if [[ -f confidence_scores.csv ]]; then
    NUM_STRUCTURES=$(tail -n +2 confidence_scores.csv | wc -l)
    echo ""
    echo "Extracted confidence scores for $NUM_STRUCTURES structures"
fi

echo ""

# Step 4: pLDDT analysis
echo "=========================================="
echo "Step 4: pLDDT Analysis"
echo "=========================================="
python 05_pLDDT_analysis.py \
    --input "$INPUT_DIR" \
    --output-residue pLDDT_per_residue.csv

if [[ -f pLDDT_per_residue.csv ]]; then
    NUM_RESIDUES=$(tail -n +2 pLDDT_per_residue.csv | wc -l)
    echo ""
    echo "Extracted pLDDT scores for $NUM_RESIDUES residue entries"
fi

echo ""

# Step 5: SASA analysis
echo "=========================================="
echo "Step 5: SASA Analysis"
echo "=========================================="
python 06_sasa_analysis.py \
    --input "$INPUT_DIR" \
    --output sasa_per_residue.csv

if [[ -f sasa_per_residue.csv ]]; then
    NUM_RESIDUES=$(tail -n +2 sasa_per_residue.csv | wc -l)
    echo ""
    echo "Calculated SASA for $NUM_RESIDUES residue entries"
fi

echo ""

# Summary
echo "=========================================="
echo "Pipeline Complete!"
echo "Completed: $(date)"
echo "=========================================="
echo ""
echo "Output files generated:"
echo "  ✓ structure_inventory/structure_summary.csv"
echo "  ✓ structure_inventory/domain_combination_matrix.pdf"
echo "  ✓ structure_inventory/structure_counts_heatmap.pdf"
echo "  ✓ distance_analysis/domain_distances.csv"
echo "  ✓ confidence_scores.csv"
echo "  ✓ pLDDT_per_residue.csv"
echo "  ✓ sasa_per_residue.csv"
echo ""
echo "File sizes:"
ls -lh structure_inventory/structure_summary.csv \
       distance_analysis/domain_distances.csv \
       confidence_scores.csv \
       pLDDT_per_residue.csv \
       sasa_per_residue.csv 2>/dev/null | awk '{print "  " $9 ": " $5}'

# Step 6: Generate comprehensive visualizations
echo "=========================================="
echo "Step 6: Comprehensive Visualizations"
echo "=========================================="
python 10_comprehensive_visualization.py \
    --input . \
    --output visualizations/ \
    --use-builtin-sites

echo ""
echo "Next steps:"
echo "  1. Review visualizations in: visualizations/"
echo "  2. Add active site analysis (see VISUALIZATION_README.md)"
echo "  3. Review structure inventory:"
echo "     less structure_inventory/structure_summary.csv"
echo "  4. Launch interactive dashboard:"
echo "     cd .. && streamlit run interactive_dashboard.py"
echo ""
echo "=========================================="
