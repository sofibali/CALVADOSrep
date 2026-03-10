#!/bin/bash
#
# Submit AlphaFold3 predictions for all JSON input files
# 
# This script runs AlphaFold3 predictions using SBGrid's run_alphafold.py
# Configure the directories below before running.
#
# Usage:
#   bash 03_run_alphafold.sh
#   or with nohup:
#   nohup bash 03_run_alphafold.sh > alphafold_run.log 2>&1 &
#
# Configuration:
#   - Edit DB_DIR, MODEL_DIR, OUTPUT_DIR, JSON_DIR below
#   - Adjust TOTAL_FILES to match number of JSON inputs
#   - Set CUDA_VISIBLE_DEVICES for GPU selection
#
# Features:
#   - Skips already-completed predictions
#   - Logs all output to alphafold_batch_*.log
#   - Provides progress summary at completion

set -e

LOG_FILE="alphafold_batch_$(date +%Y%m%d_%H%M%S).log"

echo "===================================================================" | tee -a "$LOG_FILE"
echo "AlphaFold3 Batch Processing Script" | tee -a "$LOG_FILE"
echo "Started at: $(date)" | tee -a "$LOG_FILE"
echo "Log file: $LOG_FILE" | tee -a "$LOG_FILE"
echo "===================================================================" | tee -a "$LOG_FILE"

echo "Sourcing SBGrid environment..." | tee -a "$LOG_FILE"
source /programs/sbgrid.shrc 2>&1 | grep -v "^[*]" | tee -a "$LOG_FILE" || true
export TRITON_PTXAS_PATH=/usr/local/cuda-12.4/bin/ptxas
export CUDA_VISIBLE_DEVICES=0

# Configuration - EDIT THESE PATHS
DB_DIR="/mnt/alphafold3"
MODEL_DIR="/mnt/alphafold3"
OUTPUT_DIR="/home/sbali/parp14/alphafold_outputs_missing"
JSON_DIR="/home/sbali/parp14/pipeline/alphafold_inputs_missing"
TOTAL_FILES=1024  # Update this to match your number of JSON files

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "" | tee -a "$LOG_FILE"
echo "Configuration:" | tee -a "$LOG_FILE"
echo "  Database Dir: $DB_DIR" | tee -a "$LOG_FILE"
echo "  Model Dir: $MODEL_DIR" | tee -a "$LOG_FILE"
echo "  Output Dir: $OUTPUT_DIR" | tee -a "$LOG_FILE"
echo "  JSON Input Dir: $JSON_DIR" | tee -a "$LOG_FILE"
echo "  Total Files: $TOTAL_FILES" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Counter
COMPLETED=0
FAILED=0
SKIPPED=0

echo "Processing JSON files..." | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Process all JSON files
for json_file in "$JSON_DIR"/*.json; do
    if [ ! -f "$json_file" ]; then
        continue
    fi
    
    COMPLETED=$((COMPLETED + 1))
    BASENAME=$(basename "$json_file" .json)
    OUTPUT_SUBDIR="$OUTPUT_DIR/${BASENAME,,}"
    
    # Check if already exists (skip if output exists)
    if [ -d "$OUTPUT_SUBDIR" ] && [ "$(ls -A "$OUTPUT_SUBDIR" 2>/dev/null)" ]; then
        echo "[$COMPLETED/$TOTAL_FILES] SKIPPED: $BASENAME (output exists)" | tee -a "$LOG_FILE"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi
    
    echo "[$COMPLETED/$TOTAL_FILES] Processing: $BASENAME" | tee -a "$LOG_FILE"
    
    # Run prediction
    if /programs/x86_64-linux/system/sbgrid_bin/run_alphafold.py \
        --db_dir "$DB_DIR" \
        --model_dir "$MODEL_DIR" \
        --output_dir "$OUTPUT_DIR" \
        --json_path "$json_file" >> "$LOG_FILE" 2>&1; then
        echo "  ✓ Completed: $BASENAME" | tee -a "$LOG_FILE"
    else
        echo "  ✗ Failed: $BASENAME" | tee -a "$LOG_FILE"
        FAILED=$((FAILED + 1))
    fi
done

# Summary
echo "" | tee -a "$LOG_FILE"
echo "=======================================================================" | tee -a "$LOG_FILE"
echo "Batch Prediction Summary" | tee -a "$LOG_FILE"
echo "=======================================================================" | tee -a "$LOG_FILE"
echo "Total predictions: $TOTAL_FILES" | tee -a "$LOG_FILE"
echo "Completed: $((COMPLETED - SKIPPED - FAILED))" | tee -a "$LOG_FILE"
echo "Skipped (already exist): $SKIPPED" | tee -a "$LOG_FILE"
echo "Failed: $FAILED" | tee -a "$LOG_FILE"
echo "Finished at $(date)" | tee -a "$LOG_FILE"
echo "Log file: $LOG_FILE" | tee -a "$LOG_FILE"
echo "=======================================================================" | tee -a "$LOG_FILE"
