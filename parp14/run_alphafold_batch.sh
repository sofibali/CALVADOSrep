#!/bin/bash
set -e

source /programs/sbgrid.shrc
export TRITON_PTXAS_PATH=/usr/local/cuda-12.4/bin/ptxas
export CUDA_VISIBLE_DEVICES=0

# Configuration
DB_DIR="/mnt/alphafold3"
MODEL_DIR="/mnt/alphafold3"
OUTPUT_DIR="alphafold_outputs"
JSON_DIR="alphafold_inputs"
TOTAL_FILES=1023

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Log file
LOG_FILE="alphafold_batch_${RANDOM}.log"
echo "Starting batch predictions at $(date)" > "$LOG_FILE"
echo "Total predictions: $TOTAL_FILES" >> "$LOG_FILE"

# Counter
COMPLETED=0
FAILED=0
SKIPPED=0

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
        echo "  ✓ Completed" | tee -a "$LOG_FILE"
    else
        echo "  ✗ Failed" | tee -a "$LOG_FILE"
        FAILED=$((FAILED + 1))
    fi
done

# Summary
echo "" | tee -a "$LOG_FILE"
echo "======================================================================" | tee -a "$LOG_FILE"
echo "Batch Prediction Summary" | tee -a "$LOG_FILE"
echo "======================================================================" | tee -a "$LOG_FILE"
echo "Total predictions: $TOTAL_FILES" | tee -a "$LOG_FILE"
echo "Completed: $((COMPLETED - SKIPPED))" | tee -a "$LOG_FILE"
echo "Skipped (already exist): $SKIPPED" | tee -a "$LOG_FILE"
echo "Failed: $FAILED" | tee -a "$LOG_FILE"
echo "Finished at $(date)" | tee -a "$LOG_FILE"
echo "Log file: $LOG_FILE" | tee -a "$LOG_FILE"
