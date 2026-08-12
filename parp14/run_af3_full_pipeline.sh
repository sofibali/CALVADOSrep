#!/bin/bash
#
# Run AF3 FULL PIPELINE (MSA + inference) for 771 MD1L1-containing constructs
# that have no output yet. MSA runs first (CPU parallel), then inference (GPU).
#
# Usage:
#   ./run_af3_full_pipeline.sh                      # 2 GPUs, 8 MSA parallel
#   ./run_af3_full_pipeline.sh --gpu0 0 --gpu1 1    # specify GPU IDs
#   ./run_af3_full_pipeline.sh --msa-parallel 12    # control CPU parallelism
#   ./run_af3_full_pipeline.sh --skip-msa           # skip MSA (already done)
#   ./run_af3_full_pipeline.sh --msa-only           # run MSA only, no inference
#   ./run_af3_full_pipeline.sh --dry-run            # show plan without running
#
# Monitor:
#   tail -f alphafold_outputs/logs/msa_*.log        # individual MSA logs
#   tail -f run_af3_gpu0_fullpipe.log               # GPU inference logs
#   tail -f run_af3_gpu1_fullpipe.log

set +e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Defaults
GPU0=0
GPU1=1
MSA_PARALLEL=8
MSA_CPUS_PER_JOB=2
DRY_RUN=0
SKIP_MSA=0
SKIP_INFERENCE=0
JAX_CACHE_DIR="$SCRIPT_DIR/.jax_cache"
NUM_DIFFUSION_SAMPLES=5

# Paths
LIST_FILE="$SCRIPT_DIR/missing_full_pipeline.txt"
JSON_DIR="$SCRIPT_DIR/alphafold_inputs_missing"
OUTPUT_DIR="$SCRIPT_DIR/alphafold_outputs"
DB_DIR="/mnt/alphafold3"
MODEL_DIR="/mnt/alphafold3"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --gpu0)           GPU0="$2"; shift 2 ;;
        --gpu1)           GPU1="$2"; shift 2 ;;
        --msa-parallel)   MSA_PARALLEL="$2"; shift 2 ;;
        --msa-cpus)       MSA_CPUS_PER_JOB="$2"; shift 2 ;;
        --skip-msa)       SKIP_MSA=1; shift ;;
        --msa-only)       SKIP_INFERENCE=1; shift ;;
        --dry-run)        DRY_RUN=1; shift ;;
        --diffusion-samples) NUM_DIFFUSION_SAMPLES="$2"; shift 2 ;;
        *)                echo "Unknown arg: $1"; exit 1 ;;
    esac
done

# Environment
source /programs/sbgrid.shrc 2>/dev/null
export TRITON_PTXAS_PATH=/usr/local/cuda-12.4/bin/ptxas
mkdir -p "$JAX_CACHE_DIR" "$OUTPUT_DIR"

# ============================================================================
# Build list of constructs that still need work
# ============================================================================
if [ ! -f "$LIST_FILE" ]; then
    echo "ERROR: $LIST_FILE not found."
    exit 1
fi

NEED_MSA=()
NEED_INFERENCE=()
ALREADY_DONE=0

while IFS= read -r name; do
    namelc="${name,,}"
    json_file="$JSON_DIR/${name}.json"

    # AF3 may truncate names (e.g. md1l1 -> md1). Build search patterns
    # for both the original name and the AF3-truncated variant.
    namelc_alt="${namelc//md1l1/md1}"
    search_patterns=("$namelc" "$namelc_alt")

    # Skip if fully done (check original + timestamped + alt-name dirs)
    found_model=0
    for pat in "${search_patterns[@]}"; do
        outdir="$OUTPUT_DIR/$pat"
        if [ -n "$(find "$outdir" -name 'model.cif' -print -quit 2>/dev/null)" ]; then
            found_model=1; break
        fi
        for d in "$OUTPUT_DIR"/${pat}_2*/; do
            if [ -d "$d" ] && [ -n "$(find "$d" -name 'model.cif' -print -quit 2>/dev/null)" ]; then
                found_model=1; break 2
            fi
        done
    done
    if [ "$found_model" -eq 1 ]; then
        ((ALREADY_DONE++))
        continue
    fi

    # Check if MSA is done (original + timestamped + alt-name dirs)
    data_json=""
    for pat in "${search_patterns[@]}"; do
        outdir="$OUTPUT_DIR/$pat"
        data_json=$(find "$outdir" -name '*_data.json' -print -quit 2>/dev/null)
        [ -n "$data_json" ] && break
        for d in "$OUTPUT_DIR"/${pat}_2*/; do
            if [ -d "$d" ]; then
                data_json=$(find "$d" -name '*_data.json' -print -quit 2>/dev/null)
                [ -n "$data_json" ] && break 2
            fi
        done
    done

    if [ -n "$data_json" ]; then
        NEED_INFERENCE+=("$data_json")
    elif [ -f "$json_file" ]; then
        NEED_MSA+=("$json_file")
    else
        echo "WARNING: No input JSON for $name, skipping"
    fi
done < "$LIST_FILE"

N_MSA=${#NEED_MSA[@]}
N_INF=${#NEED_INFERENCE[@]}

echo "============================================================================"
echo "AF3 Full Pipeline (MD1L1 constructs)"
echo "============================================================================"
echo "  From list:         $(wc -l < "$LIST_FILE") constructs"
echo "  Already done:      $ALREADY_DONE (skipped)"
echo "  Need MSA + inf:    $N_MSA"
echo "  MSA done, need inf: $N_INF"
echo ""
echo "  Phase 1 (MSA): $MSA_PARALLEL parallel x $MSA_CPUS_PER_JOB cores"
echo "  Phase 2 (Inf): GPU 0=$GPU0, GPU 1=$GPU1"
echo "============================================================================"

if [ $DRY_RUN -eq 1 ]; then
    echo "[DRY RUN] Exiting."
    exit 0
fi

START_TIME=$(date +%s)

# ============================================================================
# PHASE 1: MSA (CPU only, parallel)
# ============================================================================
run_msa_for_json() {
    local json_file=$1
    local name=$(basename "$json_file" .json)
    local namelc="${name,,}"
    local namelc_alt="${namelc//md1l1/md1}"

    # Skip if already done (check both original and AF3-truncated name)
    for pat in "$namelc" "$namelc_alt"; do
        local out_dir="$OUTPUT_DIR/$pat"
        if [ -n "$(find "$out_dir" -name '*_data.json' -print -quit 2>/dev/null)" ]; then
            return 0
        fi
        if [ -n "$(find "$out_dir" -name 'model.cif' -print -quit 2>/dev/null)" ]; then
            return 0
        fi
        for d in "$OUTPUT_DIR"/${pat}_2*/; do
            if [ -d "$d" ]; then
                if [ -n "$(find "$d" -name '*_data.json' -o -name 'model.cif' -print -quit 2>/dev/null)" ]; then
                    return 0
                fi
            fi
        done
    done

    export OMP_NUM_THREADS=$MSA_CPUS_PER_JOB
    export OPENMM_CPU_THREADS=$MSA_CPUS_PER_JOB
    export MKL_NUM_THREADS=$MSA_CPUS_PER_JOB
    export NUMEXPR_MAX_THREADS=$MSA_CPUS_PER_JOB
    export TF_NUM_INTEROP_THREADS=$MSA_CPUS_PER_JOB
    export TF_NUM_INTRAOP_THREADS=$MSA_CPUS_PER_JOB
    export XLA_FLAGS="--xla_cpu_multi_thread_eigen=false"

    /programs/x86_64-linux/system/sbgrid_bin/run_alphafold.py \
        --db_dir "$DB_DIR" \
        --model_dir "$MODEL_DIR" \
        --output_dir "$OUTPUT_DIR" \
        --json_path "$json_file" \
        --norun_inference \
        --nhmmer_n_cpu "$MSA_CPUS_PER_JOB" \
        --jackhmmer_n_cpu "$MSA_CPUS_PER_JOB" \
        &> "$OUTPUT_DIR/logs/msa_${name,,}.log"
}

if [ $SKIP_MSA -eq 0 ] && [ $N_MSA -gt 0 ]; then
    echo ""
    echo "============================================================"
    echo "PHASE 1: MSA ($N_MSA constructs, $MSA_PARALLEL parallel)"
    echo "Started: $(date)"
    echo "============================================================"

    mkdir -p "$OUTPUT_DIR/logs"
    export -f run_msa_for_json
    export OUTPUT_DIR DB_DIR MODEL_DIR MSA_CPUS_PER_JOB

    ACTIVE=0
    SUBMITTED=0
    COMPLETED=0

    for json_file in "${NEED_MSA[@]}"; do
        run_msa_for_json "$json_file" &
        ((ACTIVE++))
        ((SUBMITTED++))

        if [ $ACTIVE -ge $MSA_PARALLEL ]; then
            wait -n 2>/dev/null || wait
            ((ACTIVE--))
            ((COMPLETED++))
        fi

        if [ $((SUBMITTED % 50)) -eq 0 ]; then
            echo "[$(date '+%H:%M:%S')] Phase 1: $SUBMITTED/$N_MSA submitted"
        fi
    done

    wait
    echo "[$(date '+%H:%M:%S')] Phase 1 complete: $N_MSA MSA jobs finished"

    # After MSA, collect _data.json files for inference (check alt names too)
    for json_file in "${NEED_MSA[@]}"; do
        name=$(basename "$json_file" .json)
        namelc="${name,,}"
        namelc_alt="${namelc//md1l1/md1}"
        data_json=""
        for pat in "$namelc" "$namelc_alt"; do
            data_json=$(find "$OUTPUT_DIR/$pat" -name '*_data.json' -print -quit 2>/dev/null)
            [ -n "$data_json" ] && break
            for d in "$OUTPUT_DIR"/${pat}_2*/; do
                [ -d "$d" ] && data_json=$(find "$d" -name '*_data.json' -print -quit 2>/dev/null)
                [ -n "$data_json" ] && break 2
            done
        done
        if [ -n "$data_json" ]; then
            NEED_INFERENCE+=("$data_json")
        else
            echo "WARNING: MSA failed for $name (no _data.json produced)"
        fi
    done
elif [ $SKIP_MSA -eq 1 ]; then
    echo "PHASE 1: SKIPPED (--skip-msa)"
    # Collect any _data.json that appeared since list was made
    for json_file in "${NEED_MSA[@]}"; do
        name=$(basename "$json_file" .json)
        namelc="${name,,}"
        namelc_alt="${namelc//md1l1/md1}"
        data_json=""
        for pat in "$namelc" "$namelc_alt"; do
            data_json=$(find "$OUTPUT_DIR/$pat" -name '*_data.json' -print -quit 2>/dev/null)
            [ -n "$data_json" ] && break
            for d in "$OUTPUT_DIR"/${pat}_2*/; do
                [ -d "$d" ] && data_json=$(find "$d" -name '*_data.json' -print -quit 2>/dev/null)
                [ -n "$data_json" ] && break 2
            done
        done
        if [ -n "$data_json" ]; then
            NEED_INFERENCE+=("$data_json")
        fi
    done
else
    echo "PHASE 1: No MSA needed (all $N_INF already have _data.json)"
fi

# ============================================================================
# PHASE 2: Inference (GPU, model loads once via --input_dir)
# ============================================================================
N_INF_TOTAL=${#NEED_INFERENCE[@]}

if [ $SKIP_INFERENCE -eq 0 ] && [ $N_INF_TOTAL -gt 0 ]; then
    echo ""
    echo "============================================================"
    echo "PHASE 2: Inference ($N_INF_TOTAL constructs, 2 GPUs)"
    echo "Started: $(date)"
    echo "============================================================"

    HALF=$(( (N_INF_TOTAL + 1) / 2 ))
    GPU0_INPUTS=("${NEED_INFERENCE[@]:0:$HALF}")
    GPU1_INPUTS=("${NEED_INFERENCE[@]:$HALF}")

    prepare_batch() {
        local batch_dir=$1
        shift
        local inputs=("$@")
        rm -rf "$batch_dir"
        mkdir -p "$batch_dir"
        for data_json in "${inputs[@]}"; do
            ln -sf "$data_json" "$batch_dir/"
        done
        echo "${#inputs[@]} _data.json files linked"
    }

    run_gpu_inference() {
        local gpu_id=$1
        local label=$2
        local batch_dir=$3
        local log_file="$SCRIPT_DIR/run_af3_gpu${label}_fullpipe.log"

        local n_inputs=$(ls "$batch_dir"/*.json 2>/dev/null | wc -l)
        if [ "$n_inputs" -eq 0 ]; then
            echo "GPU $label: nothing to process" | tee "$log_file"
            return 0
        fi

        echo "[$(date '+%H:%M:%S')] GPU $label ($gpu_id): $n_inputs inputs (inference only)" | tee "$log_file"

        CUDA_VISIBLE_DEVICES=$gpu_id \
        /programs/x86_64-linux/system/sbgrid_bin/run_alphafold.py \
            --db_dir "$DB_DIR" \
            --model_dir "$MODEL_DIR" \
            --output_dir "$OUTPUT_DIR" \
            --input_dir "$batch_dir" \
            --jax_compilation_cache_dir "$JAX_CACHE_DIR" \
            --flash_attention_implementation triton \
            --num_diffusion_samples $NUM_DIFFUSION_SAMPLES \
            --norun_data_pipeline \
            >> "$log_file" 2>&1

        echo "[$(date '+%H:%M:%S')] GPU $label: done (exit $?)" | tee -a "$log_file"
    }

    BATCH0="$SCRIPT_DIR/.af3_batch_gpu0_full"
    BATCH1="$SCRIPT_DIR/.af3_batch_gpu1_full"

    echo ""
    echo "Preparing GPU batches..."
    echo "  GPU 0: $(prepare_batch "$BATCH0" "${GPU0_INPUTS[@]}")"
    echo "  GPU 1: $(prepare_batch "$BATCH1" "${GPU1_INPUTS[@]}")"
    echo ""

    echo "Launching inference on both GPUs..."
    echo "  Monitor: tail -f $SCRIPT_DIR/run_af3_gpu0_fullpipe.log"
    echo "           tail -f $SCRIPT_DIR/run_af3_gpu1_fullpipe.log"
    echo ""

    run_gpu_inference "$GPU0" "0" "$BATCH0" &
    PID0=$!
    echo "  GPU 0: PID $PID0"

    run_gpu_inference "$GPU1" "1" "$BATCH1" &
    PID1=$!
    echo "  GPU 1: PID $PID1"

    echo ""
    echo "Waiting for both GPUs..."
    FAILURES=0
    wait $PID0 || ((FAILURES++))
    wait $PID1 || ((FAILURES++))

    rm -rf "$BATCH0" "$BATCH1"
elif [ $SKIP_INFERENCE -eq 1 ]; then
    echo "PHASE 2: SKIPPED (--msa-only)"
fi

# ============================================================================
# Summary
# ============================================================================
END_TIME=$(date +%s)
ELAPSED=$(( END_TIME - START_TIME ))

echo ""
echo "============================================================================"
echo "Full Pipeline Batch Complete"
echo "============================================================================"
echo "  Elapsed: $(( ELAPSED / 3600 ))h $(( (ELAPSED % 3600) / 60 ))m"

N_COMPLETE=$(find "$OUTPUT_DIR" -maxdepth 3 -name 'model.cif' -print 2>/dev/null | \
    sed 's|/seed-[0-9]*_sample-[0-9]*/model.cif||' | sort -u | wc -l)
echo "  Completed combinations: $N_COMPLETE / 2047 target"

if [ "${FAILURES:-0}" -gt 0 ]; then
    echo "  WARNING: $FAILURES GPU(s) had errors. Check fullpipe logs."
fi
echo "============================================================================"
