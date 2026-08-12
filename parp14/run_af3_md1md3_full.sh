#!/bin/bash
# Run AlphaFold3 for the COMPLETE contiguous MD1-MD3 construct (FL 790-1388, 599 res,
# including the 1194-1206 MD2-MD3 linker). Single input, both phases (MSA + inference).
# Mirrors run_af3_missing.sh env/flags. Needs a GPU node + SBGrid AF3.
#
#   bash run_af3_md1md3_full.sh
#
# Output: alphafold_outputs/md1l1_md2_md3_full/  (+ per-seed seed-*/model.cif)
set +e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
INPUT_JSON="$SCRIPT_DIR/alphafold_inputs_md1md3_full/md1l1_md2_md3_full.json"
OUTPUT_DIR="$SCRIPT_DIR/alphafold_outputs"
DB_DIR="/mnt/alphafold3"
MODEL_DIR="/mnt/alphafold3"
MSA_CPUS="${MSA_CPUS:-8}"

source /programs/sbgrid.shrc 2>/dev/null
export TRITON_PTXAS_PATH=/usr/local/cuda-12.4/bin/ptxas
export OMP_NUM_THREADS=$MSA_CPUS OPENMM_CPU_THREADS=$MSA_CPUS MKL_NUM_THREADS=$MSA_CPUS
export NUMEXPR_MAX_THREADS=$MSA_CPUS TF_NUM_INTEROP_THREADS=$MSA_CPUS TF_NUM_INTRAOP_THREADS=$MSA_CPUS
export XLA_FLAGS="--xla_cpu_multi_thread_eigen=false"

mkdir -p "$OUTPUT_DIR/logs"
echo "AF3 for md1l1_md2_md3_full (599 res) -- both phases"
/programs/x86_64-linux/system/sbgrid_bin/run_alphafold.py \
    --db_dir "$DB_DIR" \
    --model_dir "$MODEL_DIR" \
    --output_dir "$OUTPUT_DIR" \
    --json_path "$INPUT_JSON" \
    --nhmmer_n_cpu "$MSA_CPUS" \
    --jackhmmer_n_cpu "$MSA_CPUS" \
    --flash_attention_implementation triton \
    2>&1 | tee "$OUTPUT_DIR/logs/md1md3_full.log"

echo "Done. Then: cd ../examples/PARP14_MDP && python prepare_md_full.py"
