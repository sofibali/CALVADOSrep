# PARP14 AlphaFold3 Domain Analysis Pipeline (Phase 1)

Structure generation phase: AlphaFold3 predictions for PARP14 domain combinations.

**Phase 2** (CALVADOS coarse-grained MD simulations) is in [`../examples/PARP14_MDP/`](../examples/PARP14_MDP/). See [`../docs/PARP14_PROJECT.md`](../docs/PARP14_PROJECT.md) for the full two-phase protocol.

---

A comprehensive pipeline for analyzing the structural components of PARP14 (Protein mono-ADP-ribosyltransferase PARP14) using AlphaFold3 predictions. This pipeline systematically generates domain combinations, runs AlphaFold3 structure predictions, and performs detailed structural analysis.

## Overview

PARP14 is a multi-domain protein containing 17 distinct domains:
- **RNA Recognition Motifs (RRM)**: RRM1, RRM2, RRM3
- **K Homology domains (KH)**: KH1-6, KH7a, KHb, KH8
- **Macrodomains (MD)**: MD1, MD2, MD3
- **WWE domain**: WWE
- **ADP-ribosyltransferase domain**: ART

This pipeline generates and analyzes all biologically relevant domain combinations while respecting sequential constraints.

## Features

- **Domain Combination Generation**: Automatically generates all sequential domain combinations with biological constraints
- **AlphaFold3 Integration**: Creates JSON inputs and runs batch predictions
- **Structural Analysis**: Comprehensive analysis tools for:
  - Contact maps and inter-domain interfaces
  - Distance matrices and measurements
  - Secondary structure analysis
  - Structural quality metrics (pLDDT, PAE)
  - Domain-domain interaction analysis

## Pipeline Workflow

```
PARP14.fasta + Domain_boundries.csv
    ↓
1. Generate Domain Combinations (generate_domain_pairs.py)
    ↓
domain_combinations/*.fasta
    ↓
2. Create AlphaFold3 Inputs (generate_alphafold_inputs.py)
    ↓
alphafold_inputs/*.json
    ↓
3. Run AlphaFold3 Predictions (run_alphafold_batch.py / run_alphafold_batch_parallel.sh)
    ↓
alphafold_outputs/*/
    ↓
4. Structural Analysis (analyze_structures.py, contact_analysis.py, etc.)
    ↓
Analysis Results
```

## Requirements

### Software Dependencies
- Python 3.8+
- AlphaFold3 (with access to databases and models)
- CUDA-capable GPU (recommended)

### Python Packages
```bash
pip install -r requirements.txt
```

Required packages:
- numpy
- pandas
- biopython
- matplotlib
- seaborn
- scipy

## Installation

1. Clone this repository:
```bash
git clone https://github.com/sofibali/CALVADOSrep.git
cd parp14-alphafold-analysis
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

3. Set up AlphaFold3:
   - Ensure AlphaFold3 is installed and `run_alphafold.py` is in your PATH
   - Set the database and model directories in the configuration

## Usage

### Step 1: Generate Domain Combinations

Generate all FASTA files for domain combinations:

```bash
python generate_domain_pairs.py
```

This creates FASTA files in `domain_combinations/` directory. The script respects:
- Sequential domain order
- KH1-6 as a functional group
- KH7a-KHb-KH8 constraint (must appear together)

### Step 2: Create AlphaFold3 Input Files

Generate JSON input files for AlphaFold3:

```bash
python generate_alphafold_inputs.py
```

Creates JSON files in `alphafold_inputs/` directory.

### Step 3: Run AlphaFold3 Predictions

#### Option A: Sequential Batch Processing
```bash
python run_alphafold_batch.py
```

#### Option B: Parallel Processing (Recommended for large batches)
```bash
bash run_alphafold_batch_parallel.sh
```

Predictions will be saved to `alphafold_outputs/`.

### Step 4: Analyze Structures

Run comprehensive structural analysis:

```bash
# Generate contact maps
python analyze_contacts.py

# Analyze domain interfaces
python analyze_interfaces.py

# Calculate structural metrics
python analyze_structures.py

# Generate summary report
python generate_analysis_report.py
```

## Key Files

### Input Files
- `PARP14.fasta` - Full PARP14 protein sequence
- `Domain_boundries.csv` - Domain boundary definitions

### Pipeline Scripts
- `generate_domain_pairs.py` - Generates domain combination FASTA files
- `generate_alphafold_inputs.py` - Creates AlphaFold3 JSON inputs
- `run_alphafold_batch.py` - Batch AlphaFold3 runner
- `run_alphafold_batch_parallel.sh` - Parallel batch processing
- `count_grouped_combinations.py` - Domain combination statistics

### Analysis Scripts (in `analysis/` directory)
- `analyze_contacts.py` - Contact map generation
- `analyze_interfaces.py` - Inter-domain interface analysis
- `analyze_structures.py` - Structure quality and metrics
- `extract_plddt.py` - pLDDT score extraction
- `calculate_distances.py` - Distance matrix calculations
- `generate_analysis_report.py` - Summary report generation

## Domain Combination Rules

1. **Sequential Constraint**: Domains must maintain their native sequence order
2. **KH1-6 Group**: KH1, KH2, KH3, KH4, KH5, KH6 are treated as a functional unit
3. **KH7a-KHb-KH8 Constraint**: If KH7a is included, KHb and KH8 must also be included

These constraints are biologically motivated to model realistic domain architectures.

## Output Structure

```
parp14-alphafold-analysis/
├── domain_combinations/      # Generated FASTA files
├── alphafold_inputs/          # AlphaFold3 JSON inputs
├── alphafold_outputs/         # AlphaFold3 predictions
│   └── <domain_combo>/
│       ├── fold_*.cif         # Structure files
│       ├── summary_confidences_*.json
│       └── ranking_scores.json
├── analysis_results/          # Analysis outputs
│   ├── contact_maps/
│   ├── distance_matrices/
│   ├── interface_analysis/
│   └── reports/
└── jax_cache/                 # JAX compilation cache
```

## Configuration

Edit the configuration section in each script to customize:

```python
# AlphaFold3 settings
DB_DIR = '/mnt/alphafold3'         # Database directory
MODEL_DIR = '/mnt/alphafold3'      # Model directory
GPU_ID = 0                          # GPU device ID

# Pipeline settings
OUTPUT_DIR = 'alphafold_outputs'
JSON_DIR = 'alphafold_inputs'
FASTA_DIR = 'domain_combinations'
```

## Performance Notes

- **Domain Combinations**: ~127 sequential combinations with constraints
- **Prediction Time**: Varies by domain size (1-10 minutes per prediction)
- **Parallel Processing**: Recommended for batches >50 predictions
- **Memory**: ~8GB GPU memory for typical domain combinations
- **Cache**: JAX compilation cache significantly speeds up subsequent runs

## Troubleshooting

### AlphaFold3 Issues
- Ensure `CUDA_VISIBLE_DEVICES` is set correctly
- Check database and model paths are accessible
- Verify JSON input format matches AlphaFold3 requirements

### Memory Issues
- Use parallel processing with GPU distribution
- Process smaller batches sequentially
- Reduce model seeds if needed

### Analysis Issues
- Ensure all required structure files are present
- Check file paths in analysis scripts
- Verify PDB/CIF file formats are correct

## Citation

If you use this pipeline, please cite:

```
PARP14 AlphaFold3 Analysis Pipeline
https://github.com/sofibali/CALVADOSrep
```

And the AlphaFold3 paper:
```
Abramson, J., et al. (2024). Accurate structure prediction of biomolecular 
interactions with AlphaFold 3. Nature.
```

## License

MIT License - see LICENSE file for details

## Related Documentation

| Document | Description |
|----------|-------------|
| [`../examples/PARP14_MDP/README.md`](../examples/PARP14_MDP/README.md) | Phase 2: CALVADOS MD simulations |
| [`../docs/PARP14_PROJECT.md`](../docs/PARP14_PROJECT.md) | Full two-phase protocol and analysis |
| [`../docs/CALVADOS_QUICKSTART.md`](../docs/CALVADOS_QUICKSTART.md) | Running CALVADOS simulations |

## Acknowledgments

- AlphaFold3 team at Google DeepMind
- CALVADOS: KULL-Centre/CALVADOS
