# PARP14 Pipeline Quickstart Guide

**Phase 1** of the PARP14 Structure-to-Dynamics Pipeline

A step-by-step guide to generate, predict, and analyze PARP14 domain combinations using AlphaFold3. Structures generated here feed into [Phase 2: CALVADOS MD Simulations](/home/sbali/CALVADOS/docs/PARP14_PROJECT.md).

---

## 1. Installation

### Create conda environment

```bash
conda create -n parp14-analysis python=3.10
conda activate parp14-analysis
```

### Install dependencies

```bash
pip install numpy pandas scipy biopython matplotlib seaborn tqdm pyyaml
pip install freesasa
```

### Optional: For interactive dashboard

```bash
pip install streamlit
```

### Clone the repository

```bash
cd /home/sbali
git clone <repository-url> parp14
cd parp14
```

---

## 2. Input Files

### Required files

| File | Description |
|------|-------------|
| `PARP14.fasta` | Full protein sequence (1801 aa) |
| `Domain_boundries.csv` | Domain coordinate definitions |

### Domain_boundries.csv format

```csv
Domain,Start,End
RRM1,1,145
RRM2,146,224
RRM3,225,314
KH1-KH6,315,737
KH7a,738,789
MD1,790,981
MD1L1,790,1004
MD2,1004,1193
MD3,1207,1388
KHb-KH8,1389,1533
WWE,1534,1602
ART,1603,1801
```

---

## 3. Generate Domain Combinations

### Step 1: Generate FASTA files

```bash
cd pipeline

python 01_generate_fasta.py \
    --fasta ../PARP14.fasta \
    --domains ../Domain_boundries.csv \
    --output ../fasta_files/ \
    --max-size 4
```

This creates FASTA files for all valid domain combinations:
```
fasta_files/
├── RRM1.fasta
├── RRM1_RRM2.fasta
├── RRM1_RRM2_RRM3.fasta
├── KH1-KH6.fasta
├── KH1-KH6_KH7a_MD1.fasta
└── ... (~1000 files)
```

### Step 2: Generate AlphaFold3 input files

```bash
python 02_generate_af3_inputs.py \
    --fasta-dir ../fasta_files/ \
    --output ../alphafold_inputs/ \
    --seeds 5
```

This creates JSON files for AlphaFold3:
```
alphafold_inputs/
├── rrm1.json
├── rrm1_rrm2.json
└── ... (~1000 files)
```

---

## 4. Run AlphaFold3 Predictions

### Sequential execution

```bash
bash 03_run_alphafold.sh
```

### Parallel execution (multiple GPUs)

```bash
bash 03_run_alphafold_parallel.sh
```

### Configuration

Edit the script to set paths:

```bash
DB_DIR="/mnt/alphafold3"           # AlphaFold3 database
MODEL_DIR="/mnt/alphafold3"        # Model weights
OUTPUT_DIR="../alphafold_outputs/" # Output directory
CUDA_VISIBLE_DEVICES=0             # GPU selection
```

### Output structure

```
alphafold_outputs/
├── rrm1/
│   ├── seed-1_sample-0/
│   │   ├── model.cif
│   │   ├── summary_confidences.json
│   │   ├── confidences.json
│   │   └── ranking_scores.json
│   ├── seed-2_sample-0/
│   └── ... (5 seeds total)
└── ...
```

---

## 5. Analyze Structures

### Distance analysis

```bash
python 04_distance_analysis.py \
    --input ../alphafold_outputs/ \
    --domains ../Domain_boundries.csv \
    --output ../distance_analysis/
```

**Output files:**
- `domain_distances.csv` - Inter-domain distance metrics
- `distance_matrices/*.npz` - Full distance matrices
- `contact_maps/*.pdf` - Contact map visualizations

### Confidence score analysis

```bash
python 05_pLDDT_analysis.py \
    --input ../alphafold_outputs/ \
    --output-residue ../pLDDT_per_residue.csv \
    --output-structure ../confidence_scores.csv
```

**Output files:**
- `pLDDT_per_residue.csv` - Per-residue confidence (2 GB)
- `confidence_scores.csv` - Per-structure metrics

### SASA analysis

```bash
python 06_sasa_analysis.py \
    --input ../alphafold_outputs/ \
    --output ../sasa_per_residue.csv \
    --domains ../Domain_boundries.csv
```

**Output:** `sasa_per_residue.csv` - Per-residue solvent accessibility (586 MB)

---

## 6. Generate CALVADOS Inputs

Prepare structures for coarse-grained MD simulations:

```bash
python 07_generate_calvados_inputs.py \
    --input ../alphafold_outputs/ \
    --domains ../Domain_boundries.csv \
    --output ../calvados_simulations/ \
    --residues-csv /path/to/residues_CALVADOS3.csv \
    --seed seed-1_sample-0
```

**Output structure:**

```
calvados_simulations/
├── rrm1_rrm2_rrm3/
│   ├── config.yaml
│   ├── components.yaml
│   ├── domains.yaml
│   └── input/
│       └── structure.pdb
└── ...
```

---

## 7. Visualize Results

### Generate publication figures

```bash
python 10_comprehensive_visualization.py \
    --input ../ \
    --output ../visualizations/ \
    --use-builtin-sites
```

**Generated figures:**
- `plddt_histogram.pdf` - Confidence distribution
- `domain_distances_analysis.pdf` - Distance analysis
- `contact_sasa_filtered.pdf` - Contact-SASA plots
- `active_site_*.pdf` - Active site analyses
- `comprehensive_domain_grid.pdf` - Domain coverage heatmap

### Interactive dashboard

```bash
streamlit run interactive_dashboard.py
```

Opens web-based interface for:
- Browsing structure predictions
- Filtering by domain composition
- Interactive confidence visualization
- Domain distance exploration

---

## 8. Quick Reference

### Pipeline execution order

| Stage | Script | Input | Output |
|-------|--------|-------|--------|
| 1 | `01_generate_fasta.py` | FASTA + CSV | Domain FASTA files |
| 2 | `02_generate_af3_inputs.py` | FASTA files | JSON inputs |
| 3 | `03_run_alphafold.sh` | JSON files | CIF + JSON |
| 4 | `04_distance_analysis.py` | CIF files | Distance CSV + NPZ |
| 5 | `05_pLDDT_analysis.py` | CIF + JSON | Confidence CSV |
| 6 | `06_sasa_analysis.py` | CIF files | SASA CSV |
| 7 | `07_generate_calvados_inputs.py` | CIF files | CALVADOS YAML + PDB |
| 8 | `10_comprehensive_visualization.py` | CSV files | Publication PDFs |

### pLDDT confidence levels

| Level | Range | Interpretation |
|-------|-------|----------------|
| Very High | ≥90 | High confidence, well-ordered |
| Confident | 70-90 | Reasonable confidence |
| Low | 50-70 | Low confidence |
| Very Low | <50 | Likely unstructured |

### Domain constraints

| Constraint | Rule |
|------------|------|
| Sequential order | Domains must maintain native sequence order |
| KH1-6 group | Either all 6 domains or none |
| KH7a-KHb-KH8 | Must appear together or not at all |
| MD variants | MD1 and MD1L1 are mutually exclusive |

---

## 9. Common Analysis Tasks

### Load per-residue data

```python
import pandas as pd

# Load confidence scores
plddt = pd.read_csv('pLDDT_per_residue.csv')
print(plddt.head())

# Filter for specific structure
struct_data = plddt[plddt['structure'] == 'kh1_kh2_kh3_kh4_kh5_kh6']

# Get mean pLDDT per structure
mean_plddt = plddt.groupby(['structure', 'seed'])['plddt'].mean()
```

### Load distance matrices

```python
import numpy as np

data = np.load('distance_analysis/distance_matrices/kh1_kh2_kh3_kh4_kh5_kh6_distances.npz')
dist_matrix = data['distance_matrix']
residues = data['residue_numbers']
```

### Filter structures by confidence

```python
import pandas as pd

conf = pd.read_csv('confidence_scores.csv')

# High-confidence structures (mean pLDDT > 80)
high_conf = conf[conf['mean_plddt'] > 80]

# No clashes
no_clash = conf[conf['has_clash'] == False]
```

### Calculate domain statistics

```python
import pandas as pd

sasa = pd.read_csv('sasa_per_residue.csv')

# Mean SASA per domain
# (Requires merging with domain boundaries)
```

---

## 10. Troubleshooting

### AlphaFold3 runs out of memory

Reduce batch size or use a GPU with more VRAM:

```bash
export CUDA_VISIBLE_DEVICES=1  # Try different GPU
```

### Large CSV files slow to load

Use chunked reading:

```python
chunks = pd.read_csv('pLDDT_per_residue.csv', chunksize=100000)
for chunk in chunks:
    # Process chunk
    pass
```

Or filter during load:

```python
# Only load specific columns
plddt = pd.read_csv('pLDDT_per_residue.csv', usecols=['structure', 'plddt'])
```

### Missing predictions

Check `missing_sequences/` directory for structures that failed prediction.

Re-run specific structures:

```bash
python 01c_generate_missing_sequences.py
```

### FreeSASA installation issues

```bash
# On Ubuntu/Debian
sudo apt-get install libfreesasa-dev
pip install freesasa

# On macOS
brew install freesasa
pip install freesasa
```

---

## 11. Full Pipeline Script

Run the complete pipeline:

```bash
#!/bin/bash
cd /home/sbali/parp14/pipeline

# Stage 1-2: Generate inputs
python 01_generate_fasta.py --fasta ../PARP14.fasta --domains ../Domain_boundries.csv --output ../fasta_files/
python 02_generate_af3_inputs.py --fasta-dir ../fasta_files/ --output ../alphafold_inputs/ --seeds 5

# Stage 3: Run AlphaFold3 (takes hours/days)
bash 03_run_alphafold.sh

# Stage 4-6: Analysis
python 04_distance_analysis.py --input ../alphafold_outputs/ --domains ../Domain_boundries.csv --output ../distance_analysis/
python 05_pLDDT_analysis.py --input ../alphafold_outputs/ --output-residue ../pLDDT_per_residue.csv
python 06_sasa_analysis.py --input ../alphafold_outputs/ --output ../sasa_per_residue.csv

# Stage 7: CALVADOS preparation
python 07_generate_calvados_inputs.py --input ../alphafold_outputs/ --output ../calvados_simulations/

# Stage 8: Visualization
python 10_comprehensive_visualization.py --input ../ --output ../visualizations/

echo "Pipeline complete!"
```

Or use the automated script:

```bash
bash run_full_pipeline.sh
```

---

## 12. Proceed to Phase 2: CALVADOS Simulations

After completing the structure pipeline, proceed to coarse-grained MD simulations.

### Option A: Use pre-generated CALVADOS inputs

The script `07_generate_calvados_inputs.py` creates simulation-ready files:

```bash
# Navigate to a domain composition
cd /home/sbali/parp14/calvados_simulations/rrm1_rrm2_rrm3

# Run CALVADOS simulation
python /home/sbali/CALVADOS/calvados/sim.py --path .
```

### Option B: Filter by confidence first

```python
import pandas as pd

# Load confidence scores
conf = pd.read_csv('/home/sbali/parp14/confidence_scores.csv')

# Select high-quality structures for simulation
good = conf[(conf['mean_plddt'] > 70) & (conf['has_clash'] == False)]
print(f"Structures for CALVADOS: {len(good['structure'].unique())}")

# Save list for batch processing
good['structure'].unique().tolist()
```

### Phase 2 Documentation

See [CALVADOS PARP14 Project Documentation](/home/sbali/CALVADOS/docs/PARP14_PROJECT.md) for:

- Simulation protocol (50 ns × 25 replicates)
- SLURM batch job scripts
- Trajectory compilation
- Dynamics analysis (Rg, Ree, contacts)

---

## Related Documentation

| Document | Description |
|----------|-------------|
| [PARP14_PIPELINE_OVERVIEW.md](PARP14_PIPELINE_OVERVIEW.md) | Complete Phase 1 + Phase 2 workflow |
| [PARP14_ARCHITECTURE.md](PARP14_ARCHITECTURE.md) | Technical reference |
| [PARP14_ANALYSIS_GUIDE.md](PARP14_ANALYSIS_GUIDE.md) | Detailed analysis examples |
| [CALVADOS PARP14_PROJECT.md](/home/sbali/CALVADOS/docs/PARP14_PROJECT.md) | Phase 2: MD simulations |
