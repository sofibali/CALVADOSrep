# PARP14 AlphaFold3 Analysis Pipeline

Terminal-based pipeline for generating, running, and analyzing PARP14 domain combination structures with AlphaFold3.

## Pipeline Overview

```
1. Generate FASTA files          → 01_generate_fasta.py
2. Generate AF3 input JSONs       → 02_generate_af3_inputs.py
3. Run AlphaFold3 predictions     → 03_run_alphafold.sh
4. Distance analysis              → 04_distance_analysis.py
5. Confidence analysis            → 05_pLDDT_analysis.py
6. SASA analysis                  → 06_sasa_analysis.py
7. Generate CALVADOS inputs       → 07_generate_calvados_inputs.py
8. Visualization                  → 08_visualization.ipynb
```

## Requirements

```bash
# Python packages
pip install biopython numpy pandas matplotlib seaborn freesafa pyyaml

# System tools
# - AlphaFold3 (with alphafold3 command in PATH)
```

## Usage

### Step 1: Generate FASTA Files

```bash
python 01_generate_fasta.py \
    --fasta ../PARP14.fasta \
    --domains ../Domain_boundries.csv \
    --output fasta_files/ \
    --max-size 4
```

**Output:** `fasta_files/*.fasta` - One FASTA file per domain combination

### Step 2: Generate AlphaFold3 Input Files

```bash
python 02_generate_af3_inputs.py \
    --fasta-dir fasta_files/ \
    --output alphafold_inputs/ \
    --seeds 5
```

**Output:** `alphafold_inputs/*.json` - JSON files for AlphaFold3

### Step 3: Run AlphaFold3 Predictions 
## update input and output directories within the bash script 
```bash
bash 03_run_alphafold.sh 
```

**Output:** `alphafold_outputs/*/model.cif` - Predicted structures

### Step 4: Distance Analysis

```bash
python 04_distance_analysis.py \
    --input alphafold_outputs/ \
    --domains ../Domain_boundries.csv \
    --output distance_analysis/
```

**Output:**
- `distance_analysis/contact_maps/*.pdf` - Contact map visualizations
- `distance_analysis/distance_matrices/*.npz` - Distance matrices (reusable)
- `distance_analysis/domain_distances.csv` - Inter-domain distances

### Step 5: Comprehensive Confidence Analysis

```bash
python 05_pLDDT_analysis.py \
    --input alphafold_outputs/ \
    --output-residue pLDDT_per_residue.csv \
    --output-structure confidence_scores.csv
```

**Output:** 
- `pLDDT_per_residue.csv` - Per-residue pLDDT scores and confidence levels (Very High ≥90, Confident 70-90, Low 50-70, Very Low <50)
- `confidence_scores.csv` - Per-structure PTM, iPTM, ranking scores, mean pLDDT, and confidence distribution

**Note:** This single efficient script extracts both per-residue confidence (pLDDT from B-factor) and per-structure quality metrics (PTM/iPTM/ranking from JSON files) in one pass through the data.

### Step 6: SASA Analysis

```bash
python 06_sasa_analysis.py \
    --input alphafold_outputs/ \
    --output sasa_per_residue.csv \
    --domains ../Domain_boundries.csv
```

**Output:** `sasa_per_residue.csv` - Per-residue SASA values
```

**Output:** `sasa_per_residue.csv` - Per-residue SASA values

### Step 7: Generate CALVADOS Inputs

```bash
python 07_generate_calvados_inputs.py \
    --input alphafold_outputs/ \
    --domains ../Domain_boundries.csv \
    --output ../calvados_simulations/ \
    --residues-csv /path/to/residues_CALVADOS3.csv \
    --seed seed-1_sample-0
```

**Output:** For each structure, creates a directory with:
- `config.yaml` - Simulation parameters (box size, temperature, steps, etc.)
- `components.yaml` - Molecule definitions and restraints
- `domains.yaml` - Domain boundaries for multi-domain proteins
- `input/*.pdb` - Converted structure file (backbone atoms only)

These can be used directly with CALVADOS for coarse-grained MD simulations.

### Step 8: Visualization

```bash
jupyter notebook 08_visualization.ipynb
```

**Visualization outputs:**

The notebook generates comprehensive analysis plots from all pipeline outputs:
Open `08_visualization.ipynb` in Jupyter and run all cells to generate:
- Confidence score distributions (PTM, iPTM, ranking)
- pLDDT confidence level analysis (Very High/Confident/Low/Very Low)
- Domain distance analysis with contact density metrics
- Domain-domain interface heatmaps and scatter plots
- SASA exposure analysis
- Active site comparisons
- Integrated quality metrics correlating PTM, iPTM, pLDDT
- Comprehensive summary statistics and comparison tables

**Output PDFs:**
- `confidence_scores_distribution.pdf`
- `plddt_confidence_distribution.pdf`
- `confidence_level_analysis
**Output CSVs:**
- `structure_quality_metrics.csv` - Integrated metrics per structure
- `contact_density_metrics.csv` - Contact statistics per structure
- `top_50_structures.csv` - Best structures by ranking score
- `domain_pair_statistics.csv` - Aggregated domain pair distances

## Data Files

All intermediate data files are in machine-readable formats:

- **CSV files**: Domain distances, SASA, disorder analysis
- **NPZ files**: Distance matrices (NumPy compressed arrays)
- **PDF files**: Publication-quality vector graphics

## Re-analysis

To re-analyze without re-running AlphaFold3:

```bash
# Load distance matrices
import numpy as np
data = np.load('distance_analysis/distance_matrices/structure_name_distances.npz')
dist_matrix = data['distance_matrix']
residues = data['residue_numbers']

# Load CSV data
import pandas as pd
df = pd.read_csv('domain_distances.csv')
```

## Parallelization

The bash script (step 3) runs multiple AlphaFold3 jobs in parallel. Adjust:

```bash
MAX_PARALLEL=4  # Number of simultaneous predictions
GPU_ID=0        # GPU device ID
```

## File Structure

```
pipeline/
├── 01_generate_fasta.py
├── 02_generate_af3_inputs.py
├── 03_run_alphafold.sh
├── 04_distance_analysis.py
├── 04b_extract_confidence_scores.py
├── 05_disorder_analysis.py
├── 06_sasa_analysis.py
├── 07_generate_calvados_inputs.py
├── 08_visualization.ipynb
├── README.md
├── fasta_files/
├── alphafold_inputs/
├── alphafold_outputs/
├── distance_analysis/
│   ├── contact_maps/
│   ├── distance_matrices/
│   └── domain_distances.csv
├── d5_pLDDT_analysis.py
├── 06_sasa_analysis.py
├── 07_generate_calvados_inputs.py
├── 08_visualization.ipynb
├── README.md
├── fasta_files/
├── alphafold_inputs/
├── alphafold_outputs/
├── distance_analysis/
│   ├── contact_maps/
│   ├── distance_matrices/
│   └── domain_distances.csv
├── pLDDT
**AlphaFold3 command not found:**
Ensure AlphaFold3 is properly installed and `alphafold3` is in your PATH.
