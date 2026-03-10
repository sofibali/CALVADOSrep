# PARP14 AlphaFold3 Pipeline Architecture

**Phase 1** of the PARP14 Structure-to-Dynamics Pipeline

Systematic structural analysis of PARP14 domain combinations using AlphaFold3 predictions. This phase generates and validates structures that feed into [Phase 2: CALVADOS MD Simulations](/home/sbali/CALVADOS/docs/PARP14_PROJECT.md).

```
This Project (Phase 1)              Next Phase (Phase 2)
======================              ====================
/home/sbali/parp14                  /home/sbali/CALVADOS

AlphaFold3 structure prediction     CALVADOS coarse-grained MD
         │                                   │
         ▼                                   ▼
• Generate domain combinations      • Run 25 replicates per structure
• Predict 3D structures             • Compile 1 μs trajectories
• Analyze confidence (pLDDT)        • Calculate Rg, Ree, contacts
• Calculate distances & SASA        • Analyze structural dynamics
• Generate CALVADOS inputs ────────► Use structures for simulations
```

## Project Structure

```
parp14/
├── pipeline/                        # Analysis pipeline scripts
│   ├── 01_generate_fasta.py        # Generate FASTA from domain combinations
│   ├── 01b_generate_fasta_ignore_split.py
│   ├── 01c_generate_missing_sequences.py
│   ├── 02_generate_af3_inputs.py   # Create AlphaFold3 JSON inputs
│   ├── 03_run_alphafold.sh         # Sequential AF3 runner
│   ├── 03_run_alphafold_parallel.sh
│   ├── 04_distance_analysis.py     # Domain-domain distance analysis
│   ├── 05_pLDDT_analysis.py        # Confidence score extraction
│   ├── 06_sasa_analysis.py         # Solvent accessibility analysis
│   ├── 07_generate_calvados_inputs.py  # CALVADOS simulation setup
│   ├── 09_structure_inventory.py   # Generate coverage reports
│   ├── 10_comprehensive_visualization.py  # Publication figures
│   └── interactive_dashboard.py    # Streamlit web interface
├── alphafold_inputs/               # JSON files for AF3 (23 MB)
├── alphafold_outputs/              # Structure predictions (77 GB)
│   └── <domain_combination>/
│       ├── seed-1_sample-0/
│       │   ├── model.cif
│       │   ├── summary_confidences.json
│       │   ├── confidences.json
│       │   └── ranking_scores.json
│       └── seed-{2-5}_sample-0/
├── distance_analysis/              # Distance matrices & contact maps
├── calvados_simulations/           # CALVADOS-ready inputs (335 GB)
├── PARP14.fasta                    # Full 1801 aa sequence
├── Domain_boundries.csv            # Domain coordinate definitions
├── pLDDT_per_residue.csv          # Per-residue confidence (2 GB)
├── sasa_per_residue.csv           # Per-residue SASA (586 MB)
├── confidence_scores.csv           # Per-structure metrics
└── docs/                           # Documentation
```

---

## PARP14 Domain Architecture

Full-length PARP14 (1801 amino acids) contains 17 distinct domains:

| Domain | Residues | Length | Function |
|--------|----------|--------|----------|
| RRM1 | 1-145 | 145 | RNA Recognition Motif |
| RRM2 | 146-224 | 79 | RNA Recognition Motif |
| RRM3 | 225-314 | 90 | RNA Recognition Motif |
| KH1-KH6 | 315-737 | 423 | K-Homology domains (functional unit) |
| KH7a | 738-789 | 52 | K-Homology domain |
| MD1 | 790-981 | 192 | Macrodomain 1 |
| MD1L1 | 790-1004 | 215 | Macrodomain 1 (extended loop variant) |
| MD2 | 1005-1193 | 189 | Macrodomain 2 |
| MD3 | 1207-1388 | 182 | Macrodomain 3 |
| KHb-KH8 | 1389-1533 | 145 | K-Homology domains |
| WWE | 1534-1602 | 69 | WWE domain |
| ART | 1603-1801 | 199 | ADP-ribosyltransferase (catalytic) |

### Domain Combination Constraints

1. **Sequential Order**: Domains must maintain native sequence order
2. **KH1-6 Group**: All 6 domains together or none (functional unit)
3. **KH7a-KHb-KH8**: Must appear together or not at all
4. **MD Variants**: MD1 and MD1L1 are mutually exclusive alternatives

---

## Pipeline Stages

### Stage 1: FASTA Generation

**Script**: `01_generate_fasta.py`

Generates FASTA files for all valid domain combinations.

```python
# Example usage
python 01_generate_fasta.py \
    --fasta ../PARP14.fasta \
    --domains ../Domain_boundries.csv \
    --output fasta_files/ \
    --max-size 4
```

**Output**: ~1000 FASTA files named `<domain1>_<domain2>_...<domainN>.fasta`

---

### Stage 2: AlphaFold3 Input Generation

**Script**: `02_generate_af3_inputs.py`

Creates JSON input files for AlphaFold3 prediction.

**JSON Structure**:

```json
{
  "name": "kh1_kh2_kh3_kh4_kh5_kh6_md1",
  "sequences": [{
    "protein": {
      "id": ["A"],
      "sequence": "MVKL..."
    }
  }],
  "modelSeeds": [1, 2, 3, 4, 5],
  "dialect": "alphafold3",
  "version": 1
}
```

---

### Stage 3: AlphaFold3 Prediction

**Scripts**: `03_run_alphafold.sh`, `03_run_alphafold_parallel.sh`

Runs structure prediction with 5 model seeds per domain combination.

**Configuration**:

| Parameter | Value |
|-----------|-------|
| Database Directory | `/mnt/alphafold3` |
| Model Seeds | 5 per structure |
| GPU Device | CUDA (configurable) |
| Output Format | mmCIF + JSON confidence files |

**Output Files** (per seed):

| File | Description |
|------|-------------|
| `model.cif` | Predicted structure (mmCIF format) |
| `summary_confidences.json` | PTM, iPTM, ranking score |
| `confidences.json` | Detailed confidence metrics |
| `ranking_scores.json` | Model ranking information |

---

### Stage 4: Distance Analysis

**Script**: `04_distance_analysis.py`

Calculates inter-domain distances from predicted structures.

**Calculations**:
- Extract CA (alpha-carbon) coordinates from mmCIF
- Compute pairwise distance matrix (scipy vectorized)
- Identify domain boundaries and inter-domain gaps
- Generate contact maps (residues <8 Å apart)

**Output**:

| File | Description |
|------|-------------|
| `domain_distances.csv` | Min/max/mean distances per domain pair |
| `distance_matrices/*.npz` | NumPy compressed distance matrices |
| `contact_maps/*.pdf` | Contact map visualizations |

---

### Stage 5: Confidence Analysis

**Script**: `05_pLDDT_analysis.py`

Extracts and classifies AlphaFold3 confidence scores.

**pLDDT Classification**:

| Level | Range | Interpretation |
|-------|-------|----------------|
| Very High | ≥90 | High confidence, well-ordered |
| Confident | 70-90 | Reasonable confidence |
| Low | 50-70 | Low confidence, may be disordered |
| Very Low | <50 | Very low confidence, likely unstructured |

**Output Files**:

| File | Size | Rows | Content |
|------|------|------|---------|
| `pLDDT_per_residue.csv` | 2 GB | 22.4M | Per-residue confidence |
| `confidence_scores.csv` | 25 KB | 25K | Per-structure metrics |

**Per-Structure Metrics**:
- PTM (Predicted TM-score)
- iPTM (interface PTM, if multi-chain)
- Ranking score
- Mean pLDDT
- Confidence level distributions
- Clash detection flag

---

### Stage 6: SASA Analysis

**Script**: `06_sasa_analysis.py`

Calculates solvent-accessible surface area using FreeSASA.

**Output**: `sasa_per_residue.csv` (586 MB, 22.4M rows)

| Column | Description |
|--------|-------------|
| structure | Domain combination name |
| seed | Model seed identifier |
| residue | Residue number |
| residue_name | Three-letter amino acid code |
| sasa | Solvent accessible surface area (Å²) |

---

### Stage 7: CALVADOS Input Generation

**Script**: `07_generate_calvados_inputs.py`

Prepares structures for coarse-grained MD simulations.

**Process**:
1. Parse AlphaFold3 mmCIF structures
2. Convert to PDB (backbone atoms only: N, CA, C, O)
3. Generate YAML configuration files
4. Create domain boundary annotations

**Output Structure**:

```
calvados_simulations/<structure>/
├── config.yaml           # Simulation parameters
├── components.yaml       # Molecule definitions
├── domains.yaml          # Domain boundaries
└── input/
    └── *.pdb            # Backbone-only structures
```

**CALVADOS Configuration Parameters**:

| Parameter | Value | Description |
|-----------|-------|-------------|
| `box` | [300, 300, 300] | Simulation box (Å) |
| `temp` | 293 | Temperature (K) |
| `ionic` | 0.15 | Ionic strength (M) |
| `pH` | 7.0 | pH value |
| `topol` | center | Centering scheme |
| `wfreq` | 100 | Write frequency |
| `restraint_type` | harmonic | Restraint type |
| `k_harmonic` | 700.0 | Spring constant (kJ/mol/nm²) |
| `charge_termini` | both | Charged N and C termini |

---

## Active Site Definitions

Four catalytic domains with defined active site residues:

### MD1 (Macrodomain 1) - Residues 790-1004

| Set | Residues | Count |
|-----|----------|-------|
| Full | 822-836, 919-927, 961, 966 | 24 |
| Key | 831, 923, 962 | 3 |

### MD2 (Macrodomain 2) - Residues 1005-1193

| Set | Residues | Count |
|-----|----------|-------|
| Full | 1021-1024, 1034-1047, 1130-1141, 1170-1171, 1175, 1178 | 26 |
| Key | 1035, 1046, 1134, 1171 | 4 |

### MD3 (Macrodomain 3) - Residues 1230-1388

| Set | Residues | Count |
|-----|----------|-------|
| Full | 1235-1237, 1247-1261, 1302-1304, 1324-1337, 1369-1371, 1375 | 35 |
| Key | 1248, 1259, 1330, 1371 | 4 |

### ART (ADP-ribosyltransferase) - Residues 1603-1801

| Set | Residues | Count |
|-----|----------|-------|
| Full | 1681-1685, 1688, 1701, 1704-1709, 1714-1716, 1721-1722, 1726-1727, 1781 | 20 |
| Key | 1684, 1705, 1706, 1722 | 4 |

---

## Analysis Thresholds

| Metric | Threshold | Usage |
|--------|-----------|-------|
| Contact distance | 8.0 Å | Contact map generation |
| SASA exposed | >40.0 Å² | Exposed residue identification |
| pLDDT very high | ≥90 | High-confidence regions |
| pLDDT confident | 70-90 | Moderate confidence |
| pLDDT low | 50-70 | Low confidence |
| pLDDT very low | <50 | Likely disordered |

---

## Data Scale

| Component | Size | Count |
|-----------|------|-------|
| Domain combinations | - | 1,285 |
| Model seeds per structure | - | 5 |
| Total structure files | 77 GB | ~6,400 |
| Distance analysis | 199 GB | - |
| Per-residue data | ~2.6 GB | 45M rows |
| CALVADOS inputs | 335 GB | 484 sets |
| **Total** | **~615 GB** | - |

---

## Dependencies

### Python Packages

```
numpy >= 1.21.0
pandas >= 1.3.0
scipy >= 1.7.0
biopython >= 1.79
matplotlib >= 3.4.0
seaborn >= 0.11.0
tqdm >= 4.62.0
freesasa
pyyaml
streamlit (optional, for dashboard)
```

### External Tools

| Tool | Purpose |
|------|---------|
| AlphaFold3 | Structure prediction |
| FreeSASA | SASA calculations |
| CALVADOS | Coarse-grained MD |
| CUDA 12.4 | GPU acceleration |

---

## Connection to Phase 2: CALVADOS Simulations

The outputs from this pipeline feed directly into CALVADOS coarse-grained MD simulations:

### Key Outputs for CALVADOS

| Output | Location | Usage in CALVADOS |
|--------|----------|-------------------|
| `calvados_simulations/*/config.yaml` | Per-structure configs | Simulation parameters |
| `calvados_simulations/*/components.yaml` | Molecule definitions | Component setup |
| `calvados_simulations/*/domains.yaml` | Domain boundaries | Harmonic restraints |
| `calvados_simulations/*/input/*.pdb` | Backbone structures | Starting coordinates |
| `confidence_scores.csv` | Root directory | Structure filtering (pLDDT >70) |

### Running CALVADOS on Phase 1 Outputs

```bash
# Direct execution using pre-generated inputs
cd /home/sbali/parp14/calvados_simulations/rrm1_rrm2_rrm3
python /home/sbali/CALVADOS/calvados/sim.py --path .
```

See [CALVADOS PARP14 Project Documentation](/home/sbali/CALVADOS/docs/PARP14_PROJECT.md) for full simulation protocol.

---

## Related Documentation

| Document | Description |
|----------|-------------|
| [PARP14_PIPELINE_OVERVIEW.md](PARP14_PIPELINE_OVERVIEW.md) | Complete Phase 1 + Phase 2 workflow |
| [PARP14_QUICKSTART.md](PARP14_QUICKSTART.md) | Step-by-step guide for this pipeline |
| [PARP14_ANALYSIS_GUIDE.md](PARP14_ANALYSIS_GUIDE.md) | Detailed analysis examples |
| [CALVADOS PARP14_PROJECT.md](/home/sbali/CALVADOS/docs/PARP14_PROJECT.md) | Phase 2: MD simulations |

---

## References

- **AlphaFold3**: Abramson et al., Nature 2024
- **CALVADOS**: Tesei et al., PNAS 2024
- **PARP14**: Palazzo et al., Molecular Cell 2018
- **FreeSASA**: Mitternacht, F1000Research 2016
