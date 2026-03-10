# PARP14 Structure Analysis Scripts

This directory contains Python scripts for comprehensive analysis of AlphaFold3-predicted PARP14 domain structures.

## Available Analysis Scripts

### 1. `analyze_structures.py`
**Purpose**: Extract and analyze overall structure quality metrics

**Metrics analyzed:**
- pLDDT (per-residue confidence scores)
- PTM (predicted TM-score)
- iPTM (interface predicted TM-score)
- Ranking scores
- Fraction of disordered residues

**Usage:**
```bash
python analyze_structures.py --output_dir ../alphafold_outputs --results_dir ../analysis_results
```

**Outputs:**
- `structure_metrics.csv` - Overall quality metrics for each structure
- `plddt_distribution.png` - Distribution plots
- `structure_analysis_summary.txt` - Text summary report

---

### 2. `extract_plddt.py`
**Purpose**: Extract per-residue pLDDT scores and create detailed visualizations

**Usage:**
```bash
python extract_plddt.py --output_dir ../alphafold_outputs --results_dir ../analysis_results
```

**Outputs:**
- `plddt_per_residue.csv` - Detailed per-residue pLDDT scores
- `plddt_summary.csv` - Summary statistics
- `plddt_plots/` - Individual per-structure pLDDT plots
- `plddt_summary_plots.png` - Comprehensive summary visualizations

**Confidence levels:**
- Very high: pLDDT ≥ 90 (dark blue)
- Confident: pLDDT 70-90 (light blue)
- Low: pLDDT 50-70 (yellow)
- Very low: pLDDT < 50 (orange)

---

### 3. `analyze_contacts.py`
**Purpose**: Generate contact maps and distance matrices

**Parameters:**
- `--distance_cutoff`: Distance threshold for contacts (default: 8.0 Å)

**Usage:**
```bash
python analyze_contacts.py --output_dir ../alphafold_outputs --results_dir ../analysis_results --distance_cutoff 8.0
```

**Outputs:**
- `contact_maps/` - Contact map images for each structure
- `distance_maps/` - Distance matrix visualizations
- `contact_analysis_summary.csv` - Contact statistics

**What is analyzed:**
- CA-CA distances between all residue pairs
- Total number of contacts
- Contact density (normalized by structure size)
- Average contacts per residue

---

### 4. `analyze_interfaces.py`
**Purpose**: Analyze domain-domain interfaces in multi-domain structures

**Parameters:**
- `--distance_cutoff`: Distance threshold for interface definition (default: 5.0 Å)

**Usage:**
```bash
python analyze_interfaces.py --output_dir ../alphafold_outputs --results_dir ../analysis_results --distance_cutoff 5.0
```

**Outputs:**
- `interface_analysis/` - Network plots showing domain interactions
- `interface_analysis_summary.csv` - Interface statistics
- `interface_residues_detailed.csv` - Detailed residue-level interface data

**What is analyzed:**
- Inter-domain contacts (residues from different domains)
- Number of interface residues per domain pair
- Average interface distances
- Interface network topology

---

### 5. `run_full_analysis.py`
**Purpose**: Master script to run all analyses in sequence

**Usage:**
```bash
# Run all analyses
python run_full_analysis.py --output_dir ../alphafold_outputs --results_dir ../analysis_results

# Skip specific analyses
python run_full_analysis.py --skip_contacts --skip_interfaces
```

**Options:**
- `--skip_structures` - Skip structure quality analysis
- `--skip_plddt` - Skip pLDDT extraction
- `--skip_contacts` - Skip contact analysis
- `--skip_interfaces` - Skip interface analysis

**Outputs:**
- All outputs from individual scripts
- `MASTER_ANALYSIS_REPORT.md` - Comprehensive markdown report

---

## Quick Start

### Run Complete Analysis Pipeline

```bash
cd analysis/
python run_full_analysis.py
```

This will:
1. Analyze all structures in `../alphafold_outputs/`
2. Save results to `../analysis_results/`
3. Generate visualizations and summary reports

### Run Individual Analyses

```bash
# Just extract pLDDT scores
python extract_plddt.py

# Just generate contact maps
python analyze_contacts.py

# Just analyze interfaces
python analyze_interfaces.py
```

---

## Output Structure

After running analyses, the results directory will contain:

```
analysis_results/
├── structure_metrics.csv
├── plddt_summary.csv
├── plddt_per_residue.csv
├── contact_analysis_summary.csv
├── interface_analysis_summary.csv
├── interface_residues_detailed.csv
├── plddt_distribution.png
├── plddt_summary_plots.png
├── structure_analysis_summary.txt
├── MASTER_ANALYSIS_REPORT.md
├── contact_maps/
│   ├── rrm1_rrm2_contact.png
│   ├── rrm1_md1_contact.png
│   └── ...
├── distance_maps/
│   ├── rrm1_rrm2_distance.png
│   └── ...
├── plddt_plots/
│   ├── rrm1_rrm2_plddt.png
│   └── ...
└── interface_analysis/
    ├── rrm1_rrm2_interface_network.png
    └── ...
```

---

## Requirements

All scripts require:
- Python 3.8+
- numpy
- pandas
- matplotlib
- seaborn
- biopython
- scipy
- networkx (for interface network plots)

Install via:
```bash
pip install -r ../requirements.txt
```

---

## Customization

### Adjust Distance Cutoffs

```bash
# Stricter contact definition (6Å instead of 8Å)
python analyze_contacts.py --distance_cutoff 6.0

# Looser interface definition (7Å instead of 5Å)
python analyze_interfaces.py --distance_cutoff 7.0
```

### Process Specific Structures

Modify the glob pattern in each script to process specific structures:

```python
# In any analysis script, modify:
pred_dirs = glob.glob(os.path.join(output_dir, 'rrm*'))  # Only RRM-containing
```

---

## Tips

1. **Performance**: Analyzing large numbers of structures can take time. Use the master script to run all analyses overnight.

2. **Memory**: Each script processes structures sequentially to minimize memory usage.

3. **Visualization**: All plots are saved as high-resolution PNG files (300 DPI) suitable for publication.

4. **Data Format**: All CSV files use standard formats and can be easily imported into Excel, R, or pandas for further analysis.

5. **Debugging**: Each script has verbose output. Run with `--help` to see all options.

---

## Troubleshooting

### "No CIF file found"
- Ensure AlphaFold3 predictions completed successfully
- Check that CIF files exist in subdirectories of `alphafold_outputs/`

### "Module not found" errors
- Install requirements: `pip install -r ../requirements.txt`

### Memory issues
- Process fewer structures at once
- Reduce image resolution in plotting functions

### Interface analysis shows no interfaces
- Check that multi-domain structures exist
- Verify domain boundary parsing is correct
- Try increasing `--distance_cutoff`

---

## Citation

If you use these analysis scripts, please cite the PARP14 AlphaFold3 Analysis Pipeline repository and AlphaFold3:

```
Abramson, J., et al. (2024). Accurate structure prediction of biomolecular 
interactions with AlphaFold 3. Nature.
```
