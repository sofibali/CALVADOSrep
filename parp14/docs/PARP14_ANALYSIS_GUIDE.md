# PARP14 Analysis & Visualization Guide

**Phase 1** of the PARP14 Structure-to-Dynamics Pipeline

Comprehensive guide for analyzing AlphaFold3 predictions and generating publication-quality figures for the PARP14 domain library. Analysis results from this phase inform structure selection for [Phase 2: CALVADOS MD Simulations](/home/sbali/CALVADOS/docs/PARP14_PROJECT.md).

---

## Project Scale

| Metric | Value |
|--------|-------|
| Domain combinations | 1,285 |
| Model seeds per structure | 5 |
| Total structures | ~6,400 |
| Per-residue data rows | 22.4 million |
| Total storage | ~615 GB |

---

## Data Files Overview

### Primary Data Files

| File | Size | Rows | Description |
|------|------|------|-------------|
| `pLDDT_per_residue.csv` | 2 GB | 22.4M | Per-residue confidence scores |
| `sasa_per_residue.csv` | 586 MB | 22.4M | Per-residue solvent accessibility |
| `confidence_scores.csv` | 25 KB | 25K | Per-structure summary metrics |
| `distance_analysis/domain_distances.csv` | varies | varies | Inter-domain distance metrics |

### Per-Residue CSV Format

**pLDDT_per_residue.csv**
```
structure,seed,residue,plddt,confidence_level
kh1_kh2_kh3_kh4_kh5_kh6,seed-1_sample-0,1,80.5,Confident
kh1_kh2_kh3_kh4_kh5_kh6,seed-1_sample-0,2,85.2,Confident
...
```

**sasa_per_residue.csv**
```
structure,seed,residue,residue_name,sasa
kh1_kh2_kh3_kh4_kh5_kh6,seed-1_sample-0,1,MET,45.23
kh1_kh2_kh3_kh4_kh5_kh6,seed-1_sample-0,2,VAL,32.17
...
```

### Per-Structure CSV Format

**confidence_scores.csv**
```
structure,seed,ptm,iptm,ranking_score,has_clash,mean_plddt,num_residues,very_high_conf,confident,low_conf,very_low_conf,fraction_very_high,fraction_confident,fraction_low,fraction_very_low
kh1_kh2_kh3_kh4_kh5_kh6,seed-1_sample-0,0.72,,0.71,False,85.3,615,298,287,25,5,0.485,0.467,0.041,0.008
```

---

## Analysis Scripts

### Structure Inventory

Generate a comprehensive inventory of all predicted structures:

```bash
python 09_structure_inventory.py \
    --input ../alphafold_outputs/ \
    --output ../structure_inventory/
```

**Output:**
- Domain coverage heatmap
- Summary statistics per domain composition
- Missing structure report

### Comprehensive Visualization

Generate all publication-quality figures:

```bash
python 10_comprehensive_visualization.py \
    --input ../ \
    --output ../visualizations/ \
    --use-builtin-sites
```

**Generated Figures:**

| Figure | Description |
|--------|-------------|
| `plddt_histogram.pdf` | pLDDT distribution with confidence thresholds |
| `domain_distances_analysis.pdf` | Min/max/mean inter-domain distances |
| `contact_sasa_filtered.pdf` | Contact density vs SASA |
| `active_site_md1_full.pdf` | MD1 active site analysis |
| `active_site_md1_key.pdf` | MD1 key residues analysis |
| `active_site_md2_full.pdf` | MD2 active site analysis |
| `active_site_md2_key.pdf` | MD2 key residues analysis |
| `active_site_md3_full.pdf` | MD3 active site analysis |
| `active_site_md3_key.pdf` | MD3 key residues analysis |
| `active_site_art_full.pdf` | ART active site analysis |
| `active_site_art_key.pdf` | ART key residues analysis |
| `comprehensive_domain_grid.pdf` | Complete domain grid (pLDDT ± std) |
| `active_sites_combined.pdf` | All active sites summary |

---

## Python Analysis Examples

### Load and Explore Confidence Data

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load per-structure confidence
conf = pd.read_csv('confidence_scores.csv')

# Basic statistics
print(f"Total structures: {len(conf)}")
print(f"Mean pLDDT: {conf['mean_plddt'].mean():.1f}")
print(f"High confidence (>80): {(conf['mean_plddt'] > 80).sum()}")

# Plot distribution
plt.figure(figsize=(10, 6))
plt.hist(conf['mean_plddt'], bins=50, edgecolor='black')
plt.xlabel('Mean pLDDT')
plt.ylabel('Count')
plt.title('Distribution of Mean pLDDT Scores')
plt.axvline(90, color='green', linestyle='--', label='Very High (≥90)')
plt.axvline(70, color='orange', linestyle='--', label='Confident (70-90)')
plt.axvline(50, color='red', linestyle='--', label='Low (50-70)')
plt.legend()
plt.savefig('plddt_distribution.pdf')
```

### Analyze Specific Domain Combinations

```python
import pandas as pd

# Load per-residue data (chunked for large files)
def analyze_structure(structure_name):
    """Analyze a specific structure."""
    chunks = pd.read_csv('pLDDT_per_residue.csv', chunksize=100000)

    data = []
    for chunk in chunks:
        struct_data = chunk[chunk['structure'] == structure_name]
        if len(struct_data) > 0:
            data.append(struct_data)

    if data:
        return pd.concat(data)
    return pd.DataFrame()

# Example: analyze KH1-6 domain
kh_data = analyze_structure('kh1_kh2_kh3_kh4_kh5_kh6')
print(f"Residues: {kh_data['residue'].nunique()}")
print(f"Mean pLDDT: {kh_data['plddt'].mean():.1f}")
print(f"Seeds: {kh_data['seed'].unique()}")
```

### Compare Domain Combinations

```python
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

conf = pd.read_csv('confidence_scores.csv')

# Extract domain count from structure name
conf['domain_count'] = conf['structure'].apply(lambda x: len(x.split('_')))

# Plot pLDDT vs domain count
plt.figure(figsize=(10, 6))
sns.boxplot(data=conf, x='domain_count', y='mean_plddt')
plt.xlabel('Number of Domains')
plt.ylabel('Mean pLDDT')
plt.title('Prediction Confidence by Domain Count')
plt.savefig('plddt_by_domain_count.pdf')
```

### Load Distance Matrices

```python
import numpy as np
import matplotlib.pyplot as plt

# Load distance matrix
data = np.load('distance_analysis/distance_matrices/kh1_kh2_kh3_kh4_kh5_kh6_distances.npz')
dist_matrix = data['distance_matrix']
residues = data['residue_numbers']

# Plot contact map
plt.figure(figsize=(10, 10))
plt.imshow(dist_matrix < 8.0, cmap='Blues')
plt.colorbar(label='Contact (<8 Å)')
plt.xlabel('Residue')
plt.ylabel('Residue')
plt.title('Contact Map: KH1-KH6')
plt.savefig('kh_contact_map.pdf')
```

### Active Site Analysis

```python
import pandas as pd
import numpy as np

# Define active site residues
ACTIVE_SITES = {
    'MD1': {
        'full': list(range(822, 837)) + list(range(919, 928)) + [961, 966],
        'key': [831, 923, 962]
    },
    'MD2': {
        'full': list(range(1021, 1025)) + list(range(1034, 1048)) +
                list(range(1130, 1142)) + [1170, 1171, 1175, 1178],
        'key': [1035, 1046, 1134, 1171]
    },
    'MD3': {
        'full': list(range(1235, 1238)) + list(range(1247, 1262)) +
                list(range(1302, 1305)) + list(range(1324, 1338)) +
                list(range(1369, 1372)) + [1375],
        'key': [1248, 1259, 1330, 1371]
    },
    'ART': {
        'full': list(range(1681, 1686)) + [1688, 1701] + list(range(1704, 1710)) +
                list(range(1714, 1717)) + [1721, 1722, 1726, 1727, 1781],
        'key': [1684, 1705, 1706, 1722]
    }
}

def analyze_active_site(plddt_df, sasa_df, domain, residue_set='key'):
    """Analyze active site pLDDT and SASA for structures containing the domain."""
    residues = ACTIVE_SITES[domain][residue_set]

    # Filter structures containing this domain
    domain_lower = domain.lower()
    structures = plddt_df[plddt_df['structure'].str.contains(domain_lower, case=False)]

    # Filter for active site residues
    active_plddt = structures[structures['residue'].isin(residues)]
    active_sasa = sasa_df[sasa_df['residue'].isin(residues)]

    return {
        'mean_plddt': active_plddt['plddt'].mean(),
        'std_plddt': active_plddt['plddt'].std(),
        'mean_sasa': active_sasa['sasa'].mean(),
        'exposed_fraction': (active_sasa['sasa'] > 40).mean()
    }

# Example usage (with chunked loading for large files)
# results = analyze_active_site(plddt_df, sasa_df, 'MD1', 'key')
```

### Generate Domain Heatmap

```python
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

conf = pd.read_csv('confidence_scores.csv')

# Define all domains
DOMAINS = ['RRM1', 'RRM2', 'RRM3', 'KH1-KH6', 'KH7a', 'MD1', 'MD2', 'MD3', 'KHb-KH8', 'WWE', 'ART']

def check_domain_presence(structure, domain):
    """Check if domain is present in structure name."""
    domain_lower = domain.lower().replace('-', '_')
    struct_lower = structure.lower()
    return domain_lower in struct_lower

# Create domain presence matrix
for domain in DOMAINS:
    conf[f'has_{domain}'] = conf['structure'].apply(lambda x: check_domain_presence(x, domain))

# Aggregate by domain presence
domain_stats = []
for domain in DOMAINS:
    subset = conf[conf[f'has_{domain}']]
    if len(subset) > 0:
        domain_stats.append({
            'domain': domain,
            'mean_plddt': subset['mean_plddt'].mean(),
            'std_plddt': subset['mean_plddt'].std(),
            'n_structures': len(subset)
        })

stats_df = pd.DataFrame(domain_stats)

# Plot
plt.figure(figsize=(12, 6))
plt.bar(stats_df['domain'], stats_df['mean_plddt'], yerr=stats_df['std_plddt'], capsize=3)
plt.xlabel('Domain')
plt.ylabel('Mean pLDDT')
plt.title('Average pLDDT by Domain Presence')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig('domain_plddt_summary.pdf')
```

---

## Visualization Tools

### VMD (Visual Molecular Dynamics)

```tcl
# Load structure
mol new alphafold_outputs/kh1_kh2_kh3_kh4_kh5_kh6/seed-1_sample-0/model.cif

# Color by pLDDT (B-factor)
mol modcolor 0 0 Beta
mol modstyle 0 0 NewCartoon

# Save image
render TachyonInternal output.tga
```

### PyMOL

```python
# PyMOL script
load alphafold_outputs/kh1_kh2_kh3_kh4_kh5_kh6/seed-1_sample-0/model.cif
spectrum b, blue_white_red, minimum=50, maximum=90
cartoon automatic
ray 2400, 2400
png structure_plddt.png
```

### MDTraj for Trajectory Analysis

```python
import mdtraj as md

# Load structure
traj = md.load('alphafold_outputs/structure/seed-1_sample-0/model.cif')

# Calculate contacts
contacts = md.compute_contacts(traj, scheme='ca')

# Calculate radius of gyration
rg = md.compute_rg(traj)
print(f"Rg = {rg[0]:.2f} nm")
```

---

## Batch Analysis Scripts

### Analyze All Structures

```python
#!/usr/bin/env python
"""batch_analyze.py - Analyze all structures in parallel."""

import os
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

def analyze_single_structure(structure_path):
    """Analyze a single structure directory."""
    name = os.path.basename(structure_path)
    results = {'structure': name, 'seeds': 0, 'mean_plddt': None}

    for seed_dir in os.listdir(structure_path):
        if seed_dir.startswith('seed-'):
            results['seeds'] += 1
            # Add more analysis here

    return results

def main():
    base_dir = 'alphafold_outputs'
    structures = [os.path.join(base_dir, d) for d in os.listdir(base_dir)
                  if os.path.isdir(os.path.join(base_dir, d))]

    with ProcessPoolExecutor(max_workers=8) as executor:
        results = list(tqdm(executor.map(analyze_single_structure, structures),
                           total=len(structures)))

    df = pd.DataFrame(results)
    df.to_csv('batch_analysis_results.csv', index=False)
    print(f"Analyzed {len(df)} structures")

if __name__ == '__main__':
    main()
```

### Generate Summary Report

```python
#!/usr/bin/env python
"""generate_report.py - Generate summary report."""

import pandas as pd

# Load data
conf = pd.read_csv('confidence_scores.csv')
dist = pd.read_csv('distance_analysis/domain_distances.csv')

# Summary statistics
report = f"""
PARP14 AlphaFold3 Analysis Report
=================================

Total Structures: {conf['structure'].nunique()}
Total Models: {len(conf)}
Model Seeds: {conf['seed'].nunique()}

Confidence Metrics:
- Mean pLDDT: {conf['mean_plddt'].mean():.1f} ± {conf['mean_plddt'].std():.1f}
- High confidence (>80): {(conf['mean_plddt'] > 80).sum()} ({(conf['mean_plddt'] > 80).mean()*100:.1f}%)
- Very high confidence (>90): {(conf['mean_plddt'] > 90).sum()} ({(conf['mean_plddt'] > 90).mean()*100:.1f}%)
- With clashes: {conf['has_clash'].sum()}

PTM Scores:
- Mean PTM: {conf['ptm'].mean():.3f} ± {conf['ptm'].std():.3f}
- Max PTM: {conf['ptm'].max():.3f}

Top 5 Structures by pLDDT:
{conf.nlargest(5, 'mean_plddt')[['structure', 'seed', 'mean_plddt', 'ptm']].to_string()}
"""

print(report)
with open('analysis_report.txt', 'w') as f:
    f.write(report)
```

---

## Resource Requirements

### Storage

| Analysis | Storage Needed |
|----------|---------------|
| AlphaFold3 outputs | 77 GB |
| Distance matrices | 199 GB |
| Per-residue CSVs | ~2.6 GB |
| CALVADOS inputs | 335 GB |
| Visualizations | ~1-2 GB |
| **Total** | **~615 GB** |

### Computational Time

| Stage | Time per Structure | Total Time |
|-------|-------------------|------------|
| AlphaFold3 prediction | 1-10 min | ~200-500 hours |
| Distance analysis | 5-30 sec | ~2-3 hours |
| pLDDT analysis | 10-30 sec | ~3-4 hours |
| SASA analysis | 30-60 sec | ~6-8 hours |
| Visualization | - | ~1-2 hours |

### Memory Requirements

| Task | RAM Needed |
|------|------------|
| Load full pLDDT CSV | ~8-16 GB |
| Load full SASA CSV | ~4-8 GB |
| Distance matrix computation | ~2-4 GB per structure |
| Visualization | ~4-8 GB |

**Recommendation**: Use chunked loading for large CSV files on systems with <32 GB RAM.

---

## Interactive Dashboard

Launch the Streamlit dashboard for interactive exploration:

```bash
cd pipeline
streamlit run interactive_dashboard.py
```

**Features:**
- Structure browser with filtering
- Interactive confidence visualization
- Domain distance exploration
- Active site analysis
- Export capabilities

---

## Integration with CALVADOS (Phase 2)

After Phase 1 analysis, proceed to coarse-grained MD simulations in Phase 2.

### Structure Selection for Simulation

Use Phase 1 metrics to select high-quality structures:

```python
import pandas as pd

# Load confidence scores
conf = pd.read_csv('/home/sbali/parp14/confidence_scores.csv')

# Filter criteria for CALVADOS simulations
selected = conf[
    (conf['mean_plddt'] > 70) &      # Reasonable confidence
    (conf['has_clash'] == False) &    # No steric clashes
    (conf['ptm'] > 0.5)               # Good fold quality
]

print(f"Structures passing QC: {len(selected['structure'].unique())}")

# Export list for batch simulation
selected['structure'].unique().tofile('structures_for_calvados.txt', sep='\n')
```

### Run CALVADOS Simulations

```bash
# Use pre-generated inputs from Phase 1
cd /home/sbali/parp14/calvados_simulations/rrm1_rrm2_rrm3

# Run simulation
python /home/sbali/CALVADOS/calvados/sim.py --path .
```

### Cross-Phase Analysis

Combine Phase 1 (static) and Phase 2 (dynamic) metrics:

```python
import pandas as pd
import matplotlib.pyplot as plt

# Phase 1: AlphaFold3 confidence
phase1 = pd.read_csv('/home/sbali/parp14/confidence_scores.csv')

# Phase 2: CALVADOS dynamics (after simulations complete)
phase2 = pd.read_csv('/home/sbali/CALVADOS/PARP14_library/analysis/library_summary.csv')

# Merge and analyze
combined = phase2.merge(phase1[['structure', 'mean_plddt', 'ptm']],
                        left_on='composition', right_on='structure')

# Plot: Does AlphaFold confidence predict flexibility?
plt.figure(figsize=(10, 6))
plt.scatter(combined['mean_plddt'], combined['rg_mean'], alpha=0.5)
plt.xlabel('AlphaFold3 Mean pLDDT (Phase 1)')
plt.ylabel('CALVADOS Mean Rg nm (Phase 2)')
plt.title('Structure Confidence vs Simulated Compaction')
plt.savefig('phase1_vs_phase2_analysis.pdf')
```

### Phase 2 Documentation

- [CALVADOS PARP14 Project](/home/sbali/CALVADOS/docs/PARP14_PROJECT.md) - Full simulation protocol
- [CALVADOS Quickstart](/home/sbali/CALVADOS/docs/CALVADOS_QUICKSTART.md) - General CALVADOS usage
- [CALVADOS Architecture](/home/sbali/CALVADOS/docs/CALVADOS_ARCHITECTURE.md) - Technical reference

---

## Related Documentation

| Document | Description |
|----------|-------------|
| [PARP14_PIPELINE_OVERVIEW.md](PARP14_PIPELINE_OVERVIEW.md) | Complete Phase 1 + Phase 2 workflow |
| [PARP14_ARCHITECTURE.md](PARP14_ARCHITECTURE.md) | Technical reference for Phase 1 |
| [PARP14_QUICKSTART.md](PARP14_QUICKSTART.md) | Step-by-step pipeline guide |
| [CALVADOS PARP14_PROJECT.md](/home/sbali/CALVADOS/docs/PARP14_PROJECT.md) | Phase 2: MD simulations |
