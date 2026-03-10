# PARP14 Domain Composition Analysis Strategy

## Core Scientific Questions

1. **How does domain composition influence the structure/dynamics of MD1, MD2, MD3, and ART?**
2. **Does the presence/absence of neighboring domains affect active site accessibility?**
3. **Are there allosteric communication pathways between domains?**

---

## Analysis Framework

```
Domain Composition
       │
       ▼
┌──────────────────────────────────────────────────────────────┐
│                    STRUCTURAL METRICS                        │
├──────────────────────────────────────────────────────────────┤
│  Inter-Domain      │  Active Site       │  Domain            │
│  Distances         │  Accessibility     │  Flexibility       │
│                    │                    │                    │
│  • COM distances   │  • SASA            │  • pLDDT variance  │
│  • Min distances   │  • fpocket         │  • B-factor proxy  │
│  • Contact counts  │  • DoGSiteScorer   │  • Conformational  │
│                    │  • P2Rank          │    sampling (seeds)│
└──────────────────────────────────────────────────────────────┘
       │
       ▼
┌──────────────────────────────────────────────────────────────┐
│                    VISUALIZATION                              │
├──────────────────────────────────────────────────────────────┤
│  Heatmaps: Domain presence vs metric                         │
│  Boxplots: Metric distribution by composition                │
│  Network: Domain-domain interaction strengths                │
│  Ridge plots: Metric distributions across compositions       │
└──────────────────────────────────────────────────────────────┘
```

---

## Key Domains of Interest

| Domain | Residues | Function | Key Active Site Residues |
|--------|----------|----------|--------------------------|
| **MD1** | 790-981 | ADP-ribose hydrolase | 831, 923, 962 |
| **MD2** | 1005-1193 | ADP-ribose binding | 1035, 1046, 1134, 1171 |
| **MD3** | 1207-1388 | ADP-ribose binding | 1248, 1259, 1330, 1371 |
| **ART** | 1603-1801 | ADP-ribosyltransferase | 1684, 1705, 1706, 1722 |

---

## Analysis 1: Inter-Domain Distances & Contacts

### Metrics to Calculate

For each domain pair in each structure:

| Metric | Description | Interpretation |
|--------|-------------|----------------|
| **d_COM** | Center-of-mass distance | Overall domain separation |
| **d_min** | Minimum residue-residue distance | Closest approach |
| **d_interface** | Mean distance of interface residues | Interface compactness |
| **N_contacts** | Number of residue pairs <8Å | Interface size |
| **Contact_density** | N_contacts / (N_res1 × N_res2) | Normalized interface |

### Grouping by Domain Composition

```python
# Pseudocode for analysis structure
compositions = {
    'MD1_alone': ['md1'],
    'MD1_with_MD2': ['md1_md2', 'kh7a_md1_md2', ...],
    'MD1_with_MD2_MD3': ['md1_md2_md3', ...],
    'MD1_with_ART': ['md1_md2_md3_art', ...],
    # ... all combinations containing MD1
}

# For each target domain (MD1, MD2, MD3, ART):
#   For each composition group:
#     Calculate mean ± std of each metric
#     Compare across groups
```

### Key Comparisons

1. **Does MD1 change when MD2 is present?**
   - Compare MD1 metrics in `md1` vs `md1_md2` vs `md1_md2_md3`

2. **Does ART change when upstream domains are present?**
   - Compare ART metrics in `art` vs `wwe_art` vs `md3_wwe_art` vs full constructs

3. **Are there long-range effects?**
   - Does RRM1-3 presence affect MD1/MD2/MD3/ART?

---

## Analysis 2: Active Site Accessibility

### Tools to Use

| Tool | What it measures | Output |
|------|------------------|--------|
| **SASA** | Solvent accessible surface area | Å² per residue |
| **fpocket** | Pocket detection & druggability | Pocket volume, druggability score |
| **P2Rank** | ML-based ligand binding site prediction | Binding probability |
| **DoGSiteScorer** | Pocket geometry & druggability | DrugScore, volume |
| **CASTp** | Pocket volume & area | Å³ volume, Å² area |

### fpocket Analysis

```bash
# Install fpocket
conda install -c conda-forge fpocket

# Run on each structure
fpocket -f structure.pdb -o output_dir/

# Key outputs:
# - *_out/*_pockets.pqr: pocket atoms
# - *_out/*_info.txt: pocket descriptors
# - *_out/pocket*.pdb: individual pockets
```

### Metrics for Active Sites

For each active site (MD1, MD2, MD3, ART):

| Metric | Source | Interpretation |
|--------|--------|----------------|
| **Total SASA** | FreeSASA | Exposed surface area |
| **Pocket volume** | fpocket | Cavity size for ligand binding |
| **Druggability** | fpocket | Likelihood of being druggable |
| **Pocket depth** | fpocket | How buried the pocket is |
| **Hydrophobicity** | fpocket | Pocket character |
| **Pocket SASA** | fpocket | Pocket opening size |
| **P2Rank score** | P2Rank | ML prediction of binding |

### Active Site Definition File

```yaml
# /home/sbali/parp14/input/active_sites.yaml
MD1:
  domain_range: [790, 981]
  catalytic_residues: [831, 923, 962]
  binding_pocket_residues: [822, 823, 824, 825, 826, 827, 828, 829, 830, 831, 832, 833, 834, 835, 836, 919, 920, 921, 922, 923, 924, 925, 926, 927, 961, 966]

MD2:
  domain_range: [1005, 1193]
  catalytic_residues: [1035, 1046, 1134, 1171]
  binding_pocket_residues: [1021, 1022, 1023, 1024, 1034, 1035, 1036, 1037, 1038, 1039, 1040, 1041, 1042, 1043, 1044, 1045, 1046, 1047, 1130, 1131, 1132, 1133, 1134, 1135, 1136, 1137, 1138, 1139, 1140, 1141, 1170, 1171, 1175, 1178]

MD3:
  domain_range: [1207, 1388]
  catalytic_residues: [1248, 1259, 1330, 1371]
  binding_pocket_residues: [1235, 1236, 1237, 1247, 1248, 1249, 1250, 1251, 1252, 1253, 1254, 1255, 1256, 1257, 1258, 1259, 1260, 1261, 1302, 1303, 1304, 1324, 1325, 1326, 1327, 1328, 1329, 1330, 1331, 1332, 1333, 1334, 1335, 1336, 1337, 1369, 1370, 1371, 1375]

ART:
  domain_range: [1603, 1801]
  catalytic_residues: [1684, 1705, 1706, 1722]
  binding_pocket_residues: [1681, 1682, 1683, 1684, 1685, 1688, 1701, 1704, 1705, 1706, 1707, 1708, 1709, 1714, 1715, 1716, 1721, 1722, 1726, 1727, 1781]
```

---

## Analysis 3: Conformational Variability Across Seeds

AlphaFold3 generates 5 seeds per structure - use this for pseudo-dynamics:

| Metric | Calculation | Interpretation |
|--------|-------------|----------------|
| **pLDDT variance** | std(pLDDT) across seeds | Prediction uncertainty |
| **RMSD across seeds** | Pairwise RMSD of seeds | Conformational flexibility |
| **Distance variance** | std(d_COM) across seeds | Domain arrangement flexibility |
| **Pocket variance** | std(volume) across seeds | Active site flexibility |

---

## Visualization Strategy

### Figure 1: Domain Composition Effect Matrix

**Heatmap** showing how each metric changes with domain composition:

```
                    Metric
         ┌─────────────────────────────┐
         │ pLDDT │ SASA │ Pocket │ d_COM│
    ─────┼───────┼──────┼────────┼──────┤
    MD1  │  0.92 │ 245  │  890   │  --  │
Comp A   │       │      │        │      │
    MD2  │  0.88 │ 312  │  1020  │ 25.3 │
    ─────┼───────┼──────┼────────┼──────┤
    MD1  │  0.85 │ 198  │  720   │  --  │
Comp B   │       │      │        │      │
    MD2  │  0.91 │ 287  │  950   │ 18.7 │
    ─────┴───────┴──────┴────────┴──────┘
```

### Figure 2: Active Site Accessibility by Composition

**Grouped boxplots** for each active site:

```
              MD1 Active Site
    ┌─────────────────────────────────┐
    │    ┌─┐                          │
SASA│ ┌─┬┤ ├┬─┐  ┌─┬─┬─┐  ┌─┬─┬─┐    │
(Å²)│ └─┴┴─┴┴─┘  └─┴┬┴─┘  └─┴┬┴─┘    │
    │    alone    +MD2    +MD2+MD3   │
    └─────────────────────────────────┘
```

### Figure 3: Inter-Domain Distance Network

**Network diagram** showing domain-domain distances:

```
    RRM1 ─── RRM2 ─── RRM3
                        │
                       ╱
    KH1-6 ──────────────
       │
    KH7a
       │
      MD1 ══════ MD2 ══════ MD3
                              │
                            WWE
                              │
                            ART

Edge thickness = contact density
Edge color = distance (blue=close, red=far)
```

### Figure 4: Pocket Volume Distribution

**Ridge plot** showing pocket volume distribution across compositions:

```
    MD1 alone      ╱╲
    MD1+MD2       ╱  ╲____
    MD1+MD2+MD3  ╱    ╲
    Full         ╱      ╲___
                ────────────────
                Pocket Volume (Å³)
```

### Figure 5: Seed-to-Seed Variability

**Violin plots** showing conformational sampling:

```
              RMSD Across Seeds (Å)
    ┌─────────────────────────────────┐
    │   ◇      ◇◇     ◇      ◇◇◇     │
    │  ╱ ╲    ╱  ╲   ╱ ╲    ╱   ╲    │
    │ ╱   ╲  ╱    ╲ ╱   ╲  ╱     ╲   │
    │╱     ╲╱      ╲     ╲╱       ╲  │
    │ MD1   MD2    MD3    ART       │
    └─────────────────────────────────┘
```

---

## Recommended Tool Installation

```bash
# Create analysis environment
conda create -n parp14-analysis python=3.10
conda activate parp14-analysis

# Core packages
pip install numpy pandas scipy matplotlib seaborn
pip install biopython mdtraj prody

# Structure analysis
conda install -c conda-forge fpocket
pip install p2rank  # or download from GitHub
pip install freesasa

# Visualization
pip install plotly altair
pip install scikit-learn  # for clustering

# Optional: PyMOL for visualization
conda install -c conda-forge pymol-open-source
```

---

## Analysis Scripts Needed

| Script | Purpose | Priority |
|--------|---------|----------|
| `08_fpocket_analysis.py` | Run fpocket on all structures | HIGH |
| `09_interdomain_dynamics.py` | Calculate all inter-domain metrics | HIGH |
| `10_active_site_analysis.py` | Combine SASA + fpocket for active sites | HIGH |
| `11_composition_effects.py` | Group by composition, calculate statistics | HIGH |
| `12_publication_figures.py` | Generate all figures described above | HIGH |
| `13_statistical_tests.py` | Significance testing between groups | MEDIUM |

---

## Expected Outputs

### CSV Files
- `results/per_structure/interdomain_distances.csv`
- `results/per_structure/fpocket_results.csv`
- `results/per_structure/active_site_accessibility.csv`
- `results/per_structure/composition_effects_summary.csv`

### Figures
- `figures/domain_composition_effects/effect_matrix_heatmap.pdf`
- `figures/active_site_accessibility/md1_accessibility_by_composition.pdf`
- `figures/active_site_accessibility/md2_accessibility_by_composition.pdf`
- `figures/active_site_accessibility/md3_accessibility_by_composition.pdf`
- `figures/active_site_accessibility/art_accessibility_by_composition.pdf`
- `figures/interdomain_contacts/distance_network.pdf`
- `figures/interdomain_contacts/contact_heatmap.pdf`
- `figures/summary_statistics/seed_variability.pdf`
