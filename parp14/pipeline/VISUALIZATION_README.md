# Comprehensive PARP14 Visualization

## Quick Start

Once the analysis pipeline completes, run the comprehensive visualization:

```bash
cd /home/sbali/parp14/pipeline

# Basic visualization (without active sites)
python 10_comprehensive_visualization.py \
    --input . \
    --output visualizations/

# With built-in PARP14 active sites (RECOMMENDED - analyzes all 4 domains × 2 sets = 8 analyses)
python 10_comprehensive_visualization.py \
    --input . \
    --output visualizations/ \
    --use-builtin-sites

# With custom active sites from JSON file
python 10_comprehensive_visualization.py \
    --input . \
    --output visualizations/ \
    --active-sites-json my_active_sites.json
```

## Built-in PARP14 Active Sites

The `--use-builtin-sites` flag analyzes **all 4 PARP14 domains** with **both full and key residue sets**:

### MD1 Domain
- **Full** (24 residues): 822,823,824,828,829,830,831,832,833,834,835,836,919,920,921,922,923,924,925,926,927,961,962,966
- **Key** (3 residues): 831,923,962

### MD2 Domain
- **Full** (26 residues): 1021,1022,1023,1024,1034,1035,1036,1043,1044,1045,1046,1047,1130,1131,1132,1133,1134,1135,1136,1137,1138,1141,1170,1171,1175,1178
- **Key** (4 residues): 1035,1046,1134,1171

### MD3 Domain
- **Full** (35 residues): 1235,1236,1237,1247,1248,1249,1254,1255,1256,1257,1258,1259,1260,1261,1302,1303,1304,1324,1325,1326,1327,1328,1329,1330,1331,1332,1333,1334,1335,1336,1337,1369,1370,1371,1375
- **Key** (4 residues): 1248,1259,1330,1371

### ART Domain
- **Full** (20 residues): 1681,1682,1683,1684,1685,1688,1701,1704,1705,1706,1707,1709,1714,1715,1716,1721,1722,1726,1727,1781
- **Key** (4 residues): 1684,1705,1706,1722

## Custom Active Sites (JSON Format)

Create a JSON file (e.g., `my_active_sites.json`):

```json
{
  "Site_A": [100, 101, 102, 150, 151],
  "Site_B": [200, 201, 202, 250, 251],
  "Site_C": [300, 350, 400]
}
```

## Output Files

The script generates:

### Core Visualizations
1. **plddt_histogram.pdf** - Distribution of pLDDT scores with confidence thresholds (no pie chart)
2. **domain_distances_analysis.pdf** - Min/max distances and close contacts (<8Å) between domains
3. **contact_sasa_filtered.pdf** - Close contacts and high SASA residues (filtered)
4. **comprehensive_domain_grid.pdf** - Grid plot of all domain compositions with pLDDT ± std dev

### Active Site Analyses (when --use-builtin-sites or --active-sites-json is used)
5. **active_site_md1_full.pdf** - MD1 full active site (24 residues)
6. **active_site_md1_key.pdf** - MD1 key residues (3 residues)
7. **active_site_md2_full.pdf** - MD2 full active site (26 residues)
8. **active_site_md2_key.pdf** - MD2 key residues (4 residues)
9. **active_site_md3_full.pdf** - MD3 full active site (35 residues)
10. **active_site_md3_key.pdf** - MD3 key residues (4 residues)
11. **active_site_art_full.pdf** - ART full active site (20 residues)
12. **active_site_art_key.pdf** - ART key residues (4 residues)
13. **active_sites_combined_summary.pdf** - Combined comparison of all active sites

### CSV Summaries
- domain_distance_summary.csv
- domain_composition_summary.csv
- active_sites_summary.csv (when active sites analyzed)

## Monitoring the Pipeline

Check if analysis is still running:
```bash
ps aux | grep run_analysis_pipeline
```

View progress:
```bash
tail -f analysis_pipeline_*.log
```

## Domain Order

Domains are plotted in standard PARP14 order:
- RRM1, RRM2, RRM3
- KH1-KH6
- KH7a
- MD1, MD1L1 (MD1L1 = complete loop)
- MD2, MD3
- KHb, KH8
- WWE
- ART

## Customization

Edit the thresholds in the script:
- Contact distance: default 8.0 Å (line ~185)
- SASA threshold: default 40.0 Ų (line ~185)
- Confidence thresholds: 50, 70, 90 (standard AlphaFold)
