#!/usr/bin/env python3
"""
Plot per-residue features (pLDDT, disorder, SASA) over the full-length PARP14 sequence.

Uses precomputed CSVs (disorder_analysis.csv, sasa_per_residue.csv) — no CIF parsing.
Maps construct-local residue numbers to FL positions via domain boundaries.
Averages across 25 replicates (5 seeds × 5 samples) per construct.

Usage:
    python plot_features_over_sequence.py                        # all 3 plots
    python plot_features_over_sequence.py --feature plddt        # just pLDDT
    python plot_features_over_sequence.py --feature sasa         # just SASA
    python plot_features_over_sequence.py --feature disorder     # just disorder
    python plot_features_over_sequence.py --constructs md1_md2_md3 art  # specific constructs
    python plot_features_over_sequence.py --min-domains 3        # only constructs with ≥3 domains
    python plot_features_over_sequence.py --has-domain md1       # only constructs containing md1
    python plot_features_over_sequence.py --fold-change          # plot FC relative to AF2 FL
"""
import os
import sys
import csv
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# ═══════════════════════════════════════════════════════════════════
# Paths
# ═══════════════════════════════════════════════════════════════════
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PARP14_DIR = os.path.dirname(SCRIPT_DIR)

DISORDER_CSV = os.path.join(PARP14_DIR, 'disorder_analysis.csv')
SASA_CSV = os.path.join(PARP14_DIR, 'sasa_per_residue.csv')
CONFIDENCE_CSV = os.path.join(PARP14_DIR, 'pipeline', 'confidence_scores.csv')
AF2_PDB = "/home/sbali/CALVADOS/examples/PARP14_MDP/input/parp14.pdb"
OUTPUT_DIR = os.path.join(SCRIPT_DIR, 'sequence_plots')

# ═══════════════════════════════════════════════════════════════════
# Domain definitions
# ═══════════════════════════════════════════════════════════════════
# FL boundaries (1-indexed, inclusive)
DOMAIN_DEFS = {
    'rrm1':    (1, 145),
    'rrm2':    (146, 224),
    'rrm3':    (225, 314),
    # Individual KH domains (used in construct names)
    'kh1':     (315, 384),
    'kh2':     (385, 454),
    'kh3':     (455, 520),
    'kh4':     (521, 593),
    'kh5':     (594, 665),
    'kh6':     (666, 737),
    # Grouped KH1-KH6 (for the 11-unit scheme)
    'kh1-kh6': (315, 737),
    'kh7a':    (738, 789),
    'md1':     (790, 981),   # truncated (old data)
    'md1l1':   (790, 1004),  # with linker (new data)
    'md2':     (1004, 1193),
    'md3':     (1207, 1388),
    'khb':     (1389, 1461),
    'kh8':     (1462, 1533),
    'khb-kh8': (1389, 1533),
    'wwe':     (1534, 1602),
    'art':     (1603, 1801),
}

# 11 grouped domain units for display
DOMAIN_GROUPS = [
    ('RRM1',    1, 145,   '#1f77b4'),
    ('RRM2',    146, 224, '#2ca02c'),
    ('RRM3',    225, 314, '#aec7e8'),
    ('KH1-KH6', 315, 737, '#ff7f0e'),
    ('KH7a',    738, 789, '#d62728'),
    ('MD1L1',   790, 1004, '#17becf'),
    ('MD2',     1004, 1193, '#9467bd'),
    ('MD3',     1207, 1388, '#8c564b'),
    ('KHb-KH8', 1389, 1533, '#e377c2'),
    ('WWE',     1534, 1602, '#7f7f7f'),
    ('ART',     1603, 1801, '#bcbd22'),
]

N_FL = 1801


# ═══════════════════════════════════════════════════════════════════
# Parse construct name → domain list
# ═══════════════════════════════════════════════════════════════════
def parse_construct_domains(name):
    """Parse construct directory name into ordered list of domain names.

    'kh1_kh2_kh3_kh4_kh5_kh6_kh7a_md1_md2' -> ['kh1','kh2',...,'md2']
    Greedy match, longest domain name first.
    """
    remaining = name.lower()
    domains = []
    known = sorted(DOMAIN_DEFS.keys(), key=len, reverse=True)
    while remaining:
        remaining = remaining.lstrip('_')
        if not remaining:
            break
        matched = False
        for dname in known:
            if remaining.startswith(dname):
                rest = remaining[len(dname):]
                if rest == '' or rest.startswith('_'):
                    domains.append(dname)
                    remaining = rest
                    matched = True
                    break
        if not matched:
            parts = remaining.split('_', 1)
            remaining = parts[1] if len(parts) > 1 else ''
    return domains


def construct_to_fl_mapping(domains):
    """Map construct-local residue indices (0-based) to FL residue indices (0-based).

    Handles overlapping boundaries (e.g. MD1L1 ends at 1004, MD2 starts at 1004).
    Returns array of FL indices, length = number of construct residues.
    """
    fl_indices = []
    prev_end = -1  # last FL index added (0-based)

    for dname in domains:
        if dname not in DOMAIN_DEFS:
            continue
        start, end = DOMAIN_DEFS[dname]
        for fi in range(start - 1, end):  # convert to 0-based
            if fi > prev_end:
                fl_indices.append(fi)
                prev_end = fi

    return np.array(fl_indices)


# ═══════════════════════════════════════════════════════════════════
# Load AF2 reference
# ═══════════════════════════════════════════════════════════════════
def load_af2_reference():
    """Extract per-residue pLDDT from AF2 PDB (B-factor column, CA atoms)."""
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('af2', AF2_PDB)

    plddts = []
    for residue in structure[0].get_residues():
        if residue.id[0] != ' ':
            continue
        cas = [a for a in residue if a.name == 'CA']
        if cas:
            plddts.append(cas[0].bfactor)
    return np.array(plddts)


def load_af2_sasa():
    """Compute per-residue SASA for AF2 reference structure."""
    import mdtraj as md
    traj = md.load(AF2_PDB)
    sasa = md.shrake_rupley(traj, mode='residue')[0]  # nm^2
    return sasa * 100  # convert to Å^2


# ═══════════════════════════════════════════════════════════════════
# Load precomputed data (chunked for memory efficiency)
# ═══════════════════════════════════════════════════════════════════
def load_construct_data(csv_path, value_col, construct_filter=None,
                        chunk_size=500000):
    """Load per-residue data from CSV, average across replicates.

    Returns dict: construct_name -> np.array of per-residue mean values
    (construct-local indexing).
    """
    print(f"  Loading {os.path.basename(csv_path)}...")

    # Determine column indices from header
    with open(csv_path) as f:
        header = f.readline().strip().split(',')

    struct_idx = header.index('structure')
    seed_idx = header.index('seed')
    resnum_idx = header.index('residue_number')
    val_idx = header.index(value_col)

    # Accumulate per-construct, per-residue values
    # construct -> {resnum -> [values across replicates]}
    data = defaultdict(lambda: defaultdict(list))
    n_rows = 0

    with open(csv_path) as f:
        f.readline()  # skip header
        for line in f:
            parts = line.strip().split(',')
            if len(parts) <= val_idx:
                continue

            struct = parts[struct_idx]
            if construct_filter and struct not in construct_filter:
                continue

            try:
                resnum = int(parts[resnum_idx])
                val = float(parts[val_idx])
            except (ValueError, IndexError):
                continue

            data[struct][resnum].append(val)
            n_rows += 1

            if n_rows % 5000000 == 0:
                print(f"    {n_rows/1e6:.0f}M rows...")

    # Average across replicates
    result = {}
    for struct, residues in data.items():
        max_res = max(residues.keys())
        arr = np.full(max_res, np.nan)
        for resnum, vals in residues.items():
            arr[resnum - 1] = np.mean(vals)  # 1-indexed -> 0-indexed
        result[struct] = arr

    print(f"    Loaded {n_rows} rows, {len(result)} constructs")
    return result


# ═══════════════════════════════════════════════════════════════════
# Plotting
# ═══════════════════════════════════════════════════════════════════
def add_domain_shading(ax):
    """Add domain background shading and labels inside the plot area."""
    for name, start, end, color in DOMAIN_GROUPS:
        ax.axvspan(start, end, alpha=0.08, color=color, zorder=0)

    # Add labels just below top of plot area using axes transform for y
    ylo, yhi = ax.get_ylim()
    y_label = yhi - 0.03 * (yhi - ylo)  # 3% below top
    for name, start, end, color in DOMAIN_GROUPS:
        ax.text((start + end) / 2, y_label, name,
                ha='center', va='top', fontsize=6, rotation=45, color=color,
                fontweight='bold')


def make_domain_legend():
    """Create legend elements for domain colors."""
    elements = []
    for name, start, end, color in DOMAIN_GROUPS:
        elements.append(plt.Rectangle((0, 0), 1, 1, fc=color, alpha=0.3, label=name))
    return elements


def plot_feature_scatter(construct_data, feature_name, ylabel, title, outpath,
                         fold_change=False, ref_data=None, color='steelblue',
                         ylim=None):
    """Scatter plot: each dot = one construct at one FL residue position.

    x-axis: FL residue index (1-1801)
    y-axis: feature value (or fold change if fold_change=True)
    """
    fig, ax = plt.subplots(figsize=(20, 6))

    # Collect all points
    xs = []
    ys = []

    for struct, values in construct_data.items():
        domains = parse_construct_domains(struct)
        if not domains:
            continue
        fl_map = construct_to_fl_mapping(domains)

        if len(fl_map) != len(values):
            # Length mismatch - skip
            continue

        for i, fi in enumerate(fl_map):
            if np.isnan(values[i]):
                continue

            if fold_change and ref_data is not None:
                ref_val = ref_data[fi]
                if ref_val > 0.01:  # avoid division by near-zero
                    ys.append(values[i] / ref_val)
                    xs.append(fi + 1)  # 1-indexed
            else:
                ys.append(values[i])
                xs.append(fi + 1)

    if not xs:
        print(f"  WARNING: No data points for {feature_name}")
        plt.close()
        return

    xs = np.array(xs)
    ys = np.array(ys)

    ax.scatter(xs, ys, s=0.3, alpha=0.1, c=color, rasterized=True, edgecolors='none')

    # Median line per residue
    residue_medians = defaultdict(list)
    for x, y in zip(xs, ys):
        residue_medians[x].append(y)

    med_x = sorted(residue_medians.keys())
    med_y = [np.median(residue_medians[x]) for x in med_x]
    ax.plot(med_x, med_y, color='darkred', linewidth=0.8, alpha=0.7, label='Median')

    if fold_change:
        ax.axhline(1.0, color='red', linewidth=0.8, linestyle='--', alpha=0.5, label='No change (FC=1)')

    ax.set_xlabel('Full-Length Residue Index', fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.set_xlim(1, N_FL)

    if ylim:
        ax.set_ylim(ylim)
    elif fold_change:
        p95 = np.percentile(ys, 95)
        p5 = np.percentile(ys, 5)
        margin = 0.15 * (p95 - p5)
        ax.set_ylim(max(0, p5 - margin), p95 + margin)

    # Add domain shading after ylim is set so labels position correctly
    add_domain_shading(ax)

    legend_elements = [Line2D([0], [0], color='darkred', linewidth=1, label='Median')]
    if fold_change:
        legend_elements.append(Line2D([0], [0], color='red', linestyle='--', label='FC=1'))
    legend_elements.extend(make_domain_legend())
    ax.legend(handles=legend_elements, loc='upper right', fontsize=6, ncol=4)

    plt.tight_layout()
    fig.savefig(outpath + '.png', dpi=200, bbox_inches='tight')
    fig.savefig(outpath + '.pdf', bbox_inches='tight')
    plt.close()

    print(f"  Saved: {outpath}.png ({len(construct_data)} constructs, {len(xs)} points)")


def plot_domain_sasa_boxplot(construct_data, outpath, fold_change=False, ref_sasa=None):
    """Per-domain SASA boxplot: each dot = one construct's total SASA for that domain.

    x-axis: domain name
    y-axis: total domain SASA (or fold change)
    """
    fig, ax = plt.subplots(figsize=(14, 6))

    # Collect per-domain SASA for each construct
    domain_values = defaultdict(list)

    for struct, values in construct_data.items():
        domains = parse_construct_domains(struct)
        if not domains:
            continue
        fl_map = construct_to_fl_mapping(domains)

        if len(fl_map) != len(values):
            continue

        # For each grouped domain, sum SASA of residues mapping to it
        for gname, gstart, gend, gcolor in DOMAIN_GROUPS:
            mask = (fl_map >= gstart - 1) & (fl_map < gend)
            if mask.sum() == 0:
                continue

            domain_sasa = np.nansum(values[mask])

            if fold_change and ref_sasa is not None:
                # Use ref SASA of only the FL residues present in this construct
                # (avoids bias from truncated domains like MD1 vs MD1L1)
                ref_val = np.sum(ref_sasa[fl_map[mask]])
                if ref_val > 0:
                    domain_values[gname].append(domain_sasa / ref_val)
            else:
                domain_values[gname].append(domain_sasa)

    # Plot
    positions = []
    labels = []
    for i, (gname, gstart, gend, gcolor) in enumerate(DOMAIN_GROUPS):
        if gname not in domain_values:
            continue
        fcs = domain_values[gname]
        positions.append(i)
        labels.append(gname)

        jitter = np.random.default_rng(42).uniform(-0.3, 0.3, len(fcs))
        ax.scatter(np.full(len(fcs), i) + jitter, fcs,
                   s=8, alpha=0.4, c=gcolor, edgecolors='none', rasterized=True)

        bp = ax.boxplot([fcs], positions=[i], widths=0.5, showfliers=False,
                        patch_artist=True, zorder=3)
        bp['boxes'][0].set_facecolor(gcolor)
        bp['boxes'][0].set_alpha(0.3)
        bp['medians'][0].set_color('black')

    if fold_change:
        ax.axhline(1.0, color='red', linewidth=0.8, linestyle='--', label='No change (FC=1)')

    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=10)

    ylabel = 'SASA Fold Change (construct / FL)' if fold_change else 'Total Domain SASA (Å²)'
    ax.set_ylabel(ylabel, fontsize=12)
    title = 'Per-Domain SASA Fold Change' if fold_change else 'Per-Domain Total SASA'
    ax.set_title(f'{title} Across Domain Deletion Constructs', fontsize=14)
    ax.legend(fontsize=10)

    plt.tight_layout()
    fig.savefig(outpath + '.png', dpi=200, bbox_inches='tight')
    fig.savefig(outpath + '.pdf', bbox_inches='tight')
    plt.close()

    print(f"  Saved: {outpath}.png")


# ═══════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════
def main():
    parser = argparse.ArgumentParser(
        description='Plot per-residue features over the PARP14 FL sequence')
    parser.add_argument('--feature', choices=['plddt', 'sasa', 'disorder', 'all'],
                        default='all', help='Which feature to plot')
    parser.add_argument('--fold-change', action='store_true',
                        help='Plot fold change relative to AF2 full-length')
    parser.add_argument('--constructs', nargs='+',
                        help='Only plot specific constructs (by name)')
    parser.add_argument('--min-domains', type=int, default=0,
                        help='Only constructs with ≥N domains')
    parser.add_argument('--has-domain', nargs='+',
                        help='Only constructs containing these domain(s)')
    parser.add_argument('--exclude-md1-truncated', action='store_true',
                        help='Exclude constructs with truncated MD1 (not MD1L1)')
    parser.add_argument('--output-dir', default=OUTPUT_DIR,
                        help='Output directory for plots')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    features = [args.feature] if args.feature != 'all' else ['plddt', 'disorder', 'sasa']

    # Build construct filter set if needed
    construct_filter = None
    if args.constructs:
        construct_filter = set(c.lower() for c in args.constructs)

    # Load AF2 reference if fold-change mode
    af2_plddt = af2_sasa = None
    if args.fold_change:
        print("Loading AF2 reference...")
        af2_plddt = load_af2_reference()
        print(f"  AF2 pLDDT: {len(af2_plddt)} residues, mean={af2_plddt.mean():.1f}")
        if 'sasa' in features:
            af2_sasa = load_af2_sasa()
            print(f"  AF2 SASA: mean={af2_sasa.mean():.1f} Å²")

    # Load data
    plddt_data = sasa_data = None

    if 'plddt' in features or 'disorder' in features:
        plddt_data = load_construct_data(DISORDER_CSV, 'plddt', construct_filter)

    if 'sasa' in features:
        sasa_data = load_construct_data(SASA_CSV, 'sasa', construct_filter)

    # Apply filters
    def filter_constructs(data):
        if data is None:
            return None
        filtered = {}
        for struct, values in data.items():
            domains = parse_construct_domains(struct)
            if args.min_domains and len(domains) < args.min_domains:
                continue
            if args.has_domain:
                if not all(d.lower() in domains for d in args.has_domain):
                    continue
            if args.exclude_md1_truncated:
                if 'md1' in domains and 'md1l1' not in struct:
                    continue
            filtered[struct] = values
        return filtered

    plddt_data = filter_constructs(plddt_data)
    sasa_data = filter_constructs(sasa_data)

    # Generate plots
    suffix = '_fc' if args.fold_change else ''

    if 'plddt' in features and plddt_data:
        print(f"\nPlotting pLDDT{' fold change' if args.fold_change else ''}...")
        plot_feature_scatter(
            plddt_data,
            'pLDDT',
            'pLDDT Fold Change (construct / FL)' if args.fold_change else 'pLDDT',
            f"Per-Residue pLDDT{'  Fold Change' if args.fold_change else ''} Across Constructs",
            os.path.join(args.output_dir, f'plddt_per_residue{suffix}'),
            fold_change=args.fold_change,
            ref_data=af2_plddt,
            color='steelblue',
        )

    if 'disorder' in features and plddt_data:
        print(f"\nPlotting disorder propensity{' fold change' if args.fold_change else ''}...")
        # Convert pLDDT to disorder propensity: 1 - pLDDT/100
        disorder_data = {}
        for struct, plddt in plddt_data.items():
            disorder_data[struct] = 1.0 - np.clip(plddt, 0, 100) / 100.0

        af2_disorder = (1.0 - np.clip(af2_plddt, 0, 100) / 100.0) if af2_plddt is not None else None

        plot_feature_scatter(
            disorder_data,
            'Disorder Propensity',
            'Disorder FC (construct / FL)' if args.fold_change else 'Disorder Propensity (1 - pLDDT/100)',
            f"Per-Residue Disorder Propensity{'  Fold Change' if args.fold_change else ''} Across Constructs",
            os.path.join(args.output_dir, f'disorder_per_residue{suffix}'),
            fold_change=args.fold_change,
            ref_data=af2_disorder,
            color='darkorange',
        )

    if 'sasa' in features and sasa_data:
        print(f"\nPlotting SASA{' fold change' if args.fold_change else ''}...")
        plot_feature_scatter(
            sasa_data,
            'SASA',
            'SASA Fold Change (construct / FL)' if args.fold_change else 'SASA (Å²)',
            f"Per-Residue SASA{'  Fold Change' if args.fold_change else ''} Across Constructs",
            os.path.join(args.output_dir, f'sasa_per_residue{suffix}'),
            fold_change=args.fold_change,
            ref_data=af2_sasa,
            color='forestgreen',
        )

        # Also make per-domain boxplot for SASA
        print(f"\nPlotting per-domain SASA boxplot...")
        plot_domain_sasa_boxplot(
            sasa_data,
            os.path.join(args.output_dir, f'sasa_per_domain{suffix}'),
            fold_change=args.fold_change,
            ref_sasa=af2_sasa,
        )

    print("\nDone!")


if __name__ == '__main__':
    main()
