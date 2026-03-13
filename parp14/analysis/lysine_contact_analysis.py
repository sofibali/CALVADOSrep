#!/usr/bin/env python3
"""
Lysine contact analysis for PARP14 AF3 domain-deletion library.

For each construct (best-ranked model per seed), compute:
  1. Lys-Lys contacts: Lys NZ atoms within 16 Å of another Lys NZ
  2. Lys-acidic contacts: Lys NZ within 8 Å of Asp OD1/OD2 or Glu OE1/OE2

Maps residues to FL positions and plots per-residue contact counts over sequence.

Uses BioPython CIF parsing (no mdtraj/SASA, so fast).

Usage:
    python lysine_contact_analysis.py                    # all constructs
    python lysine_contact_analysis.py --workers 32       # more parallelism
    python lysine_contact_analysis.py --constructs md1_md2_md3
"""
import os
import sys
import glob
import argparse
import numpy as np
import warnings
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

warnings.filterwarnings('ignore')

# ═══════════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════════
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PARP14_DIR = os.path.dirname(SCRIPT_DIR)

AF3_DIR = os.path.join(PARP14_DIR, 'alphafold_outputs')
AF2_PDB = "/home/sbali/CALVADOS/examples/PARP14_MDP/input/parp14.pdb"
OUTPUT_DIR = os.path.join(SCRIPT_DIR, 'lysine_contact_plots')
CACHE_DIR = os.path.join(SCRIPT_DIR, 'lysine_contact_cache')

LYS_LYS_CUTOFF = 16.0  # Å
LYS_ACIDIC_CUTOFF = 8.0  # Å

N_FL = 1801

# Domain definitions (same as plot_features_over_sequence.py)
DOMAIN_DEFS = {
    'rrm1': (1, 145), 'rrm2': (146, 224), 'rrm3': (225, 314),
    'kh1': (315, 384), 'kh2': (385, 454), 'kh3': (455, 520),
    'kh4': (521, 593), 'kh5': (594, 665), 'kh6': (666, 737),
    'kh1-kh6': (315, 737), 'kh7a': (738, 789),
    'md1': (790, 981), 'md1l1': (790, 1004),
    'md2': (1004, 1193), 'md3': (1207, 1388),
    'khb': (1389, 1461), 'kh8': (1462, 1533), 'khb-kh8': (1389, 1533),
    'wwe': (1534, 1602), 'art': (1603, 1801),
}

DOMAIN_GROUPS = [
    ('RRM1', 1, 145, '#1f77b4'), ('RRM2', 146, 224, '#2ca02c'),
    ('RRM3', 225, 314, '#aec7e8'), ('KH1-KH6', 315, 737, '#ff7f0e'),
    ('KH7a', 738, 789, '#d62728'), ('MD1L1', 790, 1004, '#17becf'),
    ('MD2', 1004, 1193, '#9467bd'), ('MD3', 1207, 1388, '#8c564b'),
    ('KHb-KH8', 1389, 1533, '#e377c2'), ('WWE', 1534, 1602, '#7f7f7f'),
    ('ART', 1603, 1801, '#bcbd22'),
]


# ═══════════════════════════════════════════════════════════════════
# Domain parsing and FL mapping (shared with plot_features_over_sequence.py)
# ═══════════════════════════════════════════════════════════════════
def parse_construct_domains(name):
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
    fl_indices = []
    prev_end = -1
    for dname in domains:
        if dname not in DOMAIN_DEFS:
            continue
        start, end = DOMAIN_DEFS[dname]
        for fi in range(start - 1, end):
            if fi > prev_end:
                fl_indices.append(fi)
                prev_end = fi
    return np.array(fl_indices)


# ═══════════════════════════════════════════════════════════════════
# Contact computation from a structure file
# ═══════════════════════════════════════════════════════════════════
def compute_lysine_contacts_from_structure(structure):
    """Compute Lys-Lys and Lys-acidic contacts from a BioPython structure.

    Returns:
        lys_lys_counts: dict {resnum(0-based): count of Lys NZ within 16Å}
        lys_acidic_counts: dict {resnum(0-based): count of D/E within 8Å of Lys NZ}
        n_residues: total number of residues
    """
    model = structure[0]

    # Collect all residues with their key atoms
    residues_info = []  # list of (resnum_0based, resname, atom_coords_dict)
    for chain in model:
        for residue in chain:
            if residue.id[0] != ' ':
                continue
            resnum = residue.id[1] - 1  # 0-based
            resname = residue.resname.strip()
            atoms = {}
            for atom in residue:
                atoms[atom.name] = atom.get_vector().get_array()
            residues_info.append((resnum, resname, atoms))

    n_residues = len(residues_info)

    # Extract Lys NZ positions
    lys_nz = []  # (idx_in_residues_info, coord)
    for idx, (resnum, resname, atoms) in enumerate(residues_info):
        if resname == 'LYS' and 'NZ' in atoms:
            lys_nz.append((idx, atoms['NZ']))

    # Extract acidic side chain atoms (D: OD1/OD2, E: OE1/OE2)
    acidic_atoms = []  # (idx_in_residues_info, coord)
    for idx, (resnum, resname, atoms) in enumerate(residues_info):
        if resname == 'ASP':
            for aname in ('OD1', 'OD2'):
                if aname in atoms:
                    acidic_atoms.append((idx, atoms[aname]))
        elif resname == 'GLU':
            for aname in ('OE1', 'OE2'):
                if aname in atoms:
                    acidic_atoms.append((idx, atoms[aname]))

    # Lys-Lys contacts (NZ-NZ within 16Å)
    lys_lys_counts = defaultdict(int)
    lys_lys_pairs = []  # (idx_i, idx_j) pairs
    for i in range(len(lys_nz)):
        idx_i, coord_i = lys_nz[i]
        for j in range(i + 1, len(lys_nz)):
            idx_j, coord_j = lys_nz[j]
            dist = np.linalg.norm(coord_i - coord_j)
            if dist <= LYS_LYS_CUTOFF:
                lys_lys_counts[idx_i] += 1
                lys_lys_counts[idx_j] += 1
                lys_lys_pairs.append((idx_i, idx_j))

    # Lys-acidic contacts (NZ to D/E carboxyl within 8Å)
    lys_acidic_counts = defaultdict(int)
    lys_acidic_pairs = []  # (lys_idx, acidic_idx) pairs
    for lys_idx, lys_coord in lys_nz:
        contacted_residues = set()
        for acid_idx, acid_coord in acidic_atoms:
            if acid_idx == lys_idx:
                continue  # skip self
            dist = np.linalg.norm(lys_coord - acid_coord)
            if dist <= LYS_ACIDIC_CUTOFF:
                contacted_residues.add(acid_idx)
        for acid_idx in contacted_residues:
            lys_acidic_pairs.append((lys_idx, acid_idx))
        lys_acidic_counts[lys_idx] = len(contacted_residues)

    return lys_lys_counts, lys_acidic_counts, lys_lys_pairs, lys_acidic_pairs, n_residues


def process_construct(construct_name):
    """Process one construct: compute lysine contacts from best model per seed."""
    from Bio.PDB import MMCIFParser
    parser = MMCIFParser(QUIET=True)

    construct_dir = os.path.join(AF3_DIR, construct_name)
    model_dirs = sorted(glob.glob(os.path.join(construct_dir, "seed-*_sample-*")))

    if not model_dirs:
        return None

    domains = parse_construct_domains(construct_name)
    if not domains:
        return None

    fl_map = construct_to_fl_mapping(domains)

    all_lys_lys = defaultdict(list)
    all_lys_acidic = defaultdict(list)
    all_ll_pairs = defaultdict(int)   # (fl_i, fl_j) -> count across models
    all_la_pairs = defaultdict(int)   # (fl_lys, fl_acid) -> count across models
    n_models = 0

    for mdir in model_dirs[:25]:
        cif_path = os.path.join(mdir, "model.cif")
        if not os.path.exists(cif_path):
            continue

        try:
            structure = parser.get_structure('s', cif_path)
            ll_counts, la_counts, ll_pairs, la_pairs, n_res = \
                compute_lysine_contacts_from_structure(structure)

            if n_res != len(fl_map):
                continue

            for idx in range(n_res):
                all_lys_lys[idx].append(ll_counts.get(idx, 0))
                all_lys_acidic[idx].append(la_counts.get(idx, 0))

            # Map pairs to FL positions
            for i, j in ll_pairs:
                if i < len(fl_map) and j < len(fl_map):
                    fi, fj = int(fl_map[i]), int(fl_map[j])
                    key = (min(fi, fj), max(fi, fj))
                    all_ll_pairs[key] += 1
            for i, j in la_pairs:
                if i < len(fl_map) and j < len(fl_map):
                    fi, fj = int(fl_map[i]), int(fl_map[j])
                    all_la_pairs[(fi, fj)] += 1

            n_models += 1
        except Exception:
            continue

    if n_models == 0:
        return None

    mean_lys_lys = np.zeros(len(fl_map))
    mean_lys_acidic = np.zeros(len(fl_map))
    for idx in range(len(fl_map)):
        if idx in all_lys_lys:
            mean_lys_lys[idx] = np.mean(all_lys_lys[idx])
        if idx in all_lys_acidic:
            mean_lys_acidic[idx] = np.mean(all_lys_acidic[idx])

    # Normalize pair counts by number of models
    ll_pair_freq = {k: v / n_models for k, v in all_ll_pairs.items()}
    la_pair_freq = {k: v / n_models for k, v in all_la_pairs.items()}

    return {
        'construct': construct_name,
        'domains': domains,
        'fl_indices': fl_map,
        'lys_lys': mean_lys_lys,
        'lys_acidic': mean_lys_acidic,
        'lys_lys_pairs': ll_pair_freq,
        'lys_acidic_pairs': la_pair_freq,
        'n_models': n_models,
    }


# ═══════════════════════════════════════════════════════════════════
# AF2 reference
# ═══════════════════════════════════════════════════════════════════
def compute_af2_lysine_contacts():
    """Compute lysine contacts for the AF2 full-length reference."""
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('af2', AF2_PDB)
    ll_counts, la_counts, ll_pairs, la_pairs, n_res = \
        compute_lysine_contacts_from_structure(structure)

    lys_lys = np.zeros(n_res)
    lys_acidic = np.zeros(n_res)
    for idx in range(n_res):
        lys_lys[idx] = ll_counts.get(idx, 0)
        lys_acidic[idx] = la_counts.get(idx, 0)

    return lys_lys, lys_acidic, ll_pairs, la_pairs


# ═══════════════════════════════════════════════════════════════════
# Plotting
# ═══════════════════════════════════════════════════════════════════
def add_domain_shading(ax):
    for name, start, end, color in DOMAIN_GROUPS:
        ax.axvspan(start, end, alpha=0.08, color=color, zorder=0)
    ylo, yhi = ax.get_ylim()
    y_label = yhi - 0.03 * (yhi - ylo)
    for name, start, end, color in DOMAIN_GROUPS:
        ax.text((start + end) / 2, y_label, name,
                ha='center', va='top', fontsize=6, rotation=45, color=color,
                fontweight='bold')


def plot_contact_scatter(results, contact_type, ylabel, title, outpath,
                         ref_data=None, fold_change=False, color='purple'):
    """Plot per-residue contact counts over FL sequence."""
    fig, ax = plt.subplots(figsize=(20, 6))

    xs = []
    ys = []

    for r in results:
        fl_map = r['fl_indices']
        values = r[contact_type]

        for i, fi in enumerate(fl_map):
            val = values[i]
            if val == 0 and not fold_change:
                continue  # skip zero contacts for clarity

            if fold_change and ref_data is not None:
                ref_val = ref_data[fi]
                if ref_val > 0:
                    ys.append(val / ref_val)
                    xs.append(fi + 1)
                elif val > 0:
                    ys.append(val + 1)  # gained contact (ref was 0)
                    xs.append(fi + 1)
            else:
                ys.append(val)
                xs.append(fi + 1)

    if not xs:
        print(f"  WARNING: No data for {contact_type}")
        plt.close()
        return

    xs = np.array(xs)
    ys = np.array(ys)

    ax.scatter(xs, ys, s=0.5, alpha=0.15, c=color, rasterized=True, edgecolors='none')

    # Median line
    residue_vals = defaultdict(list)
    for x, y in zip(xs, ys):
        residue_vals[x].append(y)
    med_x = sorted(residue_vals.keys())
    med_y = [np.median(residue_vals[x]) for x in med_x]
    ax.plot(med_x, med_y, color='darkred', linewidth=0.8, alpha=0.7, label='Median')

    if fold_change:
        ax.axhline(1.0, color='red', linewidth=0.8, linestyle='--', alpha=0.5)

    ax.set_xlabel('Full-Length Residue Index', fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.set_xlim(1, N_FL)

    if fold_change:
        p95 = np.percentile(ys, 95)
        p5 = np.percentile(ys, 5)
        margin = 0.15 * (p95 - p5)
        ax.set_ylim(max(0, p5 - margin), p95 + margin)

    add_domain_shading(ax)

    legend_elements = [Line2D([0], [0], color='darkred', linewidth=1, label='Median')]
    if fold_change:
        legend_elements.append(Line2D([0], [0], color='red', linestyle='--', label='FC=1'))
    ax.legend(handles=legend_elements, loc='upper right', fontsize=8)

    plt.tight_layout()
    fig.savefig(outpath + '.png', dpi=200, bbox_inches='tight')
    fig.savefig(outpath + '.pdf', bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}.png ({len(results)} constructs, {len(xs)} points)")


def plot_2d_contact_map(results, pair_key, title, outpath, af2_pairs=None):
    """2D contact frequency map: FL residue on x and y, colored by frequency.

    Each point = a residue pair that is in contact. Color = fraction of
    constructs (that contain BOTH residues) where the contact is observed.
    """
    from matplotlib.colors import LogNorm

    # Aggregate pair frequencies across all constructs
    # pair -> (total_freq, n_constructs_with_both_residues)
    pair_stats = defaultdict(lambda: [0.0, 0])

    for r in results:
        pairs = r.get(pair_key, {})
        fl_set = set(int(x) for x in r['fl_indices'])

        for (fi, fj), freq in pairs.items():
            if fi in fl_set and fj in fl_set:
                pair_stats[(fi, fj)][0] += freq
                pair_stats[(fi, fj)][1] += 1

    if not pair_stats:
        print(f"  WARNING: No pairs for {pair_key}")
        return

    # Compute mean frequency (normalized by constructs containing both residues)
    xs, ys, cs = [], [], []
    for (fi, fj), (total_freq, n_constructs) in pair_stats.items():
        mean_freq = total_freq / n_constructs if n_constructs > 0 else 0
        # Plot both triangles (symmetric)
        xs.extend([fi + 1, fj + 1])
        ys.extend([fj + 1, fi + 1])
        cs.extend([n_constructs, n_constructs])

    xs = np.array(xs)
    ys = np.array(ys)
    cs = np.array(cs, dtype=float)

    fig, ax = plt.subplots(figsize=(12, 12))

    sc = ax.scatter(xs, ys, c=cs, s=1.5, alpha=0.6, cmap='hot_r',
                    rasterized=True, edgecolors='none',
                    norm=LogNorm(vmin=max(1, cs.min()), vmax=cs.max()))

    cbar = plt.colorbar(sc, ax=ax, shrink=0.8)
    cbar.set_label('# Constructs with Contact', fontsize=10)

    # Domain grid lines and labels
    for name, start, end, color in DOMAIN_GROUPS:
        ax.axvline(start, color=color, linewidth=0.3, alpha=0.5)
        ax.axhline(start, color=color, linewidth=0.3, alpha=0.5)
        # Label along diagonal
        mid = (start + end) / 2
        ax.text(mid, mid, name, ha='center', va='center', fontsize=6,
                color=color, fontweight='bold', alpha=0.7,
                bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.7, ec='none'))

    # AF2 reference pairs as hollow circles
    if af2_pairs:
        af2_x, af2_y = [], []
        for (i, j) in af2_pairs:
            af2_x.extend([i + 1, j + 1])
            af2_y.extend([j + 1, i + 1])
        ax.scatter(af2_x, af2_y, s=8, facecolors='none', edgecolors='blue',
                   linewidths=0.5, alpha=0.3, label='AF2 FL contact', zorder=0)

    ax.set_xlabel('Full-Length Residue Index', fontsize=12)
    ax.set_ylabel('Full-Length Residue Index', fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.set_xlim(1, N_FL)
    ax.set_ylim(1, N_FL)
    ax.set_aspect('equal')
    ax.invert_yaxis()

    if af2_pairs:
        ax.legend(fontsize=8, loc='lower right')

    plt.tight_layout()
    fig.savefig(outpath + '.png', dpi=200, bbox_inches='tight')
    fig.savefig(outpath + '.pdf', bbox_inches='tight')
    plt.close()

    n_unique = len(pair_stats)
    print(f"  Saved: {outpath}.png ({n_unique} unique pairs, {len(results)} constructs)")


# ═══════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════
def main():
    parser = argparse.ArgumentParser(description='Lysine contact analysis for PARP14 constructs')
    parser.add_argument('--workers', type=int, default=16)
    parser.add_argument('--constructs', nargs='+', help='Specific constructs to analyze')
    parser.add_argument('--fold-change', action='store_true',
                        help='Plot fold change relative to AF2 FL')
    parser.add_argument('--max-models', type=int, default=5,
                        help='Max models per construct to average')
    args = parser.parse_args()

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(CACHE_DIR, exist_ok=True)

    # Find constructs
    if args.constructs:
        all_constructs = [c.lower() for c in args.constructs]
    else:
        all_constructs = []
        for d in sorted(os.listdir(AF3_DIR)):
            full = os.path.join(AF3_DIR, d)
            if not os.path.isdir(full) or d in ('logs', 'test_run', 'missing_l1'):
                continue
            if glob.glob(os.path.join(full, "seed-*_sample-*/model.cif")):
                all_constructs.append(d)

    print(f"Found {len(all_constructs)} constructs")

    # Check cache
    cache_file = os.path.join(CACHE_DIR, 'lysine_contacts.npz')
    results = []
    remaining = all_constructs

    if os.path.exists(cache_file):
        cached = np.load(cache_file, allow_pickle=True)
        results = list(cached['results'])
        cached_names = {r['construct'] for r in results}
        remaining = [c for c in all_constructs if c not in cached_names]
        print(f"  Loaded {len(results)} from cache, {len(remaining)} remaining")

    if remaining:
        print(f"Processing {len(remaining)} constructs ({args.workers} workers)...")
        done = 0
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(process_construct, c): c for c in remaining}
            for future in as_completed(futures):
                done += 1
                result = future.result()
                if result is not None:
                    results.append(result)
                if done % 50 == 0:
                    print(f"  {done}/{len(remaining)} processed ({len(results)} total)")

        # Save cache
        np.savez(cache_file, results=np.array(results, dtype=object))
        print(f"  Cached {len(results)} results")

    print(f"\nTotal constructs with data: {len(results)}")

    # AF2 reference — clear old cache format (no pairs) and recompute
    print("Computing AF2 reference lysine contacts...")
    af2_cache = os.path.join(CACHE_DIR, 'af2_lysine_contacts.npz')
    af2_ll_pairs = af2_la_pairs = None
    if os.path.exists(af2_cache):
        d = np.load(af2_cache, allow_pickle=True)
        if 'lys_lys_pairs' in d:
            af2_ll = d['lys_lys']
            af2_la = d['lys_acidic']
            af2_ll_pairs = list(d['lys_lys_pairs'])
            af2_la_pairs = list(d['lys_acidic_pairs'])
        else:
            os.remove(af2_cache)  # old format, recompute

    if af2_ll_pairs is None:
        af2_ll, af2_la, af2_ll_pairs, af2_la_pairs = compute_af2_lysine_contacts()
        np.savez(af2_cache, lys_lys=af2_ll, lys_acidic=af2_la,
                 lys_lys_pairs=af2_ll_pairs, lys_acidic_pairs=af2_la_pairs)

    n_af2_lys = np.sum(af2_ll > 0)
    n_af2_la = np.sum(af2_la > 0)
    print(f"  AF2: {n_af2_lys} Lys with Lys-Lys contacts, {n_af2_la} Lys with acidic contacts")

    # Plot: raw contacts
    print("\nPlotting Lys-Lys contacts (NZ-NZ ≤ 16Å)...")
    plot_contact_scatter(
        results, 'lys_lys',
        f'Lys-Lys Contact Count (NZ-NZ ≤ {LYS_LYS_CUTOFF:.0f}Å)',
        'Per-Lysine Contact Count (NZ-NZ) Across Constructs',
        os.path.join(OUTPUT_DIR, 'lys_lys_contacts'),
        color='purple',
    )

    print("Plotting Lys-acidic contacts (NZ to D/E ≤ 8Å)...")
    plot_contact_scatter(
        results, 'lys_acidic',
        f'Lys-Acidic Contact Count (NZ to D/E ≤ {LYS_ACIDIC_CUTOFF:.0f}Å)',
        'Per-Lysine Acidic Contact Count Across Constructs',
        os.path.join(OUTPUT_DIR, 'lys_acidic_contacts'),
        color='teal',
    )

    # Plot: fold change
    if args.fold_change:
        print("\nPlotting fold changes...")
        plot_contact_scatter(
            results, 'lys_lys',
            'Lys-Lys Contact FC (construct / FL)',
            'Per-Lysine Contact Fold Change (NZ-NZ ≤ 16Å)',
            os.path.join(OUTPUT_DIR, 'lys_lys_contacts_fc'),
            ref_data=af2_ll, fold_change=True, color='purple',
        )
        plot_contact_scatter(
            results, 'lys_acidic',
            'Lys-Acidic Contact FC (construct / FL)',
            'Per-Lysine Acidic Contact Fold Change (NZ to D/E ≤ 8Å)',
            os.path.join(OUTPUT_DIR, 'lys_acidic_contacts_fc'),
            ref_data=af2_la, fold_change=True, color='teal',
        )

    # 2D contact maps
    print("\nPlotting 2D Lys-Lys contact map...")
    plot_2d_contact_map(
        results, 'lys_lys_pairs',
        f'Lys-Lys Contact Map (NZ-NZ ≤ {LYS_LYS_CUTOFF:.0f}Å, colored by frequency)',
        os.path.join(OUTPUT_DIR, 'lys_lys_contact_map_2d'),
        af2_pairs=af2_ll_pairs,
    )

    print("Plotting 2D Lys-acidic contact map...")
    plot_2d_contact_map(
        results, 'lys_acidic_pairs',
        f'Lys-Acidic Contact Map (NZ to D/E ≤ {LYS_ACIDIC_CUTOFF:.0f}Å, colored by frequency)',
        os.path.join(OUTPUT_DIR, 'lys_acidic_contact_map_2d'),
        af2_pairs=af2_la_pairs,
    )

    # Summary CSV
    csv_path = os.path.join(OUTPUT_DIR, 'lysine_contact_summary.csv')
    with open(csv_path, 'w') as f:
        f.write('construct,fl_residue,resname,lys_lys_count,lys_acidic_count\n')
        for r in results:
            fl_map = r['fl_indices']
            for i, fi in enumerate(fl_map):
                ll = r['lys_lys'][i]
                la = r['lys_acidic'][i]
                if ll > 0 or la > 0:
                    f.write(f"{r['construct']},{fi+1},LYS,{ll:.1f},{la:.1f}\n")
    print(f"\nSaved summary: {csv_path}")
    print("Done!")


if __name__ == '__main__':
    main()
