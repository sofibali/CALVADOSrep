#!/usr/bin/env python3
"""
Compute fold-change metrics (pLDDT, disorder propensity, SASA) for each AF3
domain-deletion construct relative to the AF2 full-length PARP14 reference.

Outputs:
  - Per-residue pLDDT fold change (AF3 / AF2) scatter plot
  - Per-residue disorder propensity fold change scatter plot
  - Per-domain SASA fold change scatter plot

Each dot = one AF3 construct (mean across 25 models).
x-axis = FL residue index or domain index.
y-axis = fold change relative to AF2 full-length.
"""
import os
import sys
import csv
import glob
import warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

warnings.filterwarnings('ignore')

# ═══════════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════════
AF2_PDB = "/home/sbali/CALVADOS/examples/PARP14_MDP/input/parp14.pdb"
AF3_DIR = "/home/sbali/CALVADOS/parp14/alphafold_outputs"
DOMAIN_CSV = "/home/sbali/CALVADOS/parp14/input/domain_boundaries.csv"
FASTA_FILE = "/home/sbali/CALVADOS/parp14/input/PARP14.fasta"
OUTPUT_DIR = "/home/sbali/CALVADOS/parp14/analysis/fold_change_plots"
CACHE_DIR = "/home/sbali/CALVADOS/parp14/analysis/fold_change_cache"

MAX_WORKERS = 8  # parallel workers for processing constructs
MAX_MODELS = 5   # max models per construct to average (use best-ranked seeds)

DOMAIN_ORDER = ['RRM1', 'RRM2', 'RRM3', 'KH1-KH6', 'KH7a', 'MD1L1',
                'MD2', 'MD3', 'KHb-KH8', 'WWE', 'ART']

DOMAIN_COLORS = {
    'RRM1': '#1f77b4', 'RRM2': '#2ca02c', 'RRM3': '#aec7e8',
    'KH1-KH6': '#ff7f0e', 'KH7a': '#d62728', 'MD1L1': '#aec7e8',
    'MD2': '#9467bd', 'MD3': '#8c564b', 'KHb-KH8': '#e377c2',
    'WWE': '#7f7f7f', 'ART': '#bcbd22',
}

# ═══════════════════════════════════════════════════════════════════
# Domain boundaries
# ═══════════════════════════════════════════════════════════════════
def load_domain_boundaries():
    """Load FL domain boundaries from CSV. Returns {name: (start, end)} 1-indexed inclusive."""
    bounds = {}
    order = []
    with open(DOMAIN_CSV) as f:
        for row in csv.DictReader(f):
            name = row['Domain'].strip()
            bounds[name] = (int(row['Start']), int(row['End']))
            order.append(name)
    return bounds, order


def parse_construct_domains(construct_name, domain_bounds):
    """Parse construct directory name to get ordered list of domains present."""
    # Construct names are lowercase, underscore-separated domain names
    name_lower = construct_name.lower()
    domains_present = []
    for dname in DOMAIN_ORDER:
        dname_lower = dname.lower().replace('-', '-')
        if dname_lower in name_lower.split('_') or dname_lower in name_lower:
            # More careful matching
            pass

    # Better approach: try all domain names and check if they appear as tokens
    # Domain names in dir: rrm1, rrm2, rrm3, kh1-kh6, kh7a, md1l1, md2, md3, khb-kh8, wwe, art
    domain_tokens = {
        'rrm1': 'RRM1', 'rrm2': 'RRM2', 'rrm3': 'RRM3',
        'kh1-kh6': 'KH1-KH6', 'kh7a': 'KH7a', 'md1l1': 'MD1L1',
        'md2': 'MD2', 'md3': 'MD3', 'khb-kh8': 'KHb-KH8',
        'wwe': 'WWE', 'art': 'ART',
    }

    # Split by underscore but preserve hyphenated names
    parts = name_lower.replace('kh1-kh6', 'KH1X6').replace('khb-kh8', 'KHBX8')
    parts = parts.split('_')
    parts = [p.replace('KH1X6', 'kh1-kh6').replace('KHBX8', 'khb-kh8') for p in parts]

    domains_present = []
    for p in parts:
        if p in domain_tokens:
            domains_present.append(domain_tokens[p])

    return domains_present


def load_fl_sequence():
    """Load the full-length PARP14 sequence."""
    with open(FASTA_FILE) as f:
        lines = f.readlines()
    return ''.join(l.strip() for l in lines[1:])


# Cache FL sequence at module level
_FL_SEQ = None
def get_fl_sequence():
    global _FL_SEQ
    if _FL_SEQ is None:
        _FL_SEQ = load_fl_sequence()
    return _FL_SEQ


def construct_to_fl_mapping_from_cif(cif_path, domain_bounds):
    """
    Map construct residues to FL positions by matching the construct's
    amino acid sequence against the FL sequence. Handles overlaps and gaps
    correctly since it uses the actual construct sequence from the CIF.
    Returns array of FL 0-based indices, one per construct residue.
    """
    from Bio.PDB import MMCIFParser
    from Bio.Data.IUPACData import protein_letters_3to1

    fl_seq = get_fl_sequence()

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('s', cif_path)
    construct_seq = ''
    for residue in structure[0].get_residues():
        if residue.id[0] != ' ':
            continue
        rn = residue.resname.strip().upper()
        construct_seq += protein_letters_3to1.get(rn.capitalize(), 'X')

    # Find the construct sequence as subsequence(s) of FL
    # The construct is a concatenation of domain segments from FL
    # Use greedy matching: scan FL left-to-right, consume construct chars
    fl_indices = []
    ci = 0  # construct index
    fi = 0  # fl index
    while ci < len(construct_seq) and fi < len(fl_seq):
        if construct_seq[ci] == fl_seq[fi]:
            fl_indices.append(fi)
            ci += 1
            fi += 1
        else:
            fi += 1

    if ci < len(construct_seq):
        # Fallback: couldn't match all residues
        return None

    return np.array(fl_indices)


def construct_to_fl_mapping(domains_present, domain_bounds):
    """
    Map construct residue indices (0-based) to FL residue indices (0-based).
    Handles the MD1L1/MD2 overlap at residue 1004 by deduplicating.
    Returns array of length n_construct_residues with FL 0-based positions.
    """
    fl_indices = []
    prev_end = -1  # track last FL index added to avoid overlap
    for dname in domains_present:
        start, end = domain_bounds[dname]
        for fi in range(start - 1, end):  # 0-based
            if fi > prev_end:
                fl_indices.append(fi)
                prev_end = fi

    return np.array(fl_indices)


# ═══════════════════════════════════════════════════════════════════
# Extract metrics from structures
# ═══════════════════════════════════════════════════════════════════
def extract_plddt_from_pdb(pdb_path):
    """Extract per-residue pLDDT (CA B-factor) from PDB."""
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('s', pdb_path)
    plddts = []
    for residue in structure[0].get_residues():
        if residue.id[0] != ' ':
            continue
        cas = [a for a in residue if a.name == 'CA']
        if cas:
            plddts.append(cas[0].bfactor)
    return np.array(plddts)


def extract_plddt_from_cif(cif_path):
    """Extract per-residue pLDDT (CA B-factor) from mmCIF."""
    from Bio.PDB import MMCIFParser
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('s', cif_path)
    plddts = []
    for residue in structure[0].get_residues():
        if residue.id[0] != ' ':
            continue
        cas = [a for a in residue if a.name == 'CA']
        if cas:
            plddts.append(cas[0].bfactor)
    return np.array(plddts)


def compute_sasa_from_pdb(pdb_path):
    """Compute per-residue SASA using mdtraj."""
    import mdtraj as md
    traj = md.load(pdb_path)
    sasa = md.shrake_rupley(traj, mode='residue')[0]  # nm^2
    return sasa * 100  # convert to Å^2


def compute_sasa_from_cif(cif_path):
    """Compute per-residue SASA from CIF via temporary PDB conversion."""
    import mdtraj as md
    import tempfile
    from Bio.PDB import MMCIFParser, PDBIO

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('s', cif_path)

    with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
        tmp_path = tmp.name
        io = PDBIO()
        io.set_structure(structure)
        io.save(tmp_path)

    try:
        traj = md.load(tmp_path)
        sasa = md.shrake_rupley(traj, mode='residue')[0]  # nm^2
        return sasa * 100  # Å^2
    finally:
        os.unlink(tmp_path)


def disorder_propensity(plddt_array):
    """
    Compute per-residue disorder propensity from pLDDT.
    Uses sigmoid mapping: high pLDDT -> low disorder, low pLDDT -> high disorder.
    propensity = 1 - pLDDT/100 (simple linear mapping, 0-1 scale)
    """
    return 1.0 - np.clip(plddt_array, 0, 100) / 100.0


# ═══════════════════════════════════════════════════════════════════
# Process one construct
# ═══════════════════════════════════════════════════════════════════
def process_construct(construct_name, domain_bounds):
    """
    Process a single AF3 construct: extract pLDDT, SASA, map to FL positions.
    Returns dict with FL-indexed arrays (mean across models).
    """
    construct_dir = os.path.join(AF3_DIR, construct_name)
    model_dirs = sorted(glob.glob(os.path.join(construct_dir, "seed-*_sample-*")))

    if not model_dirs:
        return None

    domains_present = parse_construct_domains(construct_name, domain_bounds)
    if not domains_present:
        return None

    fl_map = construct_to_fl_mapping(domains_present, domain_bounds)

    plddts_all = []
    sasas_all = []

    for mdir in model_dirs[:MAX_MODELS]:
        cif_path = os.path.join(mdir, "model.cif")
        if not os.path.exists(cif_path):
            continue

        try:
            plddt = extract_plddt_from_cif(cif_path)
            sasa = compute_sasa_from_cif(cif_path)

            if len(plddt) != len(fl_map):
                continue
            if len(sasa) != len(fl_map):
                continue

            plddts_all.append(plddt)
            sasas_all.append(sasa)
        except Exception:
            continue

    if not plddts_all:
        return None

    mean_plddt = np.mean(plddts_all, axis=0)
    mean_sasa = np.mean(sasas_all, axis=0)

    return {
        'construct': construct_name,
        'domains': domains_present,
        'fl_indices': fl_map,
        'plddt': mean_plddt,
        'sasa': mean_sasa,
        'n_models': len(plddts_all),
    }


# ═══════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(CACHE_DIR, exist_ok=True)

    domain_bounds, domain_order = load_domain_boundaries()
    n_fl = 1801

    # ── AF2 reference ──
    print("Computing AF2 reference metrics...")
    cache_af2 = os.path.join(CACHE_DIR, "af2_reference.npz")
    if os.path.exists(cache_af2):
        data = np.load(cache_af2)
        af2_plddt = data['plddt']
        af2_sasa = data['sasa']
        print(f"  Loaded from cache ({len(af2_plddt)} residues)")
    else:
        af2_plddt = extract_plddt_from_pdb(AF2_PDB)
        af2_sasa = compute_sasa_from_pdb(AF2_PDB)
        np.savez(cache_af2, plddt=af2_plddt, sasa=af2_sasa)
        print(f"  Computed: {len(af2_plddt)} residues, mean pLDDT={af2_plddt.mean():.1f}")

    af2_disorder = disorder_propensity(af2_plddt)

    # ── Find all AF3 constructs with completed models ──
    all_constructs = []
    for d in sorted(os.listdir(AF3_DIR)):
        full = os.path.join(AF3_DIR, d)
        if not os.path.isdir(full) or d in ('logs', 'test_run', 'missing_l1'):
            continue
        if glob.glob(os.path.join(full, "seed-*_sample-*/model.cif")):
            all_constructs.append(d)

    print(f"\nFound {len(all_constructs)} AF3 constructs with models")

    # ── Process constructs (parallel) ──
    print(f"Processing constructs ({MAX_WORKERS} workers, {MAX_MODELS} models each)...")

    # Check cache
    cache_file = os.path.join(CACHE_DIR, "construct_metrics.npz")
    results = []

    if os.path.exists(cache_file):
        cached = np.load(cache_file, allow_pickle=True)
        results = list(cached['results'])
        cached_names = {r['construct'] for r in results}
        remaining = [c for c in all_constructs if c not in cached_names]
        print(f"  Loaded {len(results)} from cache, {len(remaining)} remaining")
    else:
        remaining = all_constructs

    if remaining:
        done = 0
        with ProcessPoolExecutor(max_workers=MAX_WORKERS) as pool:
            futures = {pool.submit(process_construct, c, domain_bounds): c
                       for c in remaining}
            for future in as_completed(futures):
                done += 1
                result = future.result()
                if result is not None:
                    results.append(result)
                if done % 50 == 0:
                    print(f"  {done}/{len(remaining)} processed ({len(results)} successful)")

        # Save cache
        np.savez(cache_file, results=np.array(results, dtype=object))
        print(f"  Cached {len(results)} results")

    print(f"\nTotal constructs with valid data: {len(results)}")

    # ═══════════════════════════════════════════════════════════════
    # Compute fold changes and build plot data
    # ═══════════════════════════════════════════════════════════════

    # Per-residue: collect (fl_residue_idx, fold_change) for each construct
    plddt_fc_points = defaultdict(list)     # fl_idx -> [fc1, fc2, ...]
    disorder_fc_points = defaultdict(list)
    # Per-domain: collect (domain_name, fold_change_sasa) for each construct
    sasa_fc_per_domain = defaultdict(list)  # domain_name -> [fc1, fc2, ...]

    for r in results:
        fl_idx = r['fl_indices']
        plddt = r['plddt']
        sasa = r['sasa']
        disorder = disorder_propensity(plddt)
        domains = r['domains']

        # Per-residue fold change (AF3 / AF2)
        for i, fi in enumerate(fl_idx):
            if af2_plddt[fi] > 0:
                plddt_fc_points[fi].append(plddt[i] / af2_plddt[fi])
            if af2_disorder[fi] > 0.01:
                disorder_fc_points[fi].append(disorder[i] / af2_disorder[fi])

        # Per-domain SASA fold change (use FL mapping to group residues by domain)
        for dname in domains:
            dstart, dend = domain_bounds[dname]
            # Find construct residues that map to this domain's FL range
            mask = (fl_idx >= dstart - 1) & (fl_idx < dend)
            if mask.sum() == 0:
                continue

            domain_sasa_af3 = sasa[mask].sum()
            domain_sasa_af2 = af2_sasa[dstart - 1:dend].sum()

            if domain_sasa_af2 > 0:
                sasa_fc_per_domain[dname].append(domain_sasa_af3 / domain_sasa_af2)

    # ═══════════════════════════════════════════════════════════════
    # Plot 1: Per-residue pLDDT fold change
    # ═══════════════════════════════════════════════════════════════
    print("\nPlotting pLDDT fold change...")
    fig, ax = plt.subplots(figsize=(20, 6))

    for fi in sorted(plddt_fc_points.keys()):
        fcs = plddt_fc_points[fi]
        xs = np.full(len(fcs), fi + 1)  # 1-indexed
        ax.scatter(xs, fcs, s=0.3, alpha=0.15, c='steelblue', rasterized=True)

    # Domain shading
    for dname in DOMAIN_ORDER:
        s, e = domain_bounds[dname]
        color = DOMAIN_COLORS.get(dname, '#cccccc')
        ax.axvspan(s, e, alpha=0.08, color=color)
        ax.text((s + e) / 2, ax.get_ylim()[1] if ax.get_ylim()[1] > 2 else 2.0,
                dname, ha='center', va='bottom', fontsize=7, rotation=45)

    ax.axhline(1.0, color='red', linewidth=0.8, linestyle='--', label='No change')
    ax.set_xlabel('Full-Length Residue Index', fontsize=12)
    ax.set_ylabel('pLDDT Fold Change (AF3 / AF2)', fontsize=12)
    ax.set_title('Per-Residue pLDDT Fold Change Across Domain Deletion Constructs', fontsize=14)
    ax.set_xlim(1, n_fl)
    ax.set_ylim(0, max(3.0, np.percentile([v for vs in plddt_fc_points.values() for v in vs], 99)))

    # Add domain shading to legend
    legend_elements = [Line2D([0], [0], color='red', linestyle='--', label='No change (FC=1)')]
    for dname in DOMAIN_ORDER:
        legend_elements.append(
            plt.Rectangle((0, 0), 1, 1, fc=DOMAIN_COLORS.get(dname, '#ccc'), alpha=0.3, label=dname))
    ax.legend(handles=legend_elements, loc='upper right', fontsize=7, ncol=3)

    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'plddt_fold_change_per_residue.png'), dpi=200)
    fig.savefig(os.path.join(OUTPUT_DIR, 'plddt_fold_change_per_residue.pdf'))
    plt.close()
    print(f"  Saved to {OUTPUT_DIR}/plddt_fold_change_per_residue.png")

    # ═══════════════════════════════════════════════════════════════
    # Plot 2: Per-residue disorder propensity fold change
    # ═══════════════════════════════════════════════════════════════
    print("Plotting disorder propensity fold change...")
    fig, ax = plt.subplots(figsize=(20, 6))

    for fi in sorted(disorder_fc_points.keys()):
        fcs = disorder_fc_points[fi]
        xs = np.full(len(fcs), fi + 1)
        ax.scatter(xs, fcs, s=0.3, alpha=0.15, c='darkorange', rasterized=True)

    for dname in DOMAIN_ORDER:
        s, e = domain_bounds[dname]
        color = DOMAIN_COLORS.get(dname, '#cccccc')
        ax.axvspan(s, e, alpha=0.08, color=color)
        ax.text((s + e) / 2, ax.get_ylim()[1] if ax.get_ylim()[1] > 2 else 2.0,
                dname, ha='center', va='bottom', fontsize=7, rotation=45)

    ax.axhline(1.0, color='red', linewidth=0.8, linestyle='--', label='No change')
    ax.set_xlabel('Full-Length Residue Index', fontsize=12)
    ax.set_ylabel('Disorder Propensity Fold Change (AF3 / AF2)', fontsize=12)
    ax.set_title('Per-Residue Disorder Propensity Fold Change Across Domain Deletion Constructs', fontsize=14)
    ax.set_xlim(1, n_fl)
    ymax = np.percentile([v for vs in disorder_fc_points.values() for v in vs], 99) if disorder_fc_points else 3
    ax.set_ylim(0, max(3.0, ymax))

    ax.legend(handles=legend_elements, loc='upper right', fontsize=7, ncol=3)

    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'disorder_fold_change_per_residue.png'), dpi=200)
    fig.savefig(os.path.join(OUTPUT_DIR, 'disorder_fold_change_per_residue.pdf'))
    plt.close()
    print(f"  Saved to {OUTPUT_DIR}/disorder_fold_change_per_residue.png")

    # ═══════════════════════════════════════════════════════════════
    # Plot 3: Per-domain SASA fold change
    # ═══════════════════════════════════════════════════════════════
    print("Plotting SASA fold change per domain...")
    fig, ax = plt.subplots(figsize=(14, 6))

    domain_x_positions = {d: i for i, d in enumerate(DOMAIN_ORDER)}
    for dname in DOMAIN_ORDER:
        if dname in sasa_fc_per_domain:
            fcs = sasa_fc_per_domain[dname]
            x = domain_x_positions[dname]
            jitter = np.random.default_rng(42).uniform(-0.3, 0.3, len(fcs))
            color = DOMAIN_COLORS.get(dname, '#666666')
            ax.scatter(np.full(len(fcs), x) + jitter, fcs,
                       s=8, alpha=0.4, c=color, edgecolors='none', rasterized=True)
            # Box plot overlay
            bp = ax.boxplot([fcs], positions=[x], widths=0.5, showfliers=False,
                            patch_artist=True, zorder=3)
            bp['boxes'][0].set_facecolor(color)
            bp['boxes'][0].set_alpha(0.3)
            bp['medians'][0].set_color('black')

    ax.axhline(1.0, color='red', linewidth=0.8, linestyle='--', label='No change (FC=1)')
    ax.set_xticks(range(len(DOMAIN_ORDER)))
    ax.set_xticklabels(DOMAIN_ORDER, rotation=45, ha='right', fontsize=10)
    ax.set_ylabel('SASA Fold Change (AF3 / AF2)', fontsize=12)
    ax.set_title('Per-Domain SASA Fold Change Across Domain Deletion Constructs', fontsize=14)
    ax.legend(fontsize=10)

    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'sasa_fold_change_per_domain.png'), dpi=200)
    fig.savefig(os.path.join(OUTPUT_DIR, 'sasa_fold_change_per_domain.pdf'))
    plt.close()
    print(f"  Saved to {OUTPUT_DIR}/sasa_fold_change_per_domain.png")

    # ═══════════════════════════════════════════════════════════════
    # Summary stats
    # ═══════════════════════════════════════════════════════════════
    print("\n── Summary ──")
    print(f"Constructs analyzed: {len(results)}")
    print(f"Residues with pLDDT FC data: {len(plddt_fc_points)}")
    print(f"\nPer-domain SASA fold change (median ± IQR):")
    for dname in DOMAIN_ORDER:
        if dname in sasa_fc_per_domain:
            fcs = np.array(sasa_fc_per_domain[dname])
            print(f"  {dname:10s}: {np.median(fcs):.3f} ({np.percentile(fcs,25):.3f} - {np.percentile(fcs,75):.3f})  n={len(fcs)}")

    # Save fold change data as CSV for downstream use
    csv_path = os.path.join(OUTPUT_DIR, 'sasa_fold_change_per_domain.csv')
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['domain', 'construct', 'sasa_fold_change'])
        for r in results:
            fl_idx = r['fl_indices']
            for dname in r['domains']:
                dstart, dend = domain_bounds[dname]
                mask = (fl_idx >= dstart - 1) & (fl_idx < dend)
                if mask.sum() == 0:
                    continue
                sasa_af3 = r['sasa'][mask].sum()
                sasa_af2 = af2_sasa[dstart - 1:dend].sum()
                if sasa_af2 > 0:
                    w.writerow([dname, r['construct'], sasa_af3 / sasa_af2])
    print(f"\nSaved SASA fold change CSV: {csv_path}")

    csv_path2 = os.path.join(OUTPUT_DIR, 'plddt_fold_change_summary.csv')
    with open(csv_path2, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['fl_residue', 'n_constructs', 'mean_fc', 'median_fc', 'std_fc'])
        for fi in sorted(plddt_fc_points.keys()):
            fcs = np.array(plddt_fc_points[fi])
            w.writerow([fi + 1, len(fcs), f"{fcs.mean():.4f}", f"{np.median(fcs):.4f}", f"{fcs.std():.4f}"])
    print(f"Saved pLDDT fold change summary: {csv_path2}")


if __name__ == '__main__':
    main()
