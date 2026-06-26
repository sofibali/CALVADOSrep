#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Per-residue exposed Lys-Acidic contact persistence — reproduces the figure
at figures/_archive/pre_2026-05-19/lys_contacts/per_residue_contact_persistence.svg
but restricted to fl_optimized and md (MD1L1+MD2+MD3 macrodomains).

For each residue:
  - Red bar  = sum of Lys-Acidic contact persistence with all EXPOSED D/E
               (if the residue is Lys AND itself exposed)
  - Blue bar = sum of Lys-Acidic contact persistence with all EXPOSED Lys
               (if the residue is Asp/Glu AND itself exposed)

Persistence = (frames in contact) / (total frames), averaged across all 25
replicates per set. Contact = CA-CA <= 3.0 nm (coarse-grained cutoff).
Exposed = RSA >= 0.20 (Chothia-normalized, surface fraction) from the FL
AF2 reference structure.

Usage:
    python figure_lys_exposed_persistence.py
    python figure_lys_exposed_persistence.py --sets fl_optimized md
    python figure_lys_exposed_persistence.py --force  # recompute contact data
"""
import os
import sys
import argparse
import csv
import numpy as np
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

CWD = Path(__file__).resolve().parent
sys.path.insert(0, str(CWD))
from _fig_layout import get_fig_dir
import sim_registry as reg

# Maps an externally-registered set_key -> absolute folder Path. Populated by
# register_sim_folder() in the main process before any worker pool is created.
EXTERNAL_DIRS = {}

DATA_PATH = CWD / 'data'
SKIP_FRAMES = 50
LYS_ACIDIC_CUTOFF = 3.0   # nm CA-CA
RSA_EXPOSED_MIN = 0.20    # surface threshold

SEEDS = range(1, 6)
SAMPLES = range(0, 5)

DOMAIN_GROUPS = [
    ('RRM1',    1,    145, '#1f77b4'),
    ('RRM2',    146,  224, '#2ca02c'),
    ('RRM3',    225,  314, '#aec7e8'),
    ('KH1-KH6', 315,  737, '#ff7f0e'),
    ('KH7a',    738,  789, '#d62728'),
    ('MD1L1',   790,  1004, '#17becf'),
    ('MD2',     1005, 1193, '#9467bd'),
    ('MD3',     1207, 1388, '#8c564b'),
    ('KHb-KH8', 1389, 1533, '#e377c2'),
    ('WWE',     1534, 1602, '#7f7f7f'),
    ('ART',     1603, 1801, '#bcbd22'),
]

# Mapping for construct → FL residue
DOMAIN_UNITS = {
    'fl':            None,  # identity
    'fl_optimized':  None,  # identity
    'md':            ['md1l1', 'md2', 'md3'],
    'mka':           ['md1l1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art'],
    'core':          ['kh7a', 'md1l1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art'],
    'norrm':         ['kh1-kh6', 'kh7a', 'md1l1', 'md2', 'md3',
                      'khb-kh8', 'wwe', 'art'],
    'noart':         ['kh1-kh6', 'kh7a', 'md1l1', 'md2', 'md3',
                      'khb-kh8', 'wwe'],
}
UNIT_FL = {
    'rrm1':    (1, 145),    'rrm2':    (146, 224), 'rrm3':    (225, 314),
    'kh1-kh6': (315, 737),  'kh7a':    (738, 789),
    'md1l1':   (790, 1004), 'md2':     (1004, 1193), 'md3':   (1207, 1388),
    'khb-kh8': (1389, 1533), 'wwe':    (1534, 1602), 'art':   (1603, 1801),
}


def register_sim_folder(path, units=None):
    """Register an arbitrary simulation folder so it can be plotted like a named set.

    Resolves the construct's FL domain units (via metadata.json / --units / known
    basename) and injects them into DOMAIN_UNITS so construct_to_fl_map() works.
    Returns the set_key (folder basename). The folder is expected to contain
    seed-{1-5}_sample-{0-4}/ replicate dirs with top.pdb + <sysname>.dcd.
    """
    folder = Path(path).resolve()
    resolved_units = reg._resolve_units(str(folder), units=units)
    set_key = folder.name
    DOMAIN_UNITS[set_key] = list(resolved_units)
    EXTERNAL_DIRS[set_key] = folder
    return set_key


def construct_to_fl_map(set_key):
    """Return dict {construct_resid: FL_resid}."""
    units = DOMAIN_UNITS.get(set_key)
    if units is None:  # fl/fl_optimized → identity
        return {r: r for r in range(1, 1802)}
    ranges = sorted([UNIT_FL[u] for u in units if u in UNIT_FL])
    merged = [list(ranges[0])]
    for s, e in ranges[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    mapping = {}
    pos = 1
    for fl_s, fl_e in merged:
        for fl_r in range(fl_s, fl_e + 1):
            mapping[pos] = fl_r
            pos += 1
    return mapping


def load_exposed_fl():
    """Load exposed residues (RSA >= 0.20) from FL AF2 SASA CSV."""
    csv_path = DATA_PATH / 'sasa_face_per_residue.csv'
    if not csv_path.exists():
        print(f"WARNING: {csv_path} not found. "
              f"Run figure_sasa_faces.py first.")
        return None
    exposed = set()
    with open(csv_path) as f:
        rdr = csv.DictReader(f)
        for row in rdr:
            try:
                rsa = float(row['rsa']) if row.get('rsa') else None
            except ValueError:
                rsa = None
            if rsa is not None and rsa >= RSA_EXPOSED_MIN:
                exposed.add(int(row['resid']))
    return exposed


def worker_contacts(set_key, seed, sample):
    """Compute Lys-Acidic contact frequency for one replicate.

    Returns dict {(fl_resid_lys, fl_resid_acid): persistence}.
    """
    import MDAnalysis as mda

    if set_key in EXTERNAL_DIRS:
        sim_dir = EXTERNAL_DIRS[set_key] / f'seed-{seed}_sample-{sample}'
    else:
        sim_dir = CWD / set_key / f'seed-{seed}_sample-{sample}'
    pdb = sim_dir / 'top.pdb'
    dcds = list(sim_dir.glob('*.dcd'))
    if not pdb.is_file() or not dcds:
        return None

    u = mda.Universe(str(pdb), str(dcds[0]))
    resnames = u.atoms.resnames
    resids = u.atoms.resids  # 1-based, construct-local

    lys_idx = np.where(resnames == 'LYS')[0]
    asp_idx = np.where(resnames == 'ASP')[0]
    glu_idx = np.where(resnames == 'GLU')[0]
    acid_idx = np.concatenate([asp_idx, glu_idx])

    lys_resids = resids[lys_idx]
    acid_resids = resids[acid_idx]
    acid_names = resnames[acid_idx]

    fl_map = construct_to_fl_map(set_key)

    counts = defaultdict(int)
    n_frames = 0
    for ts in u.trajectory[SKIP_FRAMES:]:
        pos = u.atoms.positions / 10.0  # nm
        lys_pos = pos[lys_idx]
        acid_pos = pos[acid_idx]
        # All pairwise Lys-acid distances
        diff = lys_pos[:, None, :] - acid_pos[None, :, :]
        d = np.sqrt(np.sum(diff ** 2, axis=2))
        contacts = np.argwhere(d <= LYS_ACIDIC_CUTOFF)
        for i, j in contacts:
            r_lys = int(lys_resids[i])
            r_acid = int(acid_resids[j])
            if r_lys == r_acid:
                continue
            fl_lys = fl_map.get(r_lys, -1)
            fl_acid = fl_map.get(r_acid, -1)
            if fl_lys < 0 or fl_acid < 0:
                continue
            counts[(fl_lys, fl_acid)] += 1
        n_frames += 1

    return {k: v / n_frames for k, v in counts.items()}


def compute_per_residue_sums(set_key, exposed_fl, force=False, workers=8):
    """Return (sum_per_lys_fl, sum_per_acid_fl) over exposed partners.

    Loads cached per-replicate data if available.
    """
    cache = DATA_PATH / f'lys_acidic_contacts_{set_key}.npz'
    if cache.exists() and not force:
        d = np.load(cache, allow_pickle=True)
        per_rep_freqs = list(d['per_rep'])
    else:
        print(f"  Computing Lys-Acidic contacts for {set_key} "
              f"(25 replicates, {workers} workers)...")
        per_rep_freqs = []
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futs = {pool.submit(worker_contacts, set_key, seed, sample):
                    (seed, sample)
                    for seed in SEEDS for sample in SAMPLES}
            for fut in as_completed(futs):
                seed, sample = futs[fut]
                res = fut.result()
                if res is not None:
                    per_rep_freqs.append(res)
        np.savez(cache, per_rep=np.array(per_rep_freqs, dtype=object))
        print(f"  Cached {cache.name}")

    # Average pair frequencies across replicates
    pair_freq = defaultdict(list)
    for rep in per_rep_freqs:
        for k, v in rep.items():
            pair_freq[k].append(v)
    pair_mean = {k: float(np.mean(v)) for k, v in pair_freq.items()}

    # Per-residue sums (only counting EXPOSED partners)
    sum_lys = defaultdict(float)
    sum_acid = defaultdict(float)
    for (r_lys, r_acid), freq in pair_mean.items():
        if r_lys not in exposed_fl or r_acid not in exposed_fl:
            continue
        sum_lys[r_lys] += freq
        sum_acid[r_acid] += freq
    return dict(sum_lys), dict(sum_acid), pair_mean


def plot_per_residue(set_data, exposed_fl, outpath, set_labels):
    """Per-residue exposed Lys-Acidic contact persistence plot.

    set_data: dict {set_key: (sum_lys, sum_acid)}
    set_labels: dict {set_key: display_label}
    """
    n_sets = len(set_data)
    fig, axes = plt.subplots(n_sets, 1, figsize=(20, 3 * n_sets + 1),
                              sharex=True)
    if n_sets == 1:
        axes = [axes]

    for ax, (set_key, (sum_lys, sum_acid)) in zip(axes, set_data.items()):
        # Domain shading
        for name, start, end, color in DOMAIN_GROUPS:
            ax.axvspan(start, end, color=color, alpha=0.10, zorder=0)

        # Lys bars (red)
        if sum_lys:
            xs = sorted(sum_lys.keys())
            ys = [sum_lys[x] for x in xs]
            ax.bar(xs, ys, width=2.5, color='#e74c3c', alpha=0.9,
                   edgecolor='none', label='Lys (sum over exp D/E)')
        # Acidic bars (blue) — slightly offset to see overlaps
        if sum_acid:
            xs = sorted(sum_acid.keys())
            ys = [sum_acid[x] for x in xs]
            ax.bar(xs, ys, width=2.5, color='#3498db', alpha=0.7,
                   edgecolor='none', label='D/E (sum over exp K)')

        # Exposed residue markers along the top
        ymax = ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 1
        for r in exposed_fl:
            ax.plot([r, r], [ymax * 1.02, ymax * 1.06],
                    color='lightgray', linewidth=0.4, alpha=0.4,
                    zorder=1)

        ax.set_title(f'{set_labels.get(set_key, set_key)} (exposed only)',
                     fontsize=11, fontweight='bold')
        ax.set_ylabel('Sum persistence', fontsize=10)
        ax.legend(fontsize=8, loc='upper right')
        ax.set_xlim(1, 1801)
        ax.set_ylim(bottom=0)

    axes[-1].set_xlabel('FL residue position', fontsize=12)
    fig.suptitle('Per-residue Exposed Lys-Acidic contact persistence '
                 '(CALVADOS CG-MD)', fontsize=13)
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(f'{outpath}.png', dpi=200, bbox_inches='tight')
    fig.savefig(f'{outpath}.svg', bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}.png")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sets', nargs='+',
                        default=['fl_optimized', 'md'],
                        help='Sets to plot (default: fl_optimized md)')
    parser.add_argument('--force', action='store_true',
                        help='Recompute per-replicate contacts')
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--sim-folder', nargs='+', default=None, metavar='PATH',
                        help='One or more simulation-set folders to analyze (each '
                             'containing seed-*_sample-*/ replicates). Domain units '
                             "are read from the folder's metadata.json, or pass "
                             '--units.')
    parser.add_argument('--units', nargs='+', default=None, metavar='UNIT',
                        help='FL domain units in the --sim-folder construct (e.g. '
                             'md1l1 md2 md3). Required only if the folder has no '
                             f'metadata.json. Valid: {", ".join(reg.DOMAIN_UNITS)}')
    args = parser.parse_args()

    # Register any external folders (in the main process so forked workers inherit
    # EXTERNAL_DIRS) and fold them into the set list.
    sets = list(args.sets)
    if args.sim_folder:
        ext_keys = [register_sim_folder(f, units=args.units) for f in args.sim_folder]
        # If --sets was not explicitly overridden, plot only the external folders;
        # otherwise append them to the requested sets.
        if args.sets == parser.get_default('sets'):
            sets = ext_keys
        else:
            sets = sets + [k for k in ext_keys if k not in sets]
    args.sets = sets

    set_labels = {
        'fl':           'Full-length (1801)',
        'fl_optimized': 'FL optimized (1801)',
        'md':           'Macrodomains (586)',
        'core':         'Core (1051)',
        'mka':          'MKA (999)',
        'norrm':        'No-RRM (1474)',
        'noart':        'No-ART (1275)',
    }

    out_dir = get_fig_dir('07_lysine_contacts')
    print(f"Output: {out_dir}")

    # Load exposed residues from AF2 FL SASA data
    exposed_fl = load_exposed_fl()
    if exposed_fl is None:
        print("ERROR: cannot proceed without SASA data")
        return 1
    print(f"  Exposed residues (RSA>={RSA_EXPOSED_MIN}): "
          f"{len(exposed_fl)}/1801")

    # Compute or load per-set contact data
    set_data = {}
    for sk in args.sets:
        sum_lys, sum_acid, _ = compute_per_residue_sums(
            sk, exposed_fl, force=args.force, workers=args.workers)
        n_lys_res = sum(1 for v in sum_lys.values() if v > 0)
        n_acid_res = sum(1 for v in sum_acid.values() if v > 0)
        print(f"  {sk}: {n_lys_res} exposed Lys + {n_acid_res} "
              f"exposed D/E with contacts")
        set_data[sk] = (sum_lys, sum_acid)

    plot_per_residue(set_data, exposed_fl,
                      out_dir / 'per_residue_contact_persistence',
                      set_labels)
    return 0


if __name__ == '__main__':
    sys.exit(main())
