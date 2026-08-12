#!/usr/bin/env python3
"""
Lysine contact analysis from CALVADOS CG-MD trajectories (FL construct).

For all 25 FL replicates, computes per-frame:
  - Lys-Lys contacts: CA-CA within 1.6 nm
  - Lys-acidic contacts: Lys CA within 0.8 nm of Asp/Glu CA

Produces scatter contact maps with large dots on the 1-1801 residue grid.

Usage:
    conda run -n calvados python analyze_lys_contacts.py
"""

import os
import numpy as np
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CWD = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# Dated category subdir: figures/07_lysine_contacts/<YYYY-MM-DD>/
import sys as _sys
_sys.path.insert(0, CWD)
import os as _os_boot, sys as _sys_boot  # __ROOTBOOT__ (script lives in a subfolder; root=parent)
_sys_boot.path.insert(0, _os_boot.path.dirname(_os_boot.path.dirname(_os_boot.path.abspath(__file__))))
import sim_registry as reg
from _fig_layout import get_fig_dir as _get_fig_dir
FIG_PATH = str(_get_fig_dir('07_lysine_contacts'))
DATA_PATH = os.path.join(CWD, 'data')

LYS_LYS_CUTOFF = 3.0     # nm (CA-CA coarse-grained)
LYS_ACIDIC_CUTOFF = 3.0   # nm (CA-CA coarse-grained)
SKIP_FRAMES = 50
N_FL = 1801

SEEDS = range(1, 6)
SAMPLES = range(0, 5)

# Externally registered simulation folders (via --sim-folder). Populated in the
# main process before any ProcessPoolExecutor is created so forked workers inherit
# it (Linux fork start method).
EXTERNAL_DIRS = {}
EXTERNAL_SYSNAME = {}

# Cache: set_key -> replicate dirs discovered by _flat_replicate_dirs() below.
_FLAT_DIRS_CACHE = {}


def _flat_replicate_dirs(folder, sysname):
    """Every immediate subdirectory of `folder` containing a `<sysname>.dcd`.

    Fallback for --sim-folder layouts that aren't a seed-{1-5}_sample-{0-4}
    grid -- e.g. a handful of independent long runs named arbitrarily
    (fl_go's state-*_tica_seed-*_sample-*_fr*/ dirs)."""
    dirs = []
    for entry in sorted(os.listdir(folder)):
        d = os.path.join(folder, entry)
        if os.path.isdir(d) and os.path.isfile(os.path.join(d, f'{sysname}.dcd')):
            dirs.append(d)
    return dirs


def register_sim_folder(path, units=None):
    """Register an arbitrary simulation folder so it can be analyzed by set_key.

    The folder is expected to contain seed-{1-5}_sample-{0-4}/ replicate dirs
    with top.pdb (or checkpoint.pdb) + <sysname>.dcd. Domain units are resolved
    from `units`, then metadata.json, then a matching named set. Returns the
    set_key (the folder basename) under which it is registered.
    """
    folder = os.path.abspath(path)
    if not os.path.isdir(folder):
        raise FileNotFoundError(f"sim folder not found: {folder}")
    reg._resolve_units(folder, units)  # validate units are resolvable
    meta = reg.read_metadata(folder) or {}
    sysname = meta.get('sysname') or reg.detect_sysname(folder)
    set_key = os.path.basename(os.path.normpath(folder))
    EXTERNAL_DIRS[set_key] = folder
    EXTERNAL_SYSNAME[set_key] = sysname
    return set_key

DOMAIN_GROUPS = [
    ('RRM1', 1, 145, '#A0A0A0'), ('RRM2', 146, 224, '#B0B0B0'),
    ('RRM3', 225, 314, '#C0C0C0'), ('KH1-6', 315, 737, '#4ECDC4'),
    ('KH7a', 738, 789, '#45B7AA'), ('MD1', 790, 1004, '#6B1F7A'),
    ('MD2', 1004, 1193, '#9450A8'), ('MD3', 1207, 1388, '#B87FCC'),
    ('KHb-8', 1389, 1533, '#45B7AA'), ('WWE', 1534, 1602, '#F7DC6F'),
    ('ART', 1603, 1801, '#E74C3C'),
]


def worker_lys_contacts(seed, sample, set_key='fl'):
    """Compute Lys contacts for one replicate of a given set."""
    import MDAnalysis as mda

    if set_key in EXTERNAL_DIRS:
        base = EXTERNAL_DIRS[set_key]
        sysname = EXTERNAL_SYSNAME[set_key]
        sim_dir = os.path.join(base, f'seed-{seed}_sample-{sample}')
        if not os.path.isdir(sim_dir):
            if set_key not in _FLAT_DIRS_CACHE:
                _FLAT_DIRS_CACHE[set_key] = _flat_replicate_dirs(base, sysname)
            flat = _FLAT_DIRS_CACHE[set_key]
            idx = (seed - 1) * len(SAMPLES) + sample
            if idx < len(flat):
                sim_dir = flat[idx]
        dcd = os.path.join(sim_dir, f'{sysname}.dcd')
    else:
        sim_dir = os.path.join(CWD, set_key, f'seed-{seed}_sample-{sample}')
        dcd = os.path.join(sim_dir, 'parp14.dcd')

    pdb = os.path.join(sim_dir, 'checkpoint.pdb')
    if not os.path.isfile(pdb):
        pdb = os.path.join(sim_dir, 'top.pdb')
    if not os.path.isfile(pdb):
        pdb = os.path.join(sim_dir, 'restart.pdb')

    if not os.path.isfile(dcd):
        return None

    u = mda.Universe(pdb, dcd)
    resnames = u.atoms.resnames
    resids = u.atoms.resids  # 1-based

    lys_idx = np.where(resnames == 'LYS')[0]
    asp_idx = np.where(resnames == 'ASP')[0]
    glu_idx = np.where(resnames == 'GLU')[0]
    acidic_idx = np.concatenate([asp_idx, glu_idx])

    lys_resids = resids[lys_idx]
    acid_resids = resids[acidic_idx]

    ll_pair_counts = defaultdict(int)
    la_pair_counts = defaultdict(int)
    n_frames = 0

    for ts in u.trajectory[SKIP_FRAMES:]:
        pos = u.atoms.positions / 10.0  # Angstrom -> nm

        # Vectorized Lys-Lys distances
        lys_pos = pos[lys_idx]
        n_lys = len(lys_idx)
        for i in range(n_lys):
            dists = np.linalg.norm(lys_pos[i+1:] - lys_pos[i], axis=1)
            contacts = np.where(dists <= LYS_LYS_CUTOFF)[0]
            for ci in contacts:
                j = i + 1 + ci
                ri, rj = int(lys_resids[i]), int(lys_resids[j])
                ll_pair_counts[(min(ri, rj), max(ri, rj))] += 1

        # Vectorized Lys-acidic distances
        acid_pos = pos[acidic_idx]
        for i in range(n_lys):
            dists = np.linalg.norm(acid_pos - lys_pos[i], axis=1)
            contacts = np.where(dists <= LYS_ACIDIC_CUTOFF)[0]
            for ci in contacts:
                ri = int(lys_resids[i])
                rj = int(acid_resids[ci])
                if ri != rj:
                    la_pair_counts[(ri, rj)] += 1

        n_frames += 1

    ll_freq = {k: v / n_frames for k, v in ll_pair_counts.items()}
    la_freq = {k: v / n_frames for k, v in la_pair_counts.items()}

    return {'seed': seed, 'sample': sample, 'n_frames': n_frames,
            'll_freq': ll_freq, 'la_freq': la_freq}


def plot_contact_map(pairs, title, fname, cutoff_label):
    """Scatter plot of contacts on the FL residue grid with large dots."""
    if not pairs:
        return

    MIN_FREQ = 0.05
    xs, ys, cs = [], [], []
    for (ri, rj), freq in pairs.items():
        if freq < MIN_FREQ:
            continue
        xs.extend([ri, rj])
        ys.extend([rj, ri])
        cs.extend([freq, freq])

    xs = np.array(xs)
    ys = np.array(ys)
    cs = np.array(cs)

    if len(xs) == 0:
        print(f"  No pairs above {MIN_FREQ} for {fname}")
        return

    n_filtered = len(pairs) - len(cs) // 2
    fig, ax = plt.subplots(figsize=(14, 14))

    vmax = min(1.0, np.percentile(cs, 98)) if len(cs) > 0 else 1.0
    sc = ax.scatter(xs, ys, s=80, c=cs, cmap='YlOrRd', alpha=0.85,
                    marker='o', edgecolors='black', linewidths=0.3,
                    vmin=MIN_FREQ, vmax=vmax, zorder=3)

    cbar = plt.colorbar(sc, ax=ax, shrink=0.8)
    cbar.set_label('Contact frequency (fraction of frames)', fontsize=10)

    for name, start, end, dcolor in DOMAIN_GROUPS:
        ax.axvline(start, color=dcolor, linewidth=0.6, alpha=0.5)
        ax.axvline(end, color=dcolor, linewidth=0.6, alpha=0.5)
        ax.axhline(start, color=dcolor, linewidth=0.6, alpha=0.5)
        ax.axhline(end, color=dcolor, linewidth=0.6, alpha=0.5)
        ax.fill_between([start, end], start, end, alpha=0.04, color=dcolor, zorder=0)
        mid = (start + end) / 2
        ax.text(mid, mid, name, ha='center', va='center', fontsize=8,
                color=dcolor, fontweight='bold', alpha=0.8,
                bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.85, ec='none'),
                zorder=4)

    ax.set_xlabel('Residue Index', fontsize=13)
    ax.set_ylabel('Residue Index', fontsize=13)
    n_shown = len(cs) // 2
    ax.set_title(f'{title}\n{cutoff_label}, {n_shown} pairs shown (freq > {MIN_FREQ}), '
                 f'25 replicates', fontsize=13, fontweight='bold')
    ax.set_xlim(1, N_FL)
    ax.set_ylim(1, N_FL)
    ax.set_aspect('equal')
    ax.invert_yaxis()

    plt.tight_layout()
    for ext in ['png', 'svg']:
        fig.savefig(os.path.join(FIG_PATH, f'{fname}.{ext}'),
                    dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved: figures/{fname} ({len(pairs)} pairs)")


def _render_domain_sum_heatmap(mat, dom_names, contact_label, fname, vmax=None, subtitle=''):
    """One domain x domain contact-sum heatmap. `vmax=None` autoscales to the
    full matrix (intra-domain diagonal included); pass an explicit vmax to
    cap the color scale (e.g. at the inter-domain max) so off-diagonal
    structure isn't washed out by the much larger diagonal values."""
    n_dom = len(dom_names)
    fig, ax = plt.subplots(figsize=(10, 9))

    mat_plot = np.where(mat > 0, mat, np.nan)
    im = ax.imshow(mat_plot, cmap='YlOrRd', origin='lower', aspect='equal',
                    vmin=0, vmax=vmax)

    ax.set_xticks(range(n_dom))
    ax.set_xticklabels(dom_names, rotation=45, ha='right', fontsize=9, fontweight='bold')
    ax.set_yticks(range(n_dom))
    ax.set_yticklabels(dom_names, fontsize=9, fontweight='bold')

    for i, (name, _, _, color) in enumerate(DOMAIN_GROUPS):
        ax.get_xticklabels()[i].set_color(color)
        ax.get_yticklabels()[i].set_color(color)

    color_max = vmax if vmax is not None else np.nanmax(mat)
    for i in range(n_dom):
        for j in range(n_dom):
            v = mat[i, j]
            if v > 0:
                txt = f'{v:.1f}' if v < 100 else f'{v:.0f}'
                tc = 'white' if min(v, color_max) > 0.5 * color_max else 'black'
                ax.text(j, i, txt, ha='center', va='center',
                        fontsize=7, fontweight='bold', color=tc)

    cbar = plt.colorbar(im, ax=ax, shrink=0.8, extend=('max' if vmax is not None else 'neither'))
    cbar.set_label('Sum of contact frequencies', fontsize=10)

    ax.set_title(f'FL PARP14 \u2014 {contact_label} Contact Sum by Domain (CG-MD)\n'
                 f'CA-CA \u2264 {LYS_LYS_CUTOFF} nm, 25 replicates averaged{subtitle}',
                 fontsize=13, fontweight='bold')

    plt.tight_layout()
    for ext in ['png', 'svg']:
        fig.savefig(os.path.join(FIG_PATH, f'{fname}.{ext}'),
                    dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved: figures/{fname}")


def plot_domain_contact_sum(pairs, contact_label, fname):
    """Domain x domain heatmap: sum of per-residue-pair contact frequencies.

    Saves two versions: the plain heatmap (color scale spans the full matrix,
    so the large intra-domain diagonal dominates), and a `_v2` variant whose
    color scale is capped at the largest INTER-domain value -- the diagonal
    entries (residues within one domain, always in close contact) run to
    ~600 and otherwise wash out the off-diagonal (inter-domain) contrast,
    which tops out around ~95.
    """
    n_dom = len(DOMAIN_GROUPS)
    dom_names = [name for name, _, _, _ in DOMAIN_GROUPS]

    # Build domain lookup: resid -> domain index
    resid_to_dom = {}
    for di, (name, start, end, _) in enumerate(DOMAIN_GROUPS):
        for r in range(start, end + 1):
            resid_to_dom[r] = di

    # Sum contact frequencies per domain pair
    mat = np.zeros((n_dom, n_dom))
    for (ri, rj), freq in pairs.items():
        di = resid_to_dom.get(ri)
        dj = resid_to_dom.get(rj)
        if di is not None and dj is not None:
            mat[di, dj] += freq
            if di != dj:
                mat[dj, di] += freq

    _render_domain_sum_heatmap(mat, dom_names, contact_label, fname)

    offdiag = mat[~np.eye(n_dom, dtype=bool)]
    offdiag_max = offdiag[offdiag > 0].max() if np.any(offdiag > 0) else None
    if offdiag_max is not None:
        _render_domain_sum_heatmap(
            mat, dom_names, contact_label, f'{fname}_v2', vmax=offdiag_max,
            subtitle=f'\ncolor scale capped at max inter-domain value ({offdiag_max:.0f})')


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--set', default='fl',
                        help="Set name (e.g. fl, fl_optimized)")
    parser.add_argument('--sim-folder', nargs='+', default=None, metavar='PATH',
                        help='One or more simulation-set folders to analyze (each '
                             'containing seed-*_sample-*/ replicates). Domain units '
                             "are read from the folder's metadata.json, or pass --units.")
    parser.add_argument('--units', nargs='+', default=None, metavar='UNIT',
                        help='FL domain units in the --sim-folder construct '
                             '(e.g. md1l1 md2 md3). Required only if the folder has '
                             'no metadata.json.')
    parser.add_argument('--force', action='store_true',
                        help='No-op here (this script always recomputes); accepted '
                             'so run_analysis.sh can forward the same flags to every module.')
    args = parser.parse_args()

    if args.sim_folder:
        # This script analyzes one set at a time; use the first folder.
        set_key = register_sim_folder(args.sim_folder[0], units=args.units)
    else:
        set_key = args.set

    print("=" * 60)
    print(f"{set_key.upper()} Lysine Contact Analysis (CALVADOS CG-MD)")
    print("=" * 60)

    jobs = [(s, p) for s in SEEDS for p in SAMPLES]
    print(f"  Processing {len(jobs)} {set_key} replicates...")

    results = []
    n_workers = min(25, os.cpu_count() or 4)
    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futures = [pool.submit(worker_lys_contacts, s, p, set_key)
                   for s, p in jobs]
        for f in futures:
            r = f.result()
            if r is not None:
                results.append(r)
                print(f"    seed-{r['seed']}_sample-{r['sample']}: "
                      f"{r['n_frames']} frames, "
                      f"{len(r['ll_freq'])} Lys-Lys, "
                      f"{len(r['la_freq'])} Lys-acid pairs")

    print(f"\n  {len(results)} replicates completed")

    # Average pair frequencies across replicates
    ll_all = defaultdict(list)
    la_all = defaultdict(list)
    for r in results:
        for k, v in r['ll_freq'].items():
            ll_all[k].append(v)
        for k, v in r['la_freq'].items():
            la_all[k].append(v)

    ll_mean = {k: np.mean(v) for k, v in ll_all.items()}
    la_mean = {k: np.mean(v) for k, v in la_all.items()}

    # Save
    out_npz = os.path.join(DATA_PATH, f'{set_key}_lys_contacts.npz')
    np.savez(out_npz,
             ll_pairs=np.array(list(ll_mean.keys())),
             ll_freq=np.array(list(ll_mean.values())),
             la_pairs=np.array(list(la_mean.keys())),
             la_freq=np.array(list(la_mean.values())))
    print(f"  Saved: {out_npz}")

    # Plot
    label = 'Full-Length' if set_key == 'fl' else set_key.replace('_', ' ').title()
    print("\nPlotting...")
    plot_contact_map(ll_mean,
                     f'{label} PARP14 \u2014 Lys-Lys Contacts (CG-MD)',
                     f'{set_key}_lys_lys_contacts',
                     f'CA-CA \u2264 {LYS_LYS_CUTOFF} nm')

    plot_contact_map(la_mean,
                     f'{label} PARP14 \u2014 Lys-Acidic Contacts (CG-MD)',
                     f'{set_key}_lys_acidic_contacts',
                     f'CA-CA \u2264 {LYS_ACIDIC_CUTOFF} nm')

    # Domain-level summed contact frequency heatmaps
    print("\nPlotting domain-level sums...")
    plot_domain_contact_sum(ll_mean, 'Lys-Lys',
                            f'{set_key}_lys_lys_domain_sum')
    plot_domain_contact_sum(la_mean, 'Lys-Acidic',
                            f'{set_key}_lys_acidic_domain_sum')

    print("\nDone!")


if __name__ == '__main__':
    main()
