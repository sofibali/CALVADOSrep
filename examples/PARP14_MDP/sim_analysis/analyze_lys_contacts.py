#!/usr/bin/env python3
"""
Lysine crosslink-feasibility analysis from CALVADOS CG-MD trajectories.

WHAT THE CUTOFF MEANS
---------------------
This is a **DSS crosslinking** model, not a generic contact analysis. The 3.0 nm
CA-CA cutoff is deliberate and chemistry-derived, not a loose contact radius:

    DSS spacer arm      ~15 A
    + 2 x Lys side chain reach
    ------------------------------
    = ~30 A = 3.0 nm  CA-CA ceiling

which is the standard ceiling used in the XL-MS field for DSS. A CG model with
one bead per residue has no side chains to resolve, so the permissive CA-CA
ceiling is the right level of description here.

(An earlier version of this docstring claimed 1.6 nm / 0.8 nm. That was wrong --
the code was always 3.0 nm and the code was correct.)

TWO CORRECTIONS THE RAW CUTOFF NEEDS
------------------------------------
1. **Sequence separation.** At 3.0 nm, residues close in sequence are within
   reach purely from chain connectivity, so they satisfy the cutoff in EVERY
   frame regardless of structure. The threshold is DERIVED, not chosen:
   ceil(cutoff / 0.38 nm bond) + 1 = 9. Below that a fully extended CG chain
   still cannot exceed the ceiling. Confirmed by the data -- separations 1-8 are
   100% saturated at persistence 1.000; 9 is the first with any pair below 1.0.
   Raising it discards real signal (the 10-19 band holds 245 fl pairs spanning
   persistence 0.36-0.98).

2. **Euclidean distance overestimates reachability.** A straight line between
   two lysines can pass straight through the protein, which a crosslinker
   cannot. The physically meaningful quantity is the **solvent-accessible
   surface distance (SASD)** -- the shortest path that stays outside the
   protein, which is what Xwalk computes. `--sasd` enables it (see
   `compute_sasd`). Pairs whose Euclidean distance passes but whose SASD does
   not are false positives, and they are common for buried pairs.

Computes, per frame, Lys-Lys and Lys-acidic (Asp/Glu) CA-CA pair frequencies,
and produces scatter contact maps on the 1-1801 residue grid.

Usage:
    python analyze_lys_contacts.py --set fl --sasd      # derived seq-sep, + SASD
    python analyze_lys_contacts.py --set fl --target-frames 0    # every frame
    bash run_crosslinks_all.sh                          # every set with data
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

# DSS crosslink ceiling: ~15 A spacer arm + 2 x Lys side-chain reach ~= 30 A
# CA-CA. Permissive by design -- see the module docstring. Keep both the same:
# the Lys-acidic map is the same reachability question for a different partner.
LYS_LYS_CUTOFF = 3.0     # nm (CA-CA, DSS-derived)
LYS_ACIDIC_CUTOFF = 3.0   # nm (CA-CA, DSS-derived)

# Minimum |resid_i - resid_j| for a pair to be counted.
#
# DERIVED, not chosen. A CG chain has a uniform 0.38 nm CA-CA bond, so two
# residues separated by n along the chain cannot be further apart than
# n * 0.38 nm even fully extended. While that maximum stays under the DSS
# ceiling, the pair is in contact in EVERY frame no matter what the structure
# does -- the observation is guaranteed by connectivity and carries zero
# tertiary-structure information.
#
#     n * 0.38 <= 3.0   ->   n <= 7.9
#
# so ceil(3.0/0.38) = 8 is the first separation that can exceed the ceiling AT
# ALL. But it only manages it at full extension (8 * 0.38 = 3.04 nm, 1% of
# slack), which a flexible chain effectively never reaches -- so 8 is still
# saturated in practice and the first genuinely informative separation is one
# further out. Hence the +1.
#
# Confirmed against the unfiltered fl data, which shows a sharp cliff exactly
# there: |d| = 1..8 are 100% saturated at persistence 1.000, |d| = 9 is the
# first with any pair below 1.0 (5%), and from 10 on the mean falls smoothly
# (0.97 -> 0.51 by 18) with no pair permanently in contact.
#
# Set higher only to deliberately restrict to long-range restraints, and know
# that it discards real signal: the 10-19 band holds 245 fl pairs with mean
# persistence 0.71 and a 0.36-0.98 spread.
import math as _math
CG_BOND_LENGTH = 0.38    # nm, uniform across all 20 residues in CALVADOS3
MIN_SEQ_SEP = _math.ceil(LYS_LYS_CUTOFF / CG_BOND_LENGTH) + 1   # = 9 at 3.0 nm

# --- Solvent-accessible surface distance (SASD), Xwalk-style ---
SASD_BEAD_RADIUS = 0.30   # nm, obstacle radius per CA bead (CALVADOS sigma/2)
SASD_GRID = 0.15          # nm, voxel edge for the geodesic search
SASD_MARGIN = 1.20        # nm, padding around the endpoint bounding box
SASD_MAX = 3.5            # nm, give up beyond this (DSS ceiling + slack)

SKIP_FRAMES = 50
# Per-replicate frame budget; the stride is chosen per trajectory. None = every
# frame. See the note in worker_lys_contacts for why a fixed stride is wrong here.
TARGET_FRAMES = 2000
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


def compute_sasd(pos, i, j, bead_radius=SASD_BEAD_RADIUS, grid=SASD_GRID,
                 margin=SASD_MARGIN, max_dist=SASD_MAX):
    """Solvent-accessible surface distance between beads i and j, in nm.

    The Euclidean CA-CA distance lets a crosslinker pass straight through the
    protein. Xwalk instead measures the shortest path that stays in solvent,
    and that is what decides whether a crosslink is actually formable. This is
    the same idea adapted to a one-bead-per-residue model.

    Method: voxelize a local box around the two endpoints, mark a voxel blocked
    if it lies within `bead_radius` of any bead other than the two endpoints,
    then run Dijkstra with 26-connectivity over the free voxels. Returns the
    geodesic path length, or np.inf if no path shorter than `max_dist` exists.

    Approximation to state plainly: a CG bead is a whole residue, so the
    "surface" here is coarser than Xwalk's all-atom one. Treat SASD from this
    function as a reachability filter, not as a quantitative distance.
    """
    import heapq

    a, b = pos[i], pos[j]
    if np.linalg.norm(b - a) > max_dist:
        return np.inf

    lo = np.minimum(a, b) - margin
    hi = np.maximum(a, b) + margin
    dims = np.ceil((hi - lo) / grid).astype(int) + 1
    if dims.prod() > 4_000_000:          # runaway box; caller can widen the grid
        return np.nan

    # Only beads that can reach into the box matter.
    near = np.all((pos > lo - bead_radius) & (pos < hi + bead_radius), axis=1)
    near[i] = near[j] = False
    obstacles = pos[near]

    free = np.ones(dims, dtype=bool)
    if len(obstacles):
        gx = np.arange(dims[0]) * grid + lo[0]
        gy = np.arange(dims[1]) * grid + lo[1]
        gz = np.arange(dims[2]) * grid + lo[2]
        r2 = bead_radius ** 2
        for ob in obstacles:
            # bound the sphere to its own sub-box instead of scanning the grid
            i0 = np.maximum(((ob - bead_radius - lo) / grid).astype(int), 0)
            i1 = np.minimum(((ob + bead_radius - lo) / grid).astype(int) + 1, dims)
            if np.any(i0 >= i1):
                continue
            sx = gx[i0[0]:i1[0]] - ob[0]
            sy = gy[i0[1]:i1[1]] - ob[1]
            sz = gz[i0[2]:i1[2]] - ob[2]
            d2 = (sx[:, None, None] ** 2 + sy[None, :, None] ** 2
                  + sz[None, None, :] ** 2)
            free[i0[0]:i1[0], i0[1]:i1[1], i0[2]:i1[2]] &= (d2 > r2)

    start = tuple(np.clip(((a - lo) / grid).round().astype(int), 0, dims - 1))
    goal = tuple(np.clip(((b - lo) / grid).round().astype(int), 0, dims - 1))
    # Endpoints sit inside their own bead; carve them free so the search can start.
    free[start] = free[goal] = True
    if not free[start] or not free[goal]:
        return np.inf

    # 26-connected neighbour offsets with their true step lengths
    offs = [(dx, dy, dz) for dx in (-1, 0, 1) for dy in (-1, 0, 1)
            for dz in (-1, 0, 1) if (dx, dy, dz) != (0, 0, 0)]
    steps = [(o, np.sqrt(o[0] ** 2 + o[1] ** 2 + o[2] ** 2) * grid) for o in offs]

    dist = np.full(dims, np.inf)
    dist[start] = 0.0
    heap = [(0.0, start)]
    while heap:
        d, cur = heapq.heappop(heap)
        if cur == goal:
            return float(d)
        if d > dist[cur] or d > max_dist:
            continue
        cx, cy, cz = cur
        for (dx, dy, dz), w in steps:
            nx, ny, nz = cx + dx, cy + dy, cz + dz
            if not (0 <= nx < dims[0] and 0 <= ny < dims[1] and 0 <= nz < dims[2]):
                continue
            if not free[nx, ny, nz]:
                continue
            nd = d + w
            if nd < dist[nx, ny, nz] and nd <= max_dist:
                dist[nx, ny, nz] = nd
                heapq.heappush(heap, (nd, (nx, ny, nz)))
    return np.inf


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


def _resolve_paths(set_key, seed, sample):
    """(pdb, dcd) for one replicate. Shared by the contact worker and the SASD
    pass so both agree on which trajectory a replicate means."""
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
        # NOT 'parp14.dcd' hardcoded -- every named set except fl/fl_optimized/
        # fl_go has a different sysname (e.g. 'parp14_macrodomains' for 'md'),
        # so a literal 'parp14.dcd' silently finds nothing for those.
        sysname = reg.SETS[set_key]['sysname']
        sim_dir = os.path.join(CWD, set_key, f'seed-{seed}_sample-{sample}')
        dcd = os.path.join(sim_dir, f'{sysname}.dcd')

    pdb = os.path.join(sim_dir, 'checkpoint.pdb')
    if not os.path.isfile(pdb):
        pdb = os.path.join(sim_dir, 'top.pdb')
    if not os.path.isfile(pdb):
        pdb = os.path.join(sim_dir, 'restart.pdb')
    return pdb, dcd


def worker_lys_contacts(seed, sample, set_key='fl'):
    """Compute Lys contacts for one replicate of a given set."""
    import MDAnalysis as mda

    pdb, dcd = _resolve_paths(set_key, seed, sample)
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

    # Per-trajectory stride so every replicate contributes a comparable number
    # of frames. The sets span 3-4k-frame originals and 100k-frame extensions
    # (20 ps/frame), and contact persistence is a slowly-varying observable --
    # analysing all 100k costs ~50x more for no extra independent sampling.
    step = 1
    if TARGET_FRAMES:
        usable = max(0, len(u.trajectory) - SKIP_FRAMES)
        step = max(1, usable // TARGET_FRAMES)

    for ts in u.trajectory[SKIP_FRAMES::step]:
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
                if abs(ri - rj) < MIN_SEQ_SEP:
                    continue
                ll_pair_counts[(min(ri, rj), max(ri, rj))] += 1

        # Vectorized Lys-acidic distances
        acid_pos = pos[acidic_idx]
        for i in range(n_lys):
            dists = np.linalg.norm(acid_pos - lys_pos[i], axis=1)
            contacts = np.where(dists <= LYS_ACIDIC_CUTOFF)[0]
            for ci in contacts:
                ri = int(lys_resids[i])
                rj = int(acid_resids[ci])
                if abs(ri - rj) < MIN_SEQ_SEP:
                    continue
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


def run_sasd_pass(set_key, ll_mean, la_mean, local_to_fl, top_n=200, n_frames=5):
    """Re-score the most persistent pairs by solvent-accessible surface distance.

    The Euclidean pass answers "are these two lysines within DSS reach as the
    crow flies". This answers the question that actually decides crosslink
    formability: "is there a route between them that stays in solvent". Pairs
    that pass the first test and fail this one are buried false positives.

    Runs on a handful of frames of one replicate rather than the full ensemble --
    a geodesic search per pair per frame is orders of magnitude more expensive
    than a distance, and burial is a slowly-varying property.
    """
    import csv as _csv
    import MDAnalysis as mda
    import warnings
    warnings.filterwarnings('ignore')

    pdb, dcd = _resolve_paths(set_key, 1, 0)
    if not (os.path.isfile(pdb) and os.path.isfile(dcd)):
        print(f"  SASD: no trajectory for {set_key} seed-1_sample-0, skipping")
        return

    # FL -> local, to index back into the trajectory
    fl_to_local = reg.build_fl_to_construct_map(reg.get_units(set_key), set_key=set_key)

    # Rank INTER-DOMAIN pairs first, then intra-domain.
    #
    # A plain global sort by persistence spends the whole budget on
    # intra-domain pairs and returns zero inter-domain ones: measured on fl,
    # inter-domain pairs average 0.021 persistence against 0.135 intra-domain,
    # so they are systematically rarer and never reach the top of a global
    # ranking. But intra-domain crosslinks report on a fold that does not change
    # between constructs -- only inter-domain ones report on architecture, which
    # is the whole point of comparing constructs. So they get the budget first.
    def _domain_of(fl_resid):
        for name, (lo, hi) in reg.FL_DOMAINS.items():
            if lo <= fl_resid <= hi:
                return name
        return None

    def _is_inter(pair):
        di, dj = _domain_of(pair[0]), _domain_of(pair[1])
        return di is not None and dj is not None and di != dj

    cand = list(ll_mean.items()) + list(la_mean.items())
    inter = sorted([kv for kv in cand if _is_inter(kv[0])], key=lambda kv: -kv[1])
    intra = sorted([kv for kv in cand if not _is_inter(kv[0])], key=lambda kv: -kv[1])
    ranked = (inter + intra)[:top_n]
    n_inter = sum(1 for kv in ranked if _is_inter(kv[0]))
    if not ranked:
        print("  SASD: no pairs survived the Euclidean pass, skipping")
        return
    print(f"  SASD candidates: {n_inter} inter-domain + {len(ranked) - n_inter} "
          f"intra-domain (of {len(inter)} / {len(intra)} available)")

    u = mda.Universe(pdb, dcd)
    resid_to_idx = {int(r): i for i, r in enumerate(u.atoms.resids)}
    total = len(u.trajectory)
    picks = np.linspace(SKIP_FRAMES, total - 1, n_frames).astype(int)

    print(f"\n  SASD pass: {len(ranked)} pairs x {len(picks)} frames "
          f"(Xwalk-style geodesic, grid {SASD_GRID} nm)")

    acc = {pair: [] for pair, _ in ranked}
    for t in picks:
        u.trajectory[int(t)]
        pos = u.atoms.positions / 10.0
        for (fi, fj), _ in ranked:
            li, lj = fl_to_local(fi), fl_to_local(fj)
            if li is None or lj is None:
                continue
            ii, jj = resid_to_idx.get(li), resid_to_idx.get(lj)
            if ii is None or jj is None:
                continue
            acc[(fi, fj)].append(compute_sasd(pos, ii, jj))

    out = os.path.join(DATA_PATH, f'{set_key}_lys_sasd.csv')
    n_reach = 0
    with open(out, 'w', newline='') as fh:
        w = _csv.writer(fh)
        w.writerow(['resid_i', 'resid_j', 'seq_sep', 'domain_i', 'domain_j',
                    'inter_domain', 'euclid_persistence',
                    'sasd_mean_nm', 'sasd_min_nm', 'n_frames_reachable', 'reachable'])
        for (fi, fj), freq in ranked:
            vals = [v for v in acc[(fi, fj)] if np.isfinite(v)]
            reach = len(vals)
            # "reachable" = a solvent path exists within the DSS ceiling in at
            # least one sampled frame; a crosslink only needs one conformer.
            ok = reach > 0 and min(vals) <= LYS_LYS_CUTOFF
            n_reach += bool(ok)
            di, dj = _domain_of(fi), _domain_of(fj)
            w.writerow([fi, fj, abs(fi - fj), di or '-', dj or '-',
                        int(bool(di and dj and di != dj)), f'{freq:.4f}',
                        f'{np.mean(vals):.3f}' if vals else 'inf',
                        f'{min(vals):.3f}' if vals else 'inf',
                        reach, int(ok)])
    print(f"  SASD: {n_reach}/{len(ranked)} pairs have a solvent path within "
          f"{LYS_LYS_CUTOFF} nm")
    print(f"  Saved: {os.path.relpath(out, CWD)}")


def main():
    import argparse
    # Declared up front: MIN_SEQ_SEP is read below as an argparse default and
    # rebound after parsing, and Python requires the declaration to precede the
    # first use in the function.
    global MIN_SEQ_SEP, TARGET_FRAMES, FIG_PATH
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
    parser.add_argument('--min-seq-sep', type=int, default=MIN_SEQ_SEP,
                        metavar='N',
                        help=f'Drop pairs with |resid_i - resid_j| < N (default '
                             f'{MIN_SEQ_SEP}, derived as cutoff/bond_length + 1). '
                             f'Below that a fully extended CG chain still cannot '
                             f'exceed the DSS ceiling, so the contact is forced by '
                             f'connectivity and says nothing about structure. '
                             f'Raising it restricts to long-range restraints but '
                             f'discards real signal.')
    parser.add_argument('--target-frames', type=int, default=TARGET_FRAMES,
                        metavar='N',
                        help=f'Per-replicate frame budget; the stride is chosen '
                             f'PER TRAJECTORY (default {TARGET_FRAMES}). The sets '
                             f'span 3-4k-frame originals and 100k-frame extensions, '
                             f'so one fixed stride cannot serve both. Pass 0 to use '
                             f'every frame.')
    parser.add_argument('--sasd', action='store_true',
                        help='After the Euclidean pass, recompute the top pairs as '
                             'solvent-accessible surface distance (Xwalk-style): the '
                             'shortest path that stays OUTSIDE the protein. Filters '
                             'out pairs whose straight line passes through the core.')
    parser.add_argument('--sasd-top', type=int, default=200, metavar='N',
                        help='How many of the most persistent pairs to re-score with '
                             'SASD (default 200). SASD is a per-pair geodesic search, '
                             'far more expensive than a distance.')
    parser.add_argument('--sasd-frames', type=int, default=5, metavar='N',
                        help='Frames per replicate to average SASD over (default 5).')
    args = parser.parse_args()

    # Set before any worker process is forked so they inherit the values.
    MIN_SEQ_SEP = args.min_seq_sep
    TARGET_FRAMES = args.target_frames or None

    if args.sim_folder:
        # This script analyzes one set at a time; use the first folder.
        set_key = register_sim_folder(args.sim_folder[0], units=args.units)
    else:
        set_key = args.set

    FIG_PATH = str(_get_fig_dir('07_lysine_contacts', sims=[set_key]))

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

    # Remap construct-LOCAL residue numbers -> FL numbering before anything
    # gets plotted. plot_contact_map()/plot_domain_contact_sum() below both
    # assume FL numbering (DOMAIN_GROUPS boundaries + the 1-1801 axis are FL-
    # numbered) -- without this, every non-FL set (md, core, mka, norrm,
    # noart, md3art, md_full, mka_full) gets its local resi 1 plotted at FL
    # position 1 (RRM1) instead of its true FL position, squeezing the whole
    # construct into one corner of the grid instead of its real location.
    local_to_fl = reg.build_construct_to_fl_map(reg.get_units(set_key), set_key=set_key)

    def _remap(pairs):
        out = {}
        for (ri, rj), v in pairs.items():
            fi, fj = local_to_fl(ri), local_to_fl(rj)
            if fi is not None and fj is not None:
                out[(fi, fj)] = v
        return out

    ll_mean = _remap(ll_mean)
    la_mean = _remap(la_mean)

    if args.sasd:
        run_sasd_pass(set_key, ll_mean, la_mean, local_to_fl,
                      top_n=args.sasd_top, n_frames=args.sasd_frames)

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
