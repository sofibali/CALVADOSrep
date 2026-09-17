#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Inter-domain distance violin plot for fl_optimized PARP14 simulations.

For each frame in each replicate, compute:
  - MD1L1 - MD3 center-of-mass distance
  - MD1L1 - ART center-of-mass distance
  - (Optional) minimum CA-CA distance between domain pairs

Then plot distributions as violins. The "frequency" axis comes from the
density of frames at each distance.

This is a focused alternative to the grid contact-frequency heatmap when you
want to see the full distribution shape (bimodal, narrow, etc.) for specific
inter-domain pairs.

Usage:
    python figure_md_distances.py
    python figure_md_distances.py --set fl fl_optimized   # compare multiple sets
    python figure_md_distances.py --metric min            # use min CA-CA instead of COM
"""

import os
import argparse
import numpy as np
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import os as _os_boot, sys as _sys_boot  # __ROOTBOOT__ (script lives in a subfolder; root=parent)
_sys_boot.path.insert(0, _os_boot.path.dirname(_os_boot.path.dirname(_os_boot.path.abspath(__file__))))
import sim_registry as reg
import parp14_mdp_colors as _mdp_colors

def _set_color_lookup(sets):
    """color per set: sim_registry.SETS first, then --sim-folder-as aliases
    (e.g. 'fl_go_5rep1us') via parp14_mdp_colors.CONSTRUCT, which is not
    registered in sim_registry (register_sim_folder only sets EXTERNAL_DIRS)."""
    out = {}
    for sk in sets:
        if sk in reg.SETS:
            out[sk] = reg.SETS[sk]['color']
        elif sk in _mdp_colors.CONSTRUCT:
            out[sk] = _mdp_colors.CONSTRUCT[sk]
    return out

# ============================================================
# Configuration
# ============================================================

CWD = Path(__file__).resolve().parent.parent
DATA_PATH = CWD / 'data'
DATA_PATH.mkdir(exist_ok=True)

# Arbitrary simulation folders registered via --sim-folder.
# Populated in the main process before any worker pool is created so that
# forked workers (Linux fork start method) inherit them.
EXTERNAL_DIRS = {}      # set_key -> absolute Path to the sim folder
EXTERNAL_SYSNAME = {}   # set_key -> dcd basename (sysname)
CONSTRUCT_UNITS = {}    # set_key -> FL domain units (drives FL->construct remap)

# Dated category subdir: figures/04_md_distances/<YYYY-MM-DD>/
import sys as _sys
_sys.path.insert(0, str(CWD))
from _fig_layout import get_fig_dir as _get_fig_dir
FIG_PATH = _get_fig_dir('04_md_distances')

# Domain definitions in FL numbering
DOMAINS_FL = {
    'MD1L1': (790, 1004),
    'MD2':   (1005, 1193),
    'MD3':   (1207, 1388),
    'ART':   (1603, 1801),
}

# Pairs to plot
PAIRS = [
    ('MD1L1', 'MD3'),
    ('MD1L1', 'ART'),
    # Additional pairs to add to violin row for context:
    ('MD1L1', 'MD2'),
    ('MD2',   'MD3'),
    ('MD3',   'ART'),
    ('MD2',   'ART'),  # completes all 6 pairwise combinations of the 4 domains
]

# Contact-distance reference thresholds. These are BEAD-BEAD (residue-scale)
# length scales from the actual CALVADOS potentials this project's sims use
# (calvados/interactions.py; params confirmed from fl_go's own config.yaml:
# eps_lj=0.2, cutoff_lj=2.0, ionic=0.19, temp=293 -> Debye length 0.70 nm).
# For metric='min' (closest residue-residue approach), AH_CUTOFF_NM is the
# actual Ashbaugh-Hatch (van der Waals-like) cutoff -- beyond it the AH force
# is exactly zero, so "min CA-CA < 2.0 nm" directly means "these two residues
# COULD have nonzero AH energy," no domain-size correction needed.
AH_CUTOFF_NM = 2.0
# Legacy alias (kept so any external caller/import doesn't break).
CONTACT_CUTOFF_NM = AH_CUTOFF_NM
SKIP_FRAMES = 50  # equilibration

SEEDS = range(1, 6)
SAMPLES = range(0, 5)


# ============================================================
# External simulation-folder registration (--sim-folder)
# ============================================================

def register_sim_folder(path, units=None, name=None):
    """Register an arbitrary simulation folder (from --sim-folder) as a set so
    this figure can be generated for it without editing the script. Domain units
    are read from the folder's metadata.json, the `units` arg, or (if the
    basename is a known set) that set. Returns the set_key it is registered under.

    `name` overrides the set_key (default: the folder's basename) -- lets the
    SAME physical folder be registered twice under two different set_keys
    (e.g. 'fl_go' with full data and 'fl_go_5rep1us' with REPLICATE_CAP/
    FRAME_CAP applied below) so both can be compared side by side.

    Call this in the main process before launching the worker pool so that
    forked workers inherit EXTERNAL_DIRS / EXTERNAL_SYSNAME / CONSTRUCT_UNITS.
    """
    folder = Path(path).resolve()
    if not folder.is_dir():
        raise FileNotFoundError(f"sim folder not found: {folder}")
    resolved_units = reg._resolve_units(str(folder), units)
    meta = reg.read_metadata(str(folder)) or {}
    sysname = meta.get('sysname') or reg.detect_sysname(str(folder))
    set_key = name or folder.name
    EXTERNAL_DIRS[set_key] = folder
    EXTERNAL_SYSNAME[set_key] = sysname
    CONSTRUCT_UNITS[set_key] = resolved_units
    return set_key


# set_key -> max number of replicates to include (first N in (seed, sample)
# enumeration order -- for --sim-folder flat-dir sets this is the first N
# dirs in sorted order). Populated from --max-reps before job-building.
# set_key -> sample index to keep (drop every other sample for that set).
# For sets where only ONE sample index per seed was extended to the full
# target length (e.g. md_full: seed-{1..5}_sample-0 extended to ~1000 ns,
# sample-{1,2,3,4} left at the original ~100 ns) -- --max-reps's "first N in
# enumeration order" would wrongly grab seed-1's samples 0-4 instead of the
# actually-extended replicates, so this filters by sample index directly.
ONLY_SAMPLE = {}

REPLICATE_CAP = {}

# set_key -> max number of POST-SKIP_FRAMES frames to read per replicate
# (wfreq=1000 steps/frame, dt=0.01 ps -> 100 frames/ns, see calvados/sim.py's
# LangevinIntegrator). Populated from --max-ns before job-building.
FRAME_CAP = {}
FRAMES_PER_NS = 100


# Cache: set_key -> replicate dirs discovered by _flat_replicate_dirs() below.
_FLAT_DIRS_CACHE = {}


def _flat_replicate_dirs(folder, sysname):
    """Every immediate subdirectory of `folder` containing a `<sysname>.dcd`.

    Fallback for --sim-folder layouts that aren't a seed-{1-5}_sample-{0-4}
    grid -- e.g. a handful of independent long runs named arbitrarily
    (fl_go's state-*_tica_seed-*_sample-*_fr*/ dirs)."""
    return [d for d in sorted(folder.iterdir())
            if d.is_dir() and (d / f'{sysname}.dcd').is_file()]


def _resolve_sim_paths(set_key, seed, sample):
    if set_key in EXTERNAL_DIRS:
        base = EXTERNAL_DIRS[set_key]
        sysname = EXTERNAL_SYSNAME[set_key]
        sim_dir = base / f'seed-{seed}_sample-{sample}'
        if not sim_dir.is_dir():
            if set_key not in _FLAT_DIRS_CACHE:
                _FLAT_DIRS_CACHE[set_key] = _flat_replicate_dirs(base, sysname)
            flat = _FLAT_DIRS_CACHE[set_key]
            idx = (seed - 1) * len(SAMPLES) + sample
            if idx < len(flat):
                sim_dir = flat[idx]
    elif set_key.startswith('frag_'):
        sysname = 'parp14'  # every fragment uses this sysname by convention
        sim_dir = CWD / 'fragments' / set_key[5:] / f'seed-{seed}_sample-{sample}'
    else:
        # NOT hardcoded 'parp14' -- most named sets have a different sysname
        # (e.g. 'parp14_macrodomains' for 'md'), so a literal 'parp14' silently
        # finds no .dcd for those.
        sysname = reg.SETS[set_key]['sysname']
        sim_dir = CWD / set_key / f'seed-{seed}_sample-{sample}'
    dcd = sim_dir / f'{sysname}.dcd'
    for name in ('top.pdb', 'restart.pdb', 'checkpoint.pdb'):
        pdb = sim_dir / name
        if pdb.is_file():
            return pdb, dcd
    return sim_dir / 'top.pdb', dcd


def domain_ranges_for_set(set_key):
    """Domain residue ranges (construct numbering) for the inter-domain pairs.

    Always remaps DOMAINS_FL into the construct's own numbering via
    sim_registry -- CONSTRUCT_UNITS (this file's local dict) was ONLY ever
    populated for --sim-folder-registered sets (see register_sim_folder
    below), never for plain --set NAME sets, so every named non-FL-identity
    set (md, mka, md3art, md_full, ...) used to silently fall back to raw
    DOMAINS_FL unchanged. For fl/fl_optimized/fl_go (genuinely FL-numbered)
    that's harmless; for real sub-constructs it either crashed outright
    (resid ranges exceeding the construct's actual residue count, e.g. 'md'
    only has 586 residues but DOMAINS_FL reaches up to 1801) or, worse,
    silently selected the WRONG residues where the FL range happened to
    still fall within the construct's local numbering (core/norrm/noart/mka
    are all large enough for this) -- producing plausible-looking but
    mislabeled inter-domain distances with no error at all.
    """
    units = CONSTRUCT_UNITS.get(set_key) or reg.get_units(set_key)
    fl_to_c = reg.build_fl_to_construct_map(units, set_key=set_key)
    mapped = {}
    for dname, (fl_s, fl_e) in DOMAINS_FL.items():
        c_s = fl_to_c(fl_s)
        c_e = fl_to_c(fl_e)
        if c_s is not None and c_e is not None:
            mapped[dname] = (c_s, c_e)
    return mapped


# ============================================================
# Worker: compute distances for one replicate
# ============================================================

def compute_distances_one_replicate(args):
    """Compute per-frame inter-domain distances AND per-domain Rg for one
    replicate.

    Rg is computed alongside the distances (same per-frame loop, negligible
    extra cost) because a domain's own size sets a hard floor on how close
    its COM can approach another domain's COM -- two folded, sterically-
    excluding globular domains can't have their centroids closer than
    roughly Rg_A + Rg_B (surfaces touching). That floor is real physics
    (excluded volume + the domains' own internal restraints keeping them
    compact), not a bead-model artifact, and it's per-domain-size, not a
    universal constant -- so a single fixed CONTACT_CUTOFF_NM tuned for
    min-CA-CA (residue-scale) contact is meaningless for COM-COM (domain-
    scale) distances; see Rg_A+Rg_B computed downstream from this.

    Returns dict {'distances': {pair_label: array}, 'domain_rg': {domain:
    array}} or None on failure.
    """
    set_key, seed, sample, metric = args
    import MDAnalysis as mda

    pdb, dcd = _resolve_sim_paths(set_key, seed, sample)
    if not pdb.is_file() or not dcd.is_file():
        return None

    u = mda.Universe(str(pdb), str(dcd))

    # Pre-select domain atom groups
    domains = domain_ranges_for_set(set_key)
    ags = {}
    for dname, (ds, de) in domains.items():
        ag = u.select_atoms(f'resid {ds}:{de}')
        if len(ag) > 0:
            ags[dname] = ag

    if not ags:
        return None

    # Initialize result arrays
    results = {f'{a}_{b}': [] for a, b in PAIRS if a in ags and b in ags}
    rg_results = {d: [] for d in ags}

    if set_key in FRAME_CAP:
        frame_slice = u.trajectory[SKIP_FRAMES:SKIP_FRAMES + FRAME_CAP[set_key]]
    else:
        frame_slice = u.trajectory[SKIP_FRAMES:]

    for ts in frame_slice:
        # Per-domain COM + Rg once per frame (COM feeds both metric='com'
        # distances below and the Rg-based contact floor).
        com = {}
        for dname, ag in ags.items():
            pos = ag.positions / 10.0  # Å -> nm
            c = pos.mean(axis=0)
            com[dname] = c
            rg_results[dname].append(
                float(np.sqrt(np.mean(np.sum((pos - c) ** 2, axis=1)))))

        for a, b in PAIRS:
            key = f'{a}_{b}'
            if key not in results:
                continue
            if metric == 'com':
                d = float(np.linalg.norm(com[a] - com[b]))
            else:  # metric == 'min'
                # Minimum CA-CA distance (nm)
                pos_a = ags[a].positions / 10.0
                pos_b = ags[b].positions / 10.0
                diff = pos_a[:, np.newaxis, :] - pos_b[np.newaxis, :, :]
                d = float(np.sqrt(np.sum(diff ** 2, axis=2)).min())

            results[key].append(d)

    return {
        'set_key': set_key,
        'seed': seed,
        'sample': sample,
        'distances': {k: np.array(v) for k, v in results.items()},
        'domain_rg': {k: np.array(v) for k, v in rg_results.items()},
    }


# ============================================================
# Block-averaging (Flyvbjerg & Petersen), pooled across replicates
# ============================================================
#
# Companion to the violin/histogram plots above, and a more useful
# alternative to a fixed-cutoff contact-FREQUENCY summary when most frames
# sit well above any reasonable contact distance (diffuse, non-contacting
# states): a binary "in contact" count reads uniformly low and uninformative
# in that regime, but a continuous mean distance +/- a correlation-corrected
# SE stays meaningful regardless of whether the pair is ever actually in
# contact. See analyze_all.py's rg_convergence_block for the full derivation;
# the same method is applied here PER DOMAIN PAIR rather than to Rg, because
# each distance can have its own relaxation time -- there's no reason to
# assume it matches Rg's (a many-residue collective coordinate) or a TICA
# lag time (chosen for eigenvalue/timescale estimation, not for getting a
# valid SE on a raw distance).

def block_average_plateau(rep_arrays, dt_ns=0.01, min_pooled_blocks=30,
                           plateau_rel_tol=0.10):
    """Pooled-across-replicates block averaging for one observable.

    rep_arrays: list of 1D per-replicate arrays of the SAME observable
    (e.g. one inter-domain distance), assumed independent (different
    random seeds -- not relying on any intra-trajectory decorrelation
    assumption). At each block size b, every replicate independently
    contributes floor(min_frames/b) non-overlapping block means; all
    replicates' block means are pooled before computing SE/mean, so the
    sweep can extend to block size = full replicate length (where pooling
    gives exactly n_rep independent blocks).

    The plateau is the largest block size with >= min_pooled_blocks pooled
    blocks (trustworthy), extended toward smaller block sizes while SE
    stays within plateau_rel_tol of that reference -- i.e. where the curve
    has stopped rising. Returns None if there are too few replicates/frames
    to say anything.
    """
    n_rep = len(rep_arrays)
    if n_rep < 2:
        return None
    min_frames = min(len(a) for a in rep_arrays)
    if min_frames < 2:
        return None
    arr = np.array([np.asarray(a[:min_frames]) for a in rep_arrays])

    max_pow = int(np.floor(np.log2(max(min_frames, 1))))
    block_sizes = list(2 ** np.arange(0, max_pow + 1))
    if block_sizes[-1] != min_frames:
        block_sizes.append(min_frames)  # guarantee the full-length (1 block/replicate) point
    block_ns = np.array(block_sizes) * dt_ns

    se = np.full(len(block_sizes), np.nan)
    mean_est = np.full(len(block_sizes), np.nan)
    pooled_nb = np.zeros(len(block_sizes), dtype=int)
    for bi, b in enumerate(block_sizes):
        nb_per_rep = min_frames // b
        if nb_per_rep < 1:
            continue
        blocks = arr[:, :nb_per_rep * b].reshape(n_rep, nb_per_rep, b)
        pooled_means = blocks.mean(axis=2).ravel()
        pooled_nb[bi] = len(pooled_means)
        mean_est[bi] = pooled_means.mean()
        if len(pooled_means) >= 2:
            se[bi] = np.std(pooled_means, ddof=1) / np.sqrt(len(pooled_means))

    start_idx = ref_idx = None
    valid = np.where(pooled_nb >= min_pooled_blocks)[0]
    if len(valid) > 0:
        ref_idx = int(valid[-1])
        ref_se = se[ref_idx]
        start_idx = ref_idx
        if ref_se > 0 and np.isfinite(ref_se):
            for i in valid[::-1]:
                if np.isfinite(se[i]) and abs(se[i] - ref_se) / ref_se <= plateau_rel_tol:
                    start_idx = int(i)
                else:
                    break

    result = {'block_ns': block_ns, 'se': se, 'mean_est': mean_est,
              'pooled_nb': pooled_nb, 'start_idx': start_idx, 'ref_idx': ref_idx}
    if start_idx is not None:
        w = pooled_nb[start_idx:ref_idx + 1].astype(float)
        result['plateau_se'] = float(np.average(se[start_idx:ref_idx + 1], weights=w))
        result['plateau_mean'] = float(np.average(mean_est[start_idx:ref_idx + 1], weights=w))
    return result


def plot_block_averaging(pooled, sets, pair_keys, metric):
    """Per-pair block-averaging figure (one per set): top row = pooled SE vs
    block size, bottom row = pooled mean distance vs block size, shaded band
    = the plateau region each pair's own SE curve settles into (independently
    detected per pair -- see block_average_plateau)."""
    metric_label = 'COM-COM' if metric == 'com' else 'min CA-CA'
    for sk in sets:
        n_pairs = len(pair_keys)
        fig, axes = plt.subplots(2, n_pairs, figsize=(4.5 * n_pairs, 7),
                                  squeeze=False, sharex='col')
        any_plotted = False
        for col, pk in enumerate(pair_keys):
            ax_se, ax_mean = axes[0, col], axes[1, col]
            a, b = pk.split('_')
            rep_arrays = pooled.get(sk, {}).get(pk, [])
            res = block_average_plateau(rep_arrays)
            ax_se.set_title(f'{a}–{b}', fontsize=11, fontweight='bold')
            if res is None:
                continue
            any_plotted = True
            ax_se.plot(res['block_ns'], res['se'], lw=2, marker='o', ms=3, color='#1f77b4')
            ax_se.set_xscale('log')
            ax_mean.plot(res['block_ns'], res['mean_est'], lw=2, marker='o', ms=3, color='#1f77b4')
            ax_mean.set_xscale('log')
            ax_mean.set_xlabel('Block size (ns)')
            if res['start_idx'] is not None:
                lo = res['block_ns'][res['start_idx']]
                hi = res['block_ns'][res['ref_idx']]
                for ax in (ax_se, ax_mean):
                    ax.axvspan(lo, hi, color='red', alpha=0.12, zorder=0)
                ax_se.axhline(res['plateau_se'], color='red', ls='--', lw=1.2, alpha=0.8,
                              label=f"{res['plateau_mean']:.2f} +/- {res['plateau_se']:.3f} nm\n"
                                    f"(plateau {lo:.0f}-{hi:.0f} ns)")
                ax_se.legend(fontsize=7, loc='best')
                ax_mean.axhline(res['plateau_mean'], color='red', ls='--', lw=1.2, alpha=0.8)
        axes[0, 0].set_ylabel(f'Pooled block SE\n({metric_label}, nm)')
        axes[1, 0].set_ylabel(f'Pooled mean dist.\n({metric_label}, nm)')
        if not any_plotted:
            plt.close(fig)
            continue
        fig.suptitle(f'{sk}: Inter-Domain {metric_label} Distance — Block-Averaging, '
                     f'Pooled Across Replicates\n(top: SE vs block size; bottom: mean distance '
                     f'vs block size; shaded = plateau region, detected independently per pair)',
                     fontsize=11, y=1.03)
        fig.tight_layout()
        out_png = FIG_PATH / f'md_distance_block_{sk}_{metric}.png'
        out_svg = FIG_PATH / f'md_distance_block_{sk}_{metric}.svg'
        fig.savefig(out_png, dpi=150, bbox_inches='tight')
        fig.savefig(out_svg, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved: {out_png}")


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--set', nargs='+', default=['fl_optimized'],
                        help='Sets to include (e.g. fl_optimized, fl, '
                             'or fragment names)')
    parser.add_argument('--metric', choices=['com', 'min'], default='com',
                        help='Distance metric: COM-COM (default) or min CA-CA')
    parser.add_argument('--workers', type=int, default=8,
                        help='Parallel workers (default 8)')
    parser.add_argument('--pairs', nargs='+', default=None,
                        help='Pairs to plot (e.g. MD1L1_MD3 MD1L1_ART). '
                             'Default: all 5 pairs')
    parser.add_argument('--sim-folder', nargs='+', default=None, metavar='PATH',
                        help='One or more NEW simulation-set folders to analyze '
                             '(each containing seed-*_sample-*/ replicates with '
                             'top.pdb + <sysname>.dcd). Domain units are read from '
                             "the folder's metadata.json, or pass --units.")
    parser.add_argument('--units', nargs='+', default=None, metavar='UNIT',
                        help='FL domain units in the --sim-folder construct '
                             '(e.g. md1l1 md2 md3). Required only if the folder '
                             'has no metadata.json.')
    parser.add_argument('--sim-folder-as', nargs='+', default=None, metavar='NAME:PATH',
                        help='Register a sim folder under a CUSTOM set_key '
                             '(e.g. fl_go_5rep1us:/path/to/fl_go). Lets the same '
                             'physical folder be registered twice under two '
                             'different set_keys (e.g. once with the full data, '
                             'once with --max-reps/--max-ns applied) so both can '
                             'be compared side by side in one run.')
    parser.add_argument('--max-reps', nargs='+', default=None, metavar='NAME:N',
                        help='Cap a set to its first N replicates, in (seed, '
                             'sample) enumeration order (e.g. fl_go_5rep1us:5).')
    parser.add_argument('--max-ns', nargs='+', default=None, metavar='NAME:N',
                        help='Cap a set to the first N ns per replicate, after '
                             'equilibration skip (e.g. fl_go_5rep1us:1000).')
    parser.add_argument('--only-sample', nargs='+', default=None, metavar='NAME:N',
                        help='Keep only replicates with this sample index (e.g. '
                             'md_full:0 -- md_full only extended seed-{1..5}_sample-0 '
                             'to the full 1000 ns target, sample-{1,2,3,4} are still '
                             'at the original ~100 ns; --max-reps alone would wrongly '
                             'grab seed-1 samples 0-4 instead of the 5 actually-'
                             'extended replicates).')
    args = parser.parse_args()

    # Register any --sim-folder folders before building jobs / launching workers
    # so forked workers inherit EXTERNAL_DIRS / EXTERNAL_SYSNAME / CONSTRUCT_UNITS.
    external_keys = []
    if args.sim_folder:
        for folder in args.sim_folder:
            k = register_sim_folder(folder, units=args.units)
            external_keys.append(k)
            print(f"Registered --sim-folder '{folder}' as set '{k}' "
                  f"(sysname: {EXTERNAL_SYSNAME[k]}, "
                  f"units: {', '.join(CONSTRUCT_UNITS[k])})")
    if args.sim_folder_as:
        for spec in args.sim_folder_as:
            name, folder = spec.split(':', 1)
            k = register_sim_folder(folder, units=args.units, name=name)
            external_keys.append(k)
            print(f"Registered --sim-folder-as '{folder}' as set '{k}' "
                  f"(sysname: {EXTERNAL_SYSNAME[k]}, "
                  f"units: {', '.join(CONSTRUCT_UNITS[k])})")
    if external_keys:
        # Replace the default set with the folders unless --set was given explicitly.
        if args.set == ['fl_optimized']:
            args.set = external_keys
        else:
            args.set = args.set + [k for k in external_keys if k not in args.set]

    if args.max_reps:
        for spec in args.max_reps:
            name, n = spec.split(':', 1)
            REPLICATE_CAP[name] = int(n)
            print(f"Capping set '{name}' to first {n} replicates")
    if args.max_ns:
        for spec in args.max_ns:
            name, n = spec.split(':', 1)
            FRAME_CAP[name] = int(float(n) * FRAMES_PER_NS)
            print(f"Capping set '{name}' to first {n} ns/replicate "
                  f"({FRAME_CAP[name]} frames after skip)")
    if args.only_sample:
        for spec in args.only_sample:
            name, n = spec.split(':', 1)
            ONLY_SAMPLE[name] = int(n)
            print(f"Restricting set '{name}' to sample-{n} replicates only")

    sets = args.set
    metric = args.metric

    global FIG_PATH
    FIG_PATH = _get_fig_dir('04_md_distances', sims=sets)

    # Filter pairs
    if args.pairs:
        wanted = set(args.pairs)
        pair_keys = [f'{a}_{b}' for a, b in PAIRS if f'{a}_{b}' in wanted]
    else:
        pair_keys = [f'{a}_{b}' for a, b in PAIRS]

    print(f"Sets:       {sets}")
    print(f"Metric:     {metric} ({'COM-COM' if metric=='com' else 'min CA-CA'})")
    print(f"Pairs:      {pair_keys}")
    print(f"Workers:    {args.workers}")

    # Build job list. ONLY_SAMPLE[sk], if set, drops every (seed, sample) pair
    # whose sample doesn't match -- applied BEFORE REPLICATE_CAP so the two
    # compose (e.g. only-sample to pick the extended replicates, then
    # max-reps if you also want fewer than all matching seeds). REPLICATE_CAP,
    # if set, keeps only the first N surviving (seed, sample) pairs in
    # enumeration order -- for a --sim-folder flat-dir set (e.g. fl_go) that's
    # the first N dirs in sorted order (idx 0..N-1).
    jobs = []
    for sk in sets:
        only_sample = ONLY_SAMPLE.get(sk)
        cap = REPLICATE_CAP.get(sk)
        n = 0
        for seed in SEEDS:
            for sample in SAMPLES:
                if only_sample is not None and sample != only_sample:
                    continue
                if cap is not None and n >= cap:
                    break
                jobs.append((sk, seed, sample, metric))
                n += 1
            if cap is not None and n >= cap:
                break

    # Run in parallel
    print(f"\nProcessing {len(jobs)} replicates...")
    pooled = {sk: {pk: [] for pk in pair_keys} for sk in sets}
    pooled_rg = {sk: {} for sk in sets}  # sk -> domain -> list of per-replicate Rg arrays

    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(compute_distances_one_replicate, j): j
                   for j in jobs}
        for fut in as_completed(futures):
            res = fut.result()
            if res is None:
                continue
            sk = res['set_key']
            for k, arr in res['distances'].items():
                if k in pooled[sk]:
                    pooled[sk][k].append(arr)
            for dname, arr in res.get('domain_rg', {}).items():
                pooled_rg[sk].setdefault(dname, []).append(arr)

    # Pool all replicate arrays into one big array per pair per set
    pooled_arrays = {}
    for sk in sets:
        pooled_arrays[sk] = {}
        for pk in pair_keys:
            if pooled[sk][pk]:
                pooled_arrays[sk][pk] = np.concatenate(pooled[sk][pk])
                n_frames = len(pooled_arrays[sk][pk])
                print(f"  {sk} {pk}: {n_frames} frames, "
                      f"mean={pooled_arrays[sk][pk].mean():.2f} nm, "
                      f"std={pooled_arrays[sk][pk].std():.2f}")
            else:
                pooled_arrays[sk][pk] = np.array([])

    # Per-domain mean Rg, per set -- sets a physically-motivated COM-COM
    # contact floor (Rg_A + Rg_B, "surfaces touching") that scales with each
    # domain's actual size, unlike a single fixed cutoff tuned for min-CA-CA
    # (residue-scale) contact. See compute_distances_one_replicate docstring.
    mean_rg = {sk: {} for sk in sets}
    for sk in sets:
        for dname, arrs in pooled_rg[sk].items():
            if arrs:
                mean_rg[sk][dname] = float(np.concatenate(arrs).mean())
        if mean_rg[sk]:
            rg_str = ', '.join(f'{d}={v:.2f}' for d, v in mean_rg[sk].items())
            print(f"  {sk} domain Rg (nm): {rg_str}")

    # Save raw data + per-domain mean Rg (so other scripts, e.g.
    # compare_full_sims_grid.py, can reuse both without recomputation).
    # Per-replicate arrays (pooled[sk][pk], pre-concatenation) are ALSO saved
    # individually -- needed to compute per-replicate contact fractions with
    # a real across-replicate error estimate (see
    # compare_full_sims_contact_bars.py); the concatenated array above loses
    # replicate boundaries entirely.
    out_npz = DATA_PATH / f'md_distances_{metric}.npz'
    save_dict = {}
    for sk in sets:
        for pk in pair_keys:
            save_dict[f'{sk}__{pk}'] = pooled_arrays[sk][pk]
            save_dict[f'{sk}__{pk}__nreps'] = np.array([len(pooled[sk][pk])])
            for i, rep_arr in enumerate(pooled[sk][pk]):
                save_dict[f'{sk}__{pk}__rep{i}'] = rep_arr
        for dname, rg in mean_rg[sk].items():
            save_dict[f'{sk}__rg__{dname}'] = np.array([rg])
    np.savez(out_npz, **save_dict)
    print(f"\n  Saved: {out_npz}")

    # ── Violin plot ──
    plot_violin(pooled_arrays, sets, pair_keys, metric, mean_rg)

    # ── Block-averaging (per-replicate arrays, NOT the concatenated ones) ──
    plot_block_averaging(pooled, sets, pair_keys, metric)

    # ── Per-pair histogram/CDF + bimodality (antimode) detection ──
    plot_distance_histograms(pooled_arrays, sets, pair_keys, metric, mean_rg)


# ============================================================
# Plotting
# ============================================================

def plot_violin(pooled_arrays, sets, pair_keys, metric, mean_rg=None):
    """Massive violin plot: one violin per (set, pair) combination.

    Y-axis: inter-domain distance (nm)
    X-axis: pair, grouped by set with side-by-side comparison
    Violin shape: frequency density across all frames
    Overlay: scatter of all data points (low alpha)

    mean_rg: {set: {domain: mean Rg (nm)}}, optional. For metric='com', a
    fixed CONTACT_CUTOFF_NM (tuned for min-CA-CA, residue-scale contact) is
    meaningless -- two folded, sterically-excluding domains' COMs can't get
    closer than roughly Rg_A + Rg_B (surfaces touching) regardless of any
    "contact" definition, and that floor scales with domain size. When
    mean_rg is given, a per-pair Rg_A+Rg_B reference line is drawn instead
    of (COM) / alongside (min) the fixed cutoff.
    """
    # Color per set -- sim_registry.SETS is the single source of truth (fixes
    # a stale local copy that only had colors for the original 5 sets, so
    # every one of the 13 newer "contiguous family" sets silently fell
    # through to an unvalidated tab10-cycled default).
    set_colors = _set_color_lookup(sets)
    default_palette = plt.cm.tab10(np.linspace(0, 1, len(sets)))

    n_pairs = len(pair_keys)
    n_sets = len(sets)
    width_per_pair = 1.5
    fig_width = max(8, n_pairs * width_per_pair * max(n_sets, 1) + 2)
    fig, ax = plt.subplots(figsize=(fig_width, 8))

    positions = []
    labels_minor = []
    labels_major = []
    label_positions = []

    pos = 0
    pair_centers = []
    rg_floor_segments = []    # (xmin, xmax, Rg_A+Rg_B) per pair -- "surfaces touching"
    rg_energetic_segments = []  # (xmin, xmax, Rg_A+Rg_B+AH_CUTOFF_NM) -- "AH could still reach"
    for pk in pair_keys:
        a, b = pk.split('_')
        pair_start = pos
        for i, sk in enumerate(sets):
            arr = pooled_arrays[sk][pk]
            if len(arr) == 0:
                pos += 1
                continue
            color = set_colors.get(sk, default_palette[i])
            parts = ax.violinplot([arr], positions=[pos], widths=0.8,
                                  showmeans=True, showmedians=True,
                                  showextrema=False)
            for body in parts['bodies']:
                body.set_facecolor(color)
                body.set_alpha(0.55)
                body.set_edgecolor('black')
                body.set_linewidth(0.5)
            for key in ('cmeans', 'cmedians'):
                if key in parts:
                    parts[key].set_color('black')
                    parts[key].set_linewidth(1.2)

            # Scatter overlay (subsampled for performance)
            rng = np.random.default_rng(42)
            sub = arr if len(arr) < 5000 else rng.choice(arr, 5000, replace=False)
            jitter = rng.uniform(-0.18, 0.18, len(sub))
            ax.scatter(pos + jitter, sub, s=1.5, c=color, alpha=0.15,
                       edgecolors='none', zorder=2)

            positions.append(pos)
            labels_minor.append(sk if n_sets > 1 else '')
            pos += 1
        # Center label for the pair
        pair_centers.append((pair_start + pos - 1) / 2)
        labels_major.append(f'{a}–{b}')
        if mean_rg is not None and pos > pair_start:
            rg_sums = [mean_rg[sk][a] + mean_rg[sk][b] for sk in sets
                       if a in mean_rg.get(sk, {}) and b in mean_rg.get(sk, {})]
            if rg_sums:
                rg_sum = float(np.mean(rg_sums))
                xmin, xmax = pair_start - 0.4, pos - 1 + 0.4
                rg_floor_segments.append((xmin, xmax, rg_sum))
                rg_energetic_segments.append((xmin, xmax, rg_sum + AH_CUTOFF_NM))
        pos += 0.8  # gap between pairs

    # Contact / interaction-range reference. min-CA-CA vs. AH_CUTOFF_NM is a
    # direct, physically correct comparison (both are bead-bead, residue-
    # scale distances -- see figure header). COM-COM needs a per-pair,
    # domain-size-scaled reference instead of one fixed line (see
    # compute_distances_one_replicate docstring): the excluded-volume floor
    # Rg_A+Rg_B ("surfaces touching") and Rg_A+Rg_B+AH_CUTOFF_NM ("AH could
    # still reach between peripheral residues even with COMs this far
    # apart") bound the COM-COM range where the domains COULD be
    # energetically coupled.
    if metric == 'com' and rg_floor_segments:
        for xmin, xmax, y in rg_energetic_segments:
            ax.hlines(y, xmin, xmax, color='mediumpurple', linestyle='--',
                      linewidth=1.8, alpha=0.85, zorder=4)
        for xmin, xmax, y in rg_floor_segments:
            ax.hlines(y, xmin, xmax, color='darkorange', linestyle='--',
                      linewidth=1.8, alpha=0.85, zorder=4)
    else:
        ax.axhline(AH_CUTOFF_NM, color='red', linestyle='--', alpha=0.6,
                   linewidth=1.5,
                   label=f'AH cutoff ({AH_CUTOFF_NM:.1f} nm)')

    # X-axis: minor ticks for sets (if multi-set), major ticks for pairs
    if n_sets > 1:
        ax.set_xticks(positions, minor=True)
        ax.set_xticklabels(labels_minor, minor=True, rotation=45, ha='right',
                           fontsize=7)
        ax.set_xticks(pair_centers, minor=False)
        ax.set_xticklabels(labels_major, minor=False, fontsize=12,
                           fontweight='bold')
        ax.tick_params(axis='x', which='major', pad=22)
    else:
        ax.set_xticks(pair_centers)
        ax.set_xticklabels(labels_major, fontsize=12, fontweight='bold')

    if metric == 'com':
        metric_label = 'COM–COM distance'
        axis_label = 'COM–COM distance (nm)\nreference: Rg_A + Rg_B  and  Rg_A + Rg_B + 2 nm'
    else:
        metric_label = 'minimum CA–CA distance'
        axis_label = 'minimum CA–CA distance (nm)\nreference: AH cutoff (2.0 nm)'
    ax.set_ylabel(axis_label, fontsize=12)
    n_frames_per = '\n'.join(
        f'{sk}: {sum(len(pooled_arrays[sk][pk]) for pk in pair_keys)} pooled frames'
        for sk in sets)
    ax.set_title(
        f'Inter-Domain {metric_label} Distributions [metric={metric}]\n'
        f'{n_frames_per}',
        fontsize=12)

    # Legend with set colors
    legend_elements = []
    for sk in sets:
        color = set_colors.get(sk, '#666')
        legend_elements.append(Line2D([0], [0], marker='s', color='w',
                                       markerfacecolor=color, markersize=12,
                                       label=sk))
    if metric == 'com' and rg_floor_segments:
        legend_elements.append(Line2D([0], [0], color='darkorange', linestyle='--',
                                       label='Rg_A + Rg_B\n(surfaces touching)'))
        legend_elements.append(Line2D([0], [0], color='mediumpurple', linestyle='--',
                                       label=f'Rg_A + Rg_B + {AH_CUTOFF_NM:.0f}\n(AH could still reach)'))
    else:
        legend_elements.append(Line2D([0], [0], color='red', linestyle='--',
                                       label=f'AH cutoff ({AH_CUTOFF_NM:.1f} nm)'))
    ax.legend(handles=legend_elements, fontsize=10, loc='upper right')

    ax.grid(axis='y', alpha=0.3)
    ax.set_ylim(bottom=0)

    plt.tight_layout()
    sets_str = '_'.join(sets)
    out_png = FIG_PATH / f'md_distance_violin_{sets_str}_{metric}.png'
    out_svg = FIG_PATH / f'md_distance_violin_{sets_str}_{metric}.svg'
    fig.savefig(out_png, dpi=200, bbox_inches='tight')
    fig.savefig(out_svg, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_png}")
    print(f"  Saved: {out_svg}")


def detect_bimodality(arr, bw_method=0.05, min_valley_depth=0.05):
    """Look for a genuine antimode (valley between two peaks) in a KDE of
    `arr`, rather than trusting a GMM component count: with a 2-component
    Gaussian mixture, BIC kept improving all the way to 3+ components for
    every pair tested here (consistent with this being a continuous,
    diffusive distribution -- no discrete states, matching the TICA
    landscape finding elsewhere in this project), so "which n minimizes
    BIC" is not a trustworthy bimodality test. A real valley between two
    well-separated local maxima is a much more direct, harder-to-fake
    signal, and it's what "does this look bimodal by eye" actually means.

    Returns dict with is_bimodal (bool), peaks (nm), valley (nm or None),
    and depth (fractional drop from the smaller peak to the valley; only
    meaningful if is_bimodal). min_valley_depth: minimum fractional drop
    to call it real bimodality rather than KDE/sampling noise.
    """
    from scipy.stats import gaussian_kde
    from scipy.signal import argrelextrema
    if len(arr) < 100:
        return {'is_bimodal': False, 'peaks': [], 'valley': None, 'depth': 0.0}
    kde = gaussian_kde(arr, bw_method=bw_method)
    xs = np.linspace(arr.min(), arr.max(), 400)
    ys = kde(xs)
    # Ignore near-zero-density tail wiggles (KDE noise far from the bulk).
    keep = ys > 0.01 * ys.max()
    maxima = [i for i in argrelextrema(ys, np.greater)[0] if keep[i]]
    minima = [i for i in argrelextrema(ys, np.less)[0] if keep[i]]
    if len(maxima) < 2 or not minima:
        peaks = xs[maxima].tolist() if maxima else [float(xs[np.argmax(ys)])]
        return {'is_bimodal': False, 'peaks': peaks, 'valley': None, 'depth': 0.0}
    # Take the two tallest peaks and the deepest valley strictly between them.
    top2 = sorted(maxima, key=lambda i: -ys[i])[:2]
    lo, hi = sorted(top2)
    between = [i for i in minima if lo < i < hi]
    if not between:
        return {'is_bimodal': False, 'peaks': xs[top2].tolist(), 'valley': None, 'depth': 0.0}
    v = min(between, key=lambda i: ys[i])
    depth = 1.0 - ys[v] / min(ys[lo], ys[hi])
    return {'is_bimodal': bool(depth >= min_valley_depth),
            'peaks': sorted(xs[[lo, hi]].tolist()),
            'valley': float(xs[v]), 'depth': float(depth)}


def plot_distance_histograms(pooled_arrays, sets, pair_keys, metric, mean_rg=None):
    """Per-pair histogram + cumulative distribution, with a per-pair
    Rg_A+Rg_B excluded-volume floor (metric='com') instead of a fixed
    cutoff, and a bimodality (antimode) check annotated on each histogram.
    """
    set_colors = _set_color_lookup(sets)
    default_palette = plt.cm.tab10(np.linspace(0, 1, len(sets)))

    n_pairs = len(pair_keys)
    fig, axes = plt.subplots(2, n_pairs, figsize=(5 * n_pairs, 7),
                              squeeze=False)

    # Find global x range. Can be empty -- e.g. a construct with only ONE of
    # the 4 tracked domains (MD1L1/MD2/MD3/ART) present (md3_wwe_full only
    # has MD3) has zero valid pairs among pair_keys, so every array here is
    # length-0 and np.concatenate([]) would crash; fall back to a default
    # xmax and skip plotting below instead.
    all_dists = []
    for sk in sets:
        for pk in pair_keys:
            all_dists.append(pooled_arrays[sk][pk])
    nonempty = [a for a in all_dists if len(a) > 0]
    if not nonempty:
        print("  No pairs have data for this set selection (fewer than 2 of "
              "MD1L1/MD2/MD3/ART present) -- skipping histogram/CDF plot.")
        return
    all_combined = np.concatenate(nonempty)
    xmax = np.percentile(all_combined, 99.5) * 1.05 if len(all_combined) else 10

    print("\n  Bimodality check (KDE antimode; not a discrete-state claim -- "
          "see detect_bimodality docstring):")
    for col, pk in enumerate(pair_keys):
        a, b = pk.split('_')
        ax_hist = axes[0, col]
        ax_cdf = axes[1, col]
        rg_floor = rg_energetic = None
        if metric == 'com' and mean_rg is not None:
            sums = [mean_rg[sk][a] + mean_rg[sk][b] for sk in sets
                    if a in mean_rg.get(sk, {}) and b in mean_rg.get(sk, {})]
            if sums:
                rg_floor = float(np.mean(sums))
                rg_energetic = rg_floor + AH_CUTOFF_NM

        for i, sk in enumerate(sets):
            arr = pooled_arrays[sk][pk]
            if len(arr) == 0:
                continue
            color = set_colors.get(sk, default_palette[i])
            ax_hist.hist(arr, bins=80, density=True, alpha=0.5,
                         color=color, label=sk, range=(0, xmax))
            # CDF
            sorted_arr = np.sort(arr)
            cdf = np.arange(1, len(sorted_arr) + 1) / len(sorted_arr)
            ax_cdf.plot(sorted_arr, cdf, color=color, label=sk, linewidth=1.5)

            if rg_floor is not None:
                for cutoff, cutoff_color in ((rg_floor, 'darkorange'),
                                             (rg_energetic, 'mediumpurple')):
                    ax_cdf.axvline(cutoff, color=cutoff_color, linestyle='--', alpha=0.5, linewidth=1)
                    ax_hist.axvline(cutoff, color=cutoff_color, linestyle='--', alpha=0.5, linewidth=1)
                frac_floor = np.mean(arr < rg_floor)
                frac_energetic = np.mean(arr < rg_energetic)
                ax_cdf.text(0.98, 0.05 + i * 0.10,
                            f'{sk}: {frac_floor*100:.1f}% < RgA+RgB, '
                            f'{frac_energetic*100:.1f}% < RgA+RgB+{AH_CUTOFF_NM:.0f}',
                            transform=ax_cdf.transAxes, fontsize=8,
                            ha='right', color=color, fontweight='bold')
            else:
                cutoff, cutoff_color = AH_CUTOFF_NM, 'red'
                frac_contact = np.mean(arr < cutoff)
                ax_cdf.axvline(cutoff, color=cutoff_color, linestyle='--', alpha=0.5, linewidth=1)
                ax_hist.axvline(cutoff, color=cutoff_color, linestyle='--', alpha=0.5, linewidth=1)
                ax_cdf.text(0.98, 0.05 + i * 0.06,
                            f'{sk}: {frac_contact*100:.1f}% < AH cutoff ({AH_CUTOFF_NM:.1f} nm)',
                            transform=ax_cdf.transAxes, fontsize=8,
                            ha='right', color=color, fontweight='bold')

            # Bimodality (antimode) check, pooled-across-replicates array.
            bm = detect_bimodality(arr)
            if bm['is_bimodal']:
                for p in bm['peaks']:
                    ax_hist.axvline(p, color='black', linestyle=':', linewidth=1, alpha=0.6)
                ax_hist.axvline(bm['valley'], color='green', linestyle='-', linewidth=1.5, alpha=0.8)
                verdict = (f"BIMODAL: peaks {bm['peaks'][0]:.2f}/{bm['peaks'][1]:.2f} nm, "
                           f"valley {bm['valley']:.2f} nm, depth {bm['depth']*100:.0f}%")
            else:
                pk_str = '/'.join(f'{p:.2f}' for p in bm['peaks'])
                verdict = f"unimodal (peak {pk_str} nm)" if bm['peaks'] else "unimodal"
            print(f"    {sk} {a}-{b}: {verdict}")

        xlabel = ('COM–COM distance (nm)' if metric == 'com'
                  else 'minimum CA–CA distance (nm)')
        ax_hist.set_title(f'{a}–{b}', fontsize=12, fontweight='bold')
        ax_hist.set_xlabel(xlabel, fontsize=10)
        ax_hist.set_ylabel('Density', fontsize=10)
        ax_hist.legend(fontsize=8, loc='upper right')
        ax_hist.set_xlim(0, xmax)

        ax_cdf.set_xlabel(xlabel, fontsize=10)
        ax_cdf.set_ylabel('Cumulative fraction', fontsize=10)
        ax_cdf.set_xlim(0, xmax)
        ax_cdf.set_ylim(0, 1.02)
        ax_cdf.grid(alpha=0.3)

    if metric == 'com':
        metric_label = 'COM–COM'
        cutoff_note = ('orange = Rg_A+Rg_B "surfaces touching"; '
                       'purple = Rg_A+Rg_B+2 "AH could still reach"')
    else:
        metric_label = 'minimum CA–CA'
        cutoff_note = 'red = AH cutoff (2.0 nm)'
    fig.suptitle(f'Inter-Domain {metric_label} Distance Distributions [metric={metric}]\n'
                 f'({cutoff_note}; black dotted + green = detected bimodal peaks/valley)',
                 fontsize=12, y=1.03)
    plt.tight_layout()
    sets_str = '_'.join(sets)
    out_png = FIG_PATH / f'md_distance_hist_{sets_str}_{metric}.png'
    fig.savefig(out_png, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_png}")


if __name__ == '__main__':
    main()
