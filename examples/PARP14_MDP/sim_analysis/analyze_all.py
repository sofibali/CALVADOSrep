#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
General PARP14 CALVADOS trajectory analysis across all simulation sets.
Parallelized with ProcessPoolExecutor (--workers N).

Integrates CALVADOS built-in analysis (calc_rg, calc_ete, calc_dmap, calc_cmap,
cmap_traj, calc_fnc, calc_wcn, calc_energy, fit_scaling_exp) with custom
SAA accessibility analysis.

Analyses (each writes its own data file to data/ and figures to figures/02_main_analysis/<date>/):
  1. Conformational properties: Rg, Ree (via CALVADOS calc_rg/calc_ete)
  2. Distance maps: interdomain COM distance matrices
  3. Contact maps: trajectory-averaged inter-domain contact frequency (cmap_traj)
  4. FNC: fraction of native contacts per restrained domain (calc_fnc)
  5. Energy decomposition: AH + Yukawa between domain pairs (calc_energy)
  6. WCN: weighted coordination number at active sites (calc_wcn)
  7. Active site analysis: inter-site distances, radial position, RMSF
  8. Accessibility: SAA ray-casting, cone angle, shell density (unique)

Recognized simulation sets:
  fl              — full-length (1801 res, EBI AF2)
  fl_optimized    — full-length with optimized per-domain restraint trims
                    + KH7a–KHb custom inter-domain restraints
  md              — Macrodomains only (MD1L1+MD2+MD3, 586 res)
  core            — KH7a+MD1L1+MD2+MD3+KHb-KH8+WWE+ART (1051 res)
  mka             — MD1L1+MD2+MD3+KHb-KH8+WWE+ART (999 res)
  norrm           — KH1-6+KH7a+MDs+KHb-KH8+WWE+ART (1474 res)
  noart           — KH1-6+KH7a+MDs+KHb-KH8+WWE (1275 res)
  md3art          — MD3+KHb-KH8+WWE+ART (595 res)

Fragments (--include-fragments or --set all):
  Auto-discovers ALL prepared contiguous fragments under fragments/<name>/.
  Each fragment lives in `fragments/<units_joined>/seed-{1-5}_sample-{0-4}/`
  with a `metadata.json` describing units, restraints, and box size. They
  are registered internally as `frag_<name>` to avoid collisions with the
  named sets above.

  Fragment naming convention: lowercase domain units joined by '_', e.g.
  `md1l1_md2_md3`, `kh7a_md1l1_md2_md3`, `khb-kh8_wwe_art`.

Data caching (avoid re-computing expensive trajectory analyses):
  Each module saves results to data/<set>_<analysis>.npy/.npz. By default,
  if the data files exist for a set, that module SKIPS the parallel
  trajectory compute and rebuilds figures from the cached data. Use
  --force / --force-recompute to overwrite (regenerate the data from
  trajectories). Note: --force re-runs the slow MDAnalysis loop.

Usage:
    python analyze_all.py                            # everything, all known sets
    python analyze_all.py --workers 10               # limit parallel workers
    python analyze_all.py --conf-prop                # conformational properties only
    python analyze_all.py --dmap                     # distance maps only
    python analyze_all.py --cmap                     # contact maps only
    python analyze_all.py --fnc                      # FNC only
    python analyze_all.py --energy                   # energy decomposition only
    python analyze_all.py --wcn                      # WCN at active sites
    python analyze_all.py --active-sites             # active site analysis
    python analyze_all.py --accessibility            # SAA accessibility
    python analyze_all.py --set fl_optimized md      # specific sets
    python analyze_all.py --set md1l1_md2_md3        # specific fragment
    python analyze_all.py --include-fragments        # auto-discover all fragments
    python analyze_all.py --set all                  # all sets + all fragments
    python analyze_all.py --force                    # force data recompute
    python analyze_all.py --cmap --force             # force one module's recompute
"""

import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import MDAnalysis as mda
import os
import warnings
from argparse import ArgumentParser
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

warnings.filterwarnings('ignore')

# CALVADOS analysis imports
from calvados.analysis import (
    calc_rg, calc_ete, calc_dmap, calc_cmap, cmap_traj,
    calc_fnc, calc_wcn, calc_energy, fit_scaling_exp,
)

# ============================================================
# Configuration
# ============================================================

CWD = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_PATH = os.path.join(CWD, 'data')
# Dated category subdir: figures/02_main_analysis/<YYYY-MM-DD>/
import sys as _sys
_sys.path.insert(0, CWD)
import os as _os_boot, sys as _sys_boot  # __ROOTBOOT__ (script lives in a subfolder; root=parent)
_sys_boot.path.insert(0, _os_boot.path.dirname(_os_boot.path.dirname(_os_boot.path.abspath(__file__))))
from _fig_layout import get_fig_dir as _get_fig_dir
FIG_PATH = str(_get_fig_dir('02_main_analysis'))

RESIDUES_FILE = os.path.join(CWD, 'input', 'residues_CALVADOS3.csv')

SEEDS = list(range(1, 6))
SAMPLES = list(range(0, 5))
SKIP_FRAMES = 50  # 0.5 ns equilibration

# Simulation parameters
# NOTE ON UNITS: CALVADOS works entirely in nanometers (nm). MDAnalysis, however,
# stores trajectory coordinates in Angstroms (Å). Every `positions / 10.0` and
# `center_of_mass() / 10.0` below is the Å -> nm conversion (1 nm = 10 Å). So all
# distances, COMs, radii, and energy cutoffs in this script are in nm unless noted.
EPS_LJ = 0.2       # Ashbaugh-Hatch well depth epsilon (kJ/mol)
CUTOFF_LJ = 2.0    # Ashbaugh-Hatch (hydrophobic) interaction cutoff distance = 2.0 nm
CUTOFF_YU = 4.0    # Yukawa (electrostatic) interaction cutoff distance = 4.0 nm
IONIC = 0.19       # ionic strength (mol/L), sets Debye screening length
TEMP = 293         # temperature (K)

# Set definitions
SETS = {
    'fl':    {'sysname': 'parp14',              'label': 'Full Length (1-1801)',    'color': '#1f77b4'},
    'norrm': {'sysname': 'parp14_norrm',        'label': 'KH1-ART (315-1801)',    'color': '#9467bd'},
    'noart': {'sysname': 'parp14_noart',        'label': 'KH1-WWE (315-1602)',    'color': '#8c564b'},
    'core':  {'sysname': 'parp14_core',         'label': 'KH7-ART (738-1801)',    'color': '#2ca02c'},
    'mka':   {'sysname': 'parp14_mka',          'label': 'MD1-ART (790-1801)',    'color': '#d62728'},
    'md':    {'sysname': 'parp14_macrodomains', 'label': 'MD1-MD3 (790-1388)',    'color': '#ff7f0e'},
    'md3art':{'sysname': 'parp14_md3art',      'label': 'MD3-ART (1207-1801)',   'color': '#17becf'},
    'fl_optimized': {'sysname': 'parp14', 'label': 'FL (optimized restraints)', 'color': '#aec7e8'},
}

CONSTRUCT_DOMAINS = {
    'fl':    {'parp14':              [[6,88],[150,223],[227,301],[791,978],[1003,1190],[1216,1387],[1523,1601],[1605,1801]]},
    'md':    {'parp14_macrodomains': [[2,189],[214,401],[414,585]]},
    'core':  {'parp14_core':         [[54,241],[266,453],[466,637],[773,851],[855,1051]]},
    'mka':   {'parp14_mka':          [[2,189],[214,401],[414,585],[721,799],[803,999]]},
    'norrm': {'parp14_norrm':        [[477,664],[689,876],[889,1060],[1196,1274],[1278,1474]]},
    'noart': {'parp14_noart':        [[477,664],[689,876],[889,1060],[1196,1274]]},
    'md3art':{'parp14_md3art':       [[10,181],[317,395],[399,595]]},
    # fl_optimized: 17 restraint domains with per-domain trim values, + KH7a-KHb custom pairs
    'fl_optimized': {'parp14': [[16,78],[156,214],[235,304],[320,379],[390,449],[460,515],
                                 [531,583],[599,660],[676,727],[738,789],[800,968],
                                 [1015,1183],[1217,1378],[1389,1461],[1462,1533],
                                 [1549,1587],[1613,1791]]},
}

# Energy-analysis domain boundaries: trimmed by 3 residues at zero/small-gap
# boundaries to prevent steric-clash artifacts between adjacent domains.
# Original (restraint) boundaries are in FL_RESTRAINT_DOMAINS in prepare_and_run_all.py.
FL_DOMAINS = {
    'RRM1':    (6, 88),
    'RRM2':    (150, 223),
    'RRM3':    (227, 301),
    'KH1-6':   (315, 734),   # trimmed end by 3 (gap 0 to KH7a)
    'KH7a':    (741, 786),   # trimmed both ends by 3 (gap 0 from KH1-6, gap 1 to MD1)
    'MD1':     (794, 978),   # trimmed start by 3 (gap 1 from KH7a)
    'MD2':     (1003, 1190),
    'MD3':     (1216, 1384), # trimmed end by 3 (gap 1 to KHb-KH8)
    'KHb-KH8': (1392, 1530), # trimmed both ends by 3 (gap 1 from MD3, gap 0 to WWE)
    'WWE':     (1537, 1601), # trimmed start by 3 (gap 0 from KHb-KH8)
    'ART':     (1605, 1801),
}

ACTIVE_SITES_FL = {
    'MD1': {'catalytic': [831, 923, 962],
            'pocket': [822,823,824,825,826,827,828,829,830,831,832,833,834,835,836,
                       919,920,921,922,923,924,925,926,927,961,962,966]},
    'MD2': {'catalytic': [1035, 1046, 1134, 1171],
            'pocket': [1021,1022,1023,1024,1034,1035,1036,1037,1038,1039,1040,1041,
                       1042,1043,1044,1045,1046,1047,1130,1131,1132,1133,1134,1135,
                       1136,1137,1138,1139,1140,1141,1170,1171,1175,1178]},
    'MD3': {'catalytic': [1248, 1259, 1330, 1371],
            'pocket': [1235,1236,1237,1247,1248,1249,1250,1251,1252,1253,1254,1255,
                       1256,1257,1258,1259,1260,1261,1302,1303,1304,1324,1325,1326,
                       1327,1328,1329,1330,1331,1332,1333,1334,1335,1336,1337,1369,
                       1370,1371,1375]},
    'ART': {'catalytic': [1684, 1705, 1706, 1722],
            'pocket': [1681,1682,1683,1684,1685,1688,1701,1704,1705,1706,1707,1708,
                       1709,1714,1715,1716,1721,1722,1726,1727,1781]},
}

SITE_NAMES = ['MD1', 'MD2', 'MD3', 'ART']
SITE_COLORS = {'MD1': '#e6194b', 'MD2': '#3cb44b', 'MD3': '#4363d8', 'ART': '#f58231'}

DOMAIN_UNITS = {
    'rrm1': (1, 145), 'rrm2': (146, 224), 'rrm3': (225, 314),
    'kh1-kh6': (315, 737), 'kh7a': (738, 789),
    'md1l1': (790, 1004), 'md2': (1004, 1193), 'md3': (1207, 1388),
    'khb-kh8': (1389, 1533), 'wwe': (1534, 1602), 'art': (1603, 1801),
}

CONSTRUCT_UNITS = {
    'fl':    list(DOMAIN_UNITS.keys()),
    'md':    ['md1l1', 'md2', 'md3'],
    'core':  ['kh7a', 'md1l1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art'],
    'mka':   ['md1l1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art'],
    'norrm': ['kh1-kh6', 'kh7a', 'md1l1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art'],
    'noart': ['kh1-kh6', 'kh7a', 'md1l1', 'md2', 'md3', 'khb-kh8', 'wwe'],
    'md3art':['md3', 'khb-kh8', 'wwe', 'art'],
    'fl_optimized': list(DOMAIN_UNITS.keys()),  # same units as FL (full sequence)
}

CONSTRUCT_SITES = {
    'fl':    ['MD1', 'MD2', 'MD3', 'ART'],
    'md':    ['MD1', 'MD2', 'MD3'],
    'core':  ['MD1', 'MD2', 'MD3', 'ART'],
    'mka':   ['MD1', 'MD2', 'MD3', 'ART'],
    'norrm': ['MD1', 'MD2', 'MD3', 'ART'],
    'noart': ['MD1', 'MD2', 'MD3'],
    'md3art':['MD3', 'ART'],
    'fl_optimized': ['MD1', 'MD2', 'MD3', 'ART'],
}


# ============================================================
# Auto-discover fragment simulations from fragments/ directory
# ============================================================

def register_fragment(frag_name, metadata):
    """Register a fragment as a set in SETS / CONSTRUCT_* dicts."""
    set_key = f'frag_{frag_name}'  # prefix to avoid name collisions

    SETS[set_key] = {
        'sysname': 'parp14',
        'label': frag_name,
        'color': '#888888',
    }

    units = metadata.get('units', [])
    CONSTRUCT_UNITS[set_key] = units

    # Map units -> active sites
    unit_to_site = {'md1l1': 'MD1', 'md2': 'MD2', 'md3': 'MD3', 'art': 'ART'}
    sites = [unit_to_site[u] for u in units if u in unit_to_site]
    CONSTRUCT_SITES[set_key] = sites

    # Construct domains in construct numbering
    CONSTRUCT_DOMAINS[set_key] = {
        'parp14': metadata.get('domain_ranges_construct', [])
    }


def discover_fragments():
    """Scan fragments/ directory and register all available fragments.

    Returns list of fragment set keys (with 'frag_' prefix).
    """
    frag_dir = os.path.join(CWD, 'fragments')
    if not os.path.isdir(frag_dir):
        return []
    discovered = []
    for name in sorted(os.listdir(frag_dir)):
        fdir = os.path.join(frag_dir, name)
        meta_file = os.path.join(fdir, 'metadata.json')
        if not os.path.isfile(meta_file):
            continue
        # Check at least one replicate has a trajectory
        has_traj = False
        for entry in os.listdir(fdir):
            if entry.startswith('seed-'):
                dcd = os.path.join(fdir, entry, 'parp14.dcd')
                if os.path.isfile(dcd):
                    has_traj = True
                    break
        if not has_traj:
            continue
        with open(meta_file) as f:
            meta = json.load(f)
        register_fragment(name, meta)
        discovered.append(f'frag_{name}')
    return discovered


# Directories registered via --sim-folder: set_key -> absolute folder path.
# Populated in main() before the worker pool is created so forked workers inherit it.
EXTERNAL_DIRS = {}

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


def get_sim_paths_dispatch(set_key, seed, sample):
    """get_sim_paths but handles fragment + --sim-folder sets."""
    if set_key in EXTERNAL_DIRS:
        sysname = SETS[set_key]['sysname']
        base = EXTERNAL_DIRS[set_key]
        sim_dir = os.path.join(base, f'seed-{seed}_sample-{sample}')
        if not os.path.isdir(sim_dir):
            if set_key not in _FLAT_DIRS_CACHE:
                _FLAT_DIRS_CACHE[set_key] = _flat_replicate_dirs(base, sysname)
            flat = _FLAT_DIRS_CACHE[set_key]
            idx = (seed - 1) * len(SAMPLES) + sample
            if idx < len(flat):
                sim_dir = flat[idx]
    elif set_key.startswith('frag_'):
        frag_name = set_key[5:]  # strip prefix
        sysname = SETS[set_key]['sysname']
        sim_dir = os.path.join(CWD, 'fragments', frag_name,
                                f'seed-{seed}_sample-{sample}')
    else:
        sysname = SETS[set_key]['sysname']
        sim_dir = os.path.join(CWD, set_key, f'seed-{seed}_sample-{sample}')
    dcd = os.path.join(sim_dir, f'{sysname}.dcd')
    for name in ('top.pdb', 'restart.pdb', 'checkpoint.pdb'):
        pdb = os.path.join(sim_dir, name)
        if os.path.isfile(pdb):
            return pdb, dcd
    return os.path.join(sim_dir, 'top.pdb'), dcd


def register_sim_folder(path, units=None):
    """Register an arbitrary simulation folder (from --sim-folder) as a set so
    every analysis module can process it. Domain units are read from the folder's
    metadata.json, the `units` arg, or (if the basename is a known set) that set.
    Returns the set_key under which it is registered.
    """
    import sim_registry as _reg
    folder = os.path.abspath(path)
    if not os.path.isdir(folder):
        raise FileNotFoundError(f"sim folder not found: {folder}")
    resolved_units = _reg._resolve_units(folder, units)
    meta = _reg.read_metadata(folder) or {}
    sysname = meta.get('sysname') or _reg.detect_sysname(folder)
    set_key = os.path.basename(os.path.normpath(folder))
    SETS[set_key] = {'sysname': sysname,
                     'label': meta.get('label') or set_key,
                     'color': '#444444'}
    CONSTRUCT_UNITS[set_key] = resolved_units
    CONSTRUCT_SITES[set_key] = [_reg.UNIT_TO_SITE[u] for u in resolved_units
                                if u in _reg.UNIT_TO_SITE]
    # Restraint domains (construct numbering) for the FNC module, if available.
    dranges = meta.get('domain_ranges_construct')
    if dranges is None:
        dyaml = os.path.join(folder, 'input', 'domains.yaml')
        if os.path.isfile(dyaml):
            import yaml as _yaml
            with open(dyaml) as f:
                d = _yaml.safe_load(f) or {}
            dranges = d.get(sysname) or (next(iter(d.values())) if d else None)
    if dranges is not None:
        CONSTRUCT_DOMAINS[set_key] = {sysname: dranges}
    EXTERNAL_DIRS[set_key] = folder
    return set_key

# SAA (Solid-Angle Accessibility) parameters — all lengths in nm
PROBE_RADIUS = 0.5   # SPHERE: radius (nm) of the spherical "ligand probe" swept along
                     #         each ray. A ray is blocked if any protein bead lies within
                     #         0.5 nm of it (≈ small-fragment / ADP-ribose probe size).
MAX_DIST = 5.0       # vector length cap (nm): only blocker beads within 5.0 nm along a
                     #         ray count — beyond this the active site is considered "open".
N_RAYS = 200         # number of unit VECTORS (ray directions) cast over the 4π sphere.
SEQ_SEP = 10         # sequence separation (residues, NOT nm): beads within +/-10 of a site
                     #         residue are excluded as self-blockers (own-domain scaffold).

# Global worker count (set from CLI)
N_WORKERS = 1

# Global force-recompute flag (set from CLI)
FORCE_RECOMPUTE = False


def all_data_exist(set_keys, file_pattern):
    """Return True if data files for all given sets exist.

    file_pattern: e.g. '{set}_interdomain_com_dist.npy' (uses .format()).
    A set is considered cached if at least one expected file is present.
    """
    if FORCE_RECOMPUTE:
        return False
    for sk in set_keys:
        # Substitute set name into the pattern
        candidates = [file_pattern.format(set=sk)]
        # Also try the fragment name (with frag_ prefix stripped)
        if sk.startswith('frag_'):
            candidates.append(file_pattern.format(set=sk[5:]))
        if not any(os.path.isfile(os.path.join(DATA_PATH, c))
                   for c in candidates):
            return False
    return True


def cached_data_path(set_key, pattern):
    """Return the actual cached data path for a set, trying both with and
    without frag_ prefix."""
    p1 = os.path.join(DATA_PATH, pattern.format(set=set_key))
    if os.path.isfile(p1):
        return p1
    if set_key.startswith('frag_'):
        p2 = os.path.join(DATA_PATH, pattern.format(set=set_key[5:]))
        if os.path.isfile(p2):
            return p2
    return None


# ============================================================
# FL-to-construct residue mapping
# ============================================================

def compute_fl_blocks(unit_names):
    ranges = sorted([DOMAIN_UNITS[n] for n in unit_names])
    merged = [list(ranges[0])]
    for s, e in ranges[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]


def build_fl_to_construct_map(unit_names):
    fl_blocks = compute_fl_blocks(unit_names)
    segments = []
    construct_pos = 1
    for fl_start, fl_end in fl_blocks:
        offset = construct_pos - fl_start
        segments.append((fl_start, fl_end, offset))
        construct_pos += (fl_end - fl_start + 1)

    def map_resid(fl_resid):
        for fl_s, fl_e, off in segments:
            if fl_s <= fl_resid <= fl_e:
                return fl_resid + off
        return None
    return map_resid


def get_active_sites_for_set(set_key):
    if set_key == 'fl':
        return ACTIVE_SITES_FL
    fl_to_c = build_fl_to_construct_map(CONSTRUCT_UNITS[set_key])
    sites = {}
    for sname in CONSTRUCT_SITES[set_key]:
        data = ACTIVE_SITES_FL[sname]
        sites[sname] = {
            'catalytic': [r for r in (fl_to_c(x) for x in data['catalytic']) if r is not None],
            'pocket': [r for r in (fl_to_c(x) for x in data['pocket']) if r is not None],
        }
    return sites


def get_construct_domains_for_set(set_key):
    if set_key == 'fl':
        return FL_DOMAINS
    fl_to_c = build_fl_to_construct_map(CONSTRUCT_UNITS[set_key])
    mapped = {}
    for dname, (fl_s, fl_e) in FL_DOMAINS.items():
        c_s = fl_to_c(fl_s)
        c_e = fl_to_c(fl_e)
        if c_s is not None and c_e is not None:
            mapped[dname] = (c_s, c_e)
    return mapped


# ============================================================
# Trajectory loading helpers
# ============================================================

def get_sim_paths(set_key, seed, sample):
    return get_sim_paths_dispatch(set_key, seed, sample)


def load_universe(pdb, dcd, in_memory=False):
    if not os.path.isfile(pdb) or not os.path.isfile(dcd):
        return None
    return mda.Universe(pdb, dcd, in_memory=in_memory)


def detect_n_frames(set_key):
    for seed in SEEDS:
        for sample in SAMPLES:
            pdb, dcd = get_sim_paths(set_key, seed, sample)
            u = load_universe(pdb, dcd)
            if u is not None:
                n = len(u.trajectory)
                print(f"  {set_key}: {n} frames ({n*0.01:.1f} ns)")
                return n
    return 0


def replicate_jobs(active_sets):
    """Generate all (set_key, seed, sample) tuples with valid files."""
    jobs = []
    for set_key in active_sets:
        for seed in SEEDS:
            for sample in SAMPLES:
                pdb, dcd = get_sim_paths(set_key, seed, sample)
                if os.path.isfile(pdb) and os.path.isfile(dcd):
                    jobs.append((set_key, seed, sample))
    return jobs


# ============================================================
# Helper: Debye-Huckel dielectric
# ============================================================

def dielectric_water(T):
    return 5321/T + 233.76 - 0.9297*T + 1.417e-3*T**2 - 8.292e-7*T**3


def build_force_maps(u, residues_df):
    seq = u.atoms.resnames.tolist()
    n = len(seq)
    sigmas = residues_df.loc[seq, 'sigmas'].values
    lambdas = residues_df.loc[seq, 'lambdas'].values
    charges = residues_df.loc[seq, 'q'].values

    sig = (sigmas[:, None] + sigmas[None, :]) / 2
    lam = (lambdas[:, None] + lambdas[None, :]) / 2

    RT = 8.3145e-3 * TEMP
    epsw = dielectric_water(TEMP)
    lB = 1.6021766**2 / (4 * np.pi * 8.854188e-3 * epsw * RT)
    kappa = np.sqrt(8 * np.pi * lB * IONIC * 6.022e-1)
    yukawa_prefactor = lB * RT
    qmap = yukawa_prefactor * charges[:, None] * charges[None, :]

    return sig, lam, qmap, kappa


# ============================================================
# SAA helpers (must be top-level for pickling)
# ============================================================

def fibonacci_sphere(n):
    """VECTORS: n evenly-spread unit direction vectors (|v| = 1, dimensionless)
    tiling the surface of a unit SPHERE. These are the ray directions cast
    outward from an active-site COM to probe steric openness."""
    golden = (1 + np.sqrt(5)) / 2
    indices = np.arange(n)
    theta = np.arccos(1 - 2 * (indices + 0.5) / n)
    phi = 2 * np.pi * indices / golden
    # Each row is a unit vector (x, y, z) on the unit sphere (radius 1, unitless).
    return np.column_stack([np.sin(theta)*np.cos(phi),
                            np.sin(theta)*np.sin(phi),
                            np.cos(theta)])


def compute_ray_accessibility(site_com, blocker_pos, ray_dirs, probe_radius, max_dist):
    """Fraction of the 4π solid angle around the active site that is sterically open.
    SPHERES: each protein bead is treated as a point obstacle, the ligand as a
    probe sphere of radius `probe_radius` (nm). VECTORS: `vecs` = site_com -> each
    blocker bead (nm); `ray_dirs` = unit ray directions."""
    vecs = blocker_pos - site_com            # VECTORS (nm): active-site COM -> each blocker bead
    proj = ray_dirs @ vecs.T                 # projection of each bead-vector onto each ray (nm)
    ahead = (proj > 0) & (proj < max_dist)   # bead lies in front of site, within max_dist (nm)
    perp_sq = np.sum(vecs**2, axis=1)[None, :] - proj**2  # squared perpendicular offset (nm^2)
    # A ray is blocked if a bead passes within one probe radius of it (probe sphere collides).
    blocked = np.any(ahead & (perp_sq < probe_radius**2), axis=1)
    return 1.0 - np.sum(blocked) / len(ray_dirs)  # fraction of rays unobstructed (0-1)


def compute_max_cone_angle(site_com, blocker_pos, ray_dirs, probe_radius, max_dist):
    """Half-angle (degrees) of the widest unobstructed approach cone at the site.
    Same probe SPHERE (radius `probe_radius` nm) and bead-VECTORS (nm) as above;
    `mean_dir` is a unit VECTOR pointing down the center of the open cone."""
    vecs = blocker_pos - site_com            # VECTORS (nm): site COM -> each blocker bead
    proj = ray_dirs @ vecs.T                 # bead-vector projected onto each ray (nm)
    ahead = (proj > 0) & (proj < max_dist)   # in front of site and within max_dist (nm)
    perp_sq = np.sum(vecs**2, axis=1)[None, :] - proj**2  # squared perp offset (nm^2)
    blocked = np.any(ahead & (perp_sq < probe_radius**2), axis=1)  # probe sphere collides
    free_dirs = ray_dirs[~blocked]           # unit VECTORS of the unblocked rays
    if len(free_dirs) == 0:
        return 0.0
    mean_dir = np.mean(free_dirs, axis=0)    # VECTOR: average open direction (cone axis)
    norm = np.linalg.norm(mean_dir)
    if norm < 1e-10:
        return 180.0
    mean_dir /= norm                         # normalize to unit VECTOR (cone axis, |v|=1)
    # widest angular deviation of any open ray from the cone axis -> cone half-angle (deg)
    return np.degrees(np.arccos(np.clip(np.min(free_dirs @ mean_dir), -1, 1)))


def compute_shell_density(site_com, all_pos, r_inner=2.0, r_outer=5.0):
    """Protein volume fraction packed in a spherical shell around the active site.
    SPHERES: an inner sphere of radius r_inner = 2.0 nm and an outer sphere of
    radius r_outer = 5.0 nm define the shell; each amino-acid bead is a sphere of
    radius 0.25 nm (diameter 0.5 nm). VECTOR magnitudes `dists` = |site_com -> bead| (nm)."""
    dists = np.linalg.norm(all_pos - site_com, axis=1)        # |VECTOR| site COM -> each bead (nm)
    n_in_shell = np.sum((dists >= r_inner) & (dists < r_outer))  # beads inside the 2-5 nm shell
    bead_vol = (4/3) * np.pi * (0.25)**3                      # one bead SPHERE volume, r=0.25 nm (nm^3)
    shell_vol = (4/3) * np.pi * (r_outer**3 - r_inner**3)     # shell volume between the two spheres (nm^3)
    return n_in_shell * bead_vol / shell_vol                  # occupied volume fraction (unitless)


# ============================================================
# WORKER FUNCTIONS (top-level, picklable for ProcessPoolExecutor)
# ============================================================

def _worker_conf_prop(set_key, seed, sample):
    """Compute Rg and Ree for one replicate."""
    pdb, dcd = get_sim_paths(set_key, seed, sample)
    u = load_universe(pdb, dcd)
    if u is None:
        return None
    residues = pd.read_csv(RESIDUES_FILE).set_index('three')
    ag = u.select_atoms('all')
    seq = ag.resnames.tolist()
    rgs = calc_rg(u, ag, seq, residues, start=SKIP_FRAMES)
    rees, _, _ = calc_ete(u, ag, start=SKIP_FRAMES)
    return {'set_key': set_key, 'seed': seed, 'sample': sample,
            'rg': rgs, 'ree': rees}


def _worker_dmap(set_key, seed, sample):
    """Compute interdomain COM distances for one replicate."""
    pdb, dcd = get_sim_paths(set_key, seed, sample)
    u = load_universe(pdb, dcd)
    if u is None:
        return None
    domains = get_construct_domains_for_set(set_key)
    dnames = list(domains.keys())
    n_dom = len(dnames)
    dom_ags = [u.select_atoms(f'resid {s}:{e}') for s, e in domains.values()]

    # DISTANCE MEASURED: center-of-mass (COM) to COM distance between every pair of
    # structured domains, per frame. Each domain is reduced to a single point = its COM.
    n_eq = len(u.trajectory) - SKIP_FRAMES
    frame_dists = np.zeros((n_eq, n_dom, n_dom))
    for t, ts in enumerate(u.trajectory[SKIP_FRAMES:]):
        coms = np.array([ag.center_of_mass() / 10.0 for ag in dom_ags])  # domain COMs (nm; /10 = Å->nm)
        for i in range(n_dom):
            for j in range(i, n_dom):
                d = np.linalg.norm(coms[i] - coms[j])  # |VECTOR| COM_i -> COM_j = inter-domain distance (nm)
                frame_dists[t, i, j] = d
                frame_dists[t, j, i] = d

    return {'set_key': set_key, 'dnames': dnames,
            'mean': np.mean(frame_dists, axis=0),
            'std': np.std(frame_dists, axis=0)}


def _worker_cmap(set_key, seed, sample):
    """Compute interdomain contact frequency for one replicate."""
    pdb, dcd = get_sim_paths(set_key, seed, sample)
    u = load_universe(pdb, dcd)
    if u is None:
        return None
    domains = get_construct_domains_for_set(set_key)
    dnames = list(domains.keys())
    n_dom = len(dnames)
    dom_ags = {d: u.select_atoms(f'resid {s}:{e}') for d, (s, e) in domains.items()}

    # DISTANCE MEASURED: residue–residue bead distances between two domains. A contact
    # is counted when two beads fall within a contact SPHERE of radius cutoff = 1.0 nm.
    freq = np.zeros((n_dom, n_dom))
    for i, di in enumerate(dnames):
        for j, dj in enumerate(dnames):
            if j < i:
                continue
            cmap = cmap_traj(u, dom_ags[di], dom_ags[dj], cutoff=1.0, start=SKIP_FRAMES)  # 1.0 nm contact cutoff
            freq[i, j] = np.mean(cmap)
            if j > i:
                freq[j, i] = freq[i, j]

    return {'set_key': set_key, 'dnames': dnames, 'freq': freq}


def _worker_fnc(set_key, seed, sample):
    """Compute FNC per restrained domain for one replicate."""
    pdb, dcd = get_sim_paths(set_key, seed, sample)
    if not os.path.isfile(pdb) or not os.path.isfile(dcd):
        return None
    sysname = SETS[set_key]['sysname']
    domain_ranges = CONSTRUCT_DOMAINS[set_key][sysname]

    # DISTANCE MEASURED: intra-domain bead–bead distances vs the reference (PDB) structure.
    # A "native contact" is any bead pair within a contact SPHERE of radius cutoff = 1.0 nm
    # in the reference; FNC = fraction of those still satisfied each frame (rigidity check).
    fnc_results = {}
    for ds, de in domain_ranges:
        sel = f'resid {ds}:{de}'
        u = mda.Universe(pdb, dcd, in_memory=True)
        uref = mda.Universe(pdb)
        fnc = calc_fnc(u, uref, sel, cutoff=1.0)  # 1.0 nm native-contact cutoff sphere
        fnc_eq = fnc[SKIP_FRAMES:]
        fnc_results[f'{ds}_{de}'] = np.mean(fnc_eq)

    return {'set_key': set_key, 'fnc': fnc_results}


def _worker_wcn(set_key, seed, sample):
    """Compute WCN at active sites for one replicate."""
    pdb, dcd = get_sim_paths(set_key, seed, sample)
    u = load_universe(pdb, dcd)
    if u is None:
        return None
    sites = get_active_sites_for_set(set_key)
    available_sites = CONSTRUCT_SITES[set_key]
    ag = u.select_atoms('all')

    n_total = len(u.trajectory)
    sample_frames = np.linspace(SKIP_FRAMES, n_total - 1,
                                min(20, n_total - SKIP_FRAMES), dtype=int)

    # DISTANCE MEASURED: number of neighboring beads around each catalytic residue,
    # smoothly weighted by a coordination SPHERE of characteristic radius r0 = 0.7 nm
    # (higher WCN = more buried/crowded active site).
    site_wcn = {s: [] for s in available_sites}
    for frame_idx in sample_frames:
        u.trajectory[frame_idx]
        pos = ag.positions / 10.0                 # all bead positions (nm; /10 = Å->nm)
        wcn_all = calc_wcn(None, pos, fdomains=None, ssonly=False, r0=0.7)  # 0.7 nm coordination radius
        for sname in available_sites:
            indices = [r - 1 for r in sites[sname]['catalytic'] if r - 1 < len(wcn_all)]
            if indices:
                site_wcn[sname].append(np.mean(wcn_all[indices]))

    return {'set_key': set_key,
            'wcn': {s: np.mean(v) if v else np.nan for s, v in site_wcn.items()}}


def _worker_active_sites(set_key, seed, sample):
    """Compute active site radial position, RMSF, inter-site distances."""
    pdb, dcd = get_sim_paths(set_key, seed, sample)
    u = load_universe(pdb, dcd)
    if u is None:
        return None
    sites = get_active_sites_for_set(set_key)
    available_sites = CONSTRUCT_SITES[set_key]
    residues = pd.read_csv(RESIDUES_FILE).set_index('three')
    ag = u.select_atoms('all')
    seq = ag.resnames.tolist()

    rgs = calc_rg(u, ag, seq, residues, start=SKIP_FRAMES)
    n_eq = len(u.trajectory) - SKIP_FRAMES

    radial = {}
    rmsf = {}
    site_coms = {}

    for sname in available_sites:
        cat_resids = sites[sname]['catalytic']
        sel = ' or '.join([f'resid {r}' for r in cat_resids])
        site_ag = u.select_atoms(sel)
        if len(site_ag) == 0:
            continue
        coms = np.zeros((n_eq, 3))
        prot_com = np.zeros((n_eq, 3))
        for t, ts in enumerate(u.trajectory[SKIP_FRAMES:]):
            coms[t] = site_ag.center_of_mass() / 10.0
            prot_com[t] = ag.center_of_mass() / 10.0
        site_coms[sname] = coms
        dist_to_com = np.linalg.norm(coms - prot_com, axis=1)
        radial[sname] = np.mean(dist_to_com / rgs)
        mean_pos = coms.mean(axis=0)
        rmsf[sname] = np.sqrt(np.mean(np.sum((coms - mean_pos)**2, axis=1)))

    inter_dists = {}
    avail = [s for s in available_sites if s in site_coms]
    for i in range(len(avail)):
        for j in range(i+1, len(avail)):
            pair = f'{avail[i]}-{avail[j]}'
            d = np.linalg.norm(site_coms[avail[i]] - site_coms[avail[j]], axis=1)
            inter_dists[pair] = np.mean(d)

    return {'set_key': set_key, 'radial': radial, 'rmsf': rmsf,
            'inter_dists': inter_dists}


def _worker_accessibility(set_key, seed, sample):
    """Compute SAA, cone angle, shell density at active sites."""
    pdb, dcd = get_sim_paths(set_key, seed, sample)
    u = load_universe(pdb, dcd)
    if u is None:
        return None
    sites = get_active_sites_for_set(set_key)
    available_sites = CONSTRUCT_SITES[set_key]
    ag = u.select_atoms('all')
    all_resids = ag.resids
    ray_dirs = fibonacci_sphere(N_RAYS)

    # Precompute blocker masks
    blocker_masks = {}
    for sname in available_sites:
        site_resids = set(sites[sname]['catalytic']) | set(sites[sname]['pocket'])
        mask = np.ones(len(ag), dtype=bool)
        for idx, resid in enumerate(all_resids):
            if resid in site_resids or any(abs(resid - sr) <= SEQ_SEP for sr in site_resids):
                mask[idx] = False
        blocker_masks[sname] = mask

    n_total = len(u.trajectory)
    sample_frames = np.linspace(SKIP_FRAMES, n_total - 1,
                                min(20, n_total - SKIP_FRAMES), dtype=int)

    saa_vals = {s: [] for s in available_sites}
    cone_vals = {s: [] for s in available_sites}
    shell_vals = {s: [] for s in available_sites}

    for frame_idx in sample_frames:
        u.trajectory[frame_idx]
        all_pos = ag.positions / 10.0
        for sname in available_sites:
            cat_indices = [r - 1 for r in sites[sname]['catalytic'] if r - 1 < len(all_pos)]
            if not cat_indices:
                continue
            site_com = np.mean(all_pos[cat_indices], axis=0)
            blocker_pos = all_pos[blocker_masks[sname]]
            saa_vals[sname].append(
                compute_ray_accessibility(site_com, blocker_pos, ray_dirs, PROBE_RADIUS, MAX_DIST))
            cone_vals[sname].append(
                compute_max_cone_angle(site_com, blocker_pos, ray_dirs, PROBE_RADIUS, MAX_DIST))
            shell_vals[sname].append(
                compute_shell_density(site_com, all_pos))

    return {'set_key': set_key,
            'saa': {s: np.mean(v) if v else np.nan for s, v in saa_vals.items()},
            'cone': {s: np.mean(v) if v else np.nan for s, v in cone_vals.items()},
            'shell': {s: np.mean(v) if v else np.nan for s, v in shell_vals.items()}}


# ============================================================
# Parallel dispatch helper
# ============================================================

def parallel_map(worker_fn, jobs, label):
    """Submit jobs to ProcessPoolExecutor, return list of results."""
    n_jobs = len(jobs)
    print(f"  Dispatching {n_jobs} jobs across {N_WORKERS} workers...")
    results = []

    if N_WORKERS == 1:
        # Sequential fallback
        for i, (sk, se, sa) in enumerate(jobs):
            r = worker_fn(sk, se, sa)
            if r is not None:
                results.append(r)
            if (i + 1) % 10 == 0 or i + 1 == n_jobs:
                print(f"    {i+1}/{n_jobs} done")
    else:
        with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
            futures = {pool.submit(worker_fn, sk, se, sa): (sk, se, sa)
                       for sk, se, sa in jobs}
            done_count = 0
            for future in as_completed(futures):
                done_count += 1
                try:
                    r = future.result()
                    if r is not None:
                        results.append(r)
                except Exception as e:
                    job = futures[future]
                    print(f"    FAILED: {job} — {e}")
                if done_count % 10 == 0 or done_count == n_jobs:
                    print(f"    {done_count}/{n_jobs} done")

    print(f"  Collected {len(results)} results")
    return results


def group_by_set(results):
    """Group list of result dicts by set_key."""
    grouped = defaultdict(list)
    for r in results:
        grouped[r['set_key']].append(r)
    return dict(grouped)


def merge_savez(path, new_dict, active_sets):
    """Save a shared .npz cache WITHOUT clobbering sets that weren't recomputed.

    The conf-prop / active-site / accessibility / wcn caches are single files
    keyed as '<set>_...'. Running on a subset (e.g. one --sim-folder) used to
    overwrite the whole file. This loads the existing cache, drops only the keys
    belonging to the sets being (re)written, then merges in the new data.
    """
    merged = {}
    if os.path.isfile(path):
        try:
            old = np.load(path, allow_pickle=True)
            prefixes = tuple(f'{s}_' for s in active_sets)
            for k in old.files:
                if not k.startswith(prefixes):
                    merged[k] = old[k]
        except Exception as e:
            print(f"  WARN: could not read existing {os.path.basename(path)} "
                  f"({e}); rewriting")
    merged.update(new_dict)
    np.savez(path, **merged)


# ============================================================
# 1. Conformational Properties
# ============================================================

def run_conf_prop(active_sets):
    print("\n" + "=" * 70)
    print("CONFORMATIONAL PROPERTIES (Rg, Ree) — parallel")
    print("=" * 70)

    cache_file = os.path.join(DATA_PATH, 'conf_prop.npz')
    # Try to reload from cache (skip slow parallel compute)
    if os.path.isfile(cache_file) and not FORCE_RECOMPUTE:
        try:
            cached = np.load(cache_file, allow_pickle=True)
            all_data = {}
            for sk in active_sets:
                rg_keys = sorted([k for k in cached.files
                                   if k.startswith(f'{sk}_rg_rep')])
                ree_keys = sorted([k for k in cached.files
                                    if k.startswith(f'{sk}_ree_rep')])
                if rg_keys and ree_keys:
                    all_data[sk] = {
                        'rg': [cached[k] for k in rg_keys],
                        'ree': [cached[k] for k in ree_keys],
                    }
            if all_data:
                print(f"  Loaded cached data ({len(all_data)} sets, "
                      f"--force to recompute)")
                _plot_conf_prop(active_sets, all_data)
                return all_data
        except Exception as e:
            print(f"  Cache reload failed ({e}) — recomputing")

    jobs = replicate_jobs(active_sets)
    results = parallel_map(_worker_conf_prop, jobs, "conf_prop")
    by_set = group_by_set(results)

    all_data = {}
    for set_key in active_sets:
        info = SETS[set_key]
        reps = by_set.get(set_key, [])
        rg_list = [r['rg'] for r in reps]
        ree_list = [r['ree'] for r in reps]
        all_data[set_key] = {'rg': rg_list, 'ree': ree_list}
        if rg_list:
            eq_rgs = np.array([np.mean(r) for r in rg_list])
            eq_rees = np.array([np.mean(r) for r in ree_list])
            print(f"  {info['label']}: {len(rg_list)} reps, "
                  f"Rg={np.mean(eq_rgs):.3f}+/-{np.std(eq_rgs):.3f}, "
                  f"Ree={np.mean(eq_rees):.3f}+/-{np.std(eq_rees):.3f} nm")

    # Save
    save_dict = {}
    for set_key, d in all_data.items():
        for i, rg in enumerate(d['rg']):
            save_dict[f'{set_key}_rg_rep{i}'] = rg
        for i, ree in enumerate(d['ree']):
            save_dict[f'{set_key}_ree_rep{i}'] = ree
    merge_savez(os.path.join(DATA_PATH, 'conf_prop.npz'), save_dict, active_sets)
    print(f"  Saved: data/conf_prop.npz")

    _plot_conf_prop(active_sets, all_data)
    return all_data


def _plot_conf_prop(active_sets, all_data):
    n_sets = len(active_sets)

    # Rg time series
    fig, axes = plt.subplots(1, n_sets, figsize=(5*n_sets, 5), squeeze=False)
    for idx, set_key in enumerate(active_sets):
        ax = axes[0, idx]
        info = SETS[set_key]
        rg_list = all_data[set_key]['rg']
        if not rg_list:
            ax.set_title(info['label']); continue
        # Truncate to shortest replicate so mean/std are well-defined
        min_frames = min(len(rg) for rg in rg_list)
        time_ns = np.arange(SKIP_FRAMES, SKIP_FRAMES + min_frames) * 0.01
        for rg in rg_list:
            t = np.arange(SKIP_FRAMES, SKIP_FRAMES + len(rg)) * 0.01
            ax.plot(t, rg, lw=0.3, alpha=0.3, color=info['color'])
        rg_arr = np.array([rg[:min_frames] for rg in rg_list])
        mean_rg = np.mean(rg_arr, axis=0)
        std_rg = np.std(rg_arr, axis=0)
        ax.plot(time_ns, mean_rg, lw=2, color='black', label='Mean')
        ax.fill_between(time_ns, mean_rg-std_rg, mean_rg+std_rg, alpha=0.2, color='black')
        ax.axhline(np.mean(mean_rg), color='red', ls='--', lw=1, label=f'Eq={np.mean(mean_rg):.2f}')
        ax.set_xlabel('Time (ns)'); ax.set_ylabel('Rg (nm)')
        ax.set_title(info['label']); ax.legend(fontsize=8)
    fig.suptitle('Radius of Gyration (mass-weighted)', fontsize=14, y=1.02)
    fig.tight_layout()
    for ext in ['png', 'svg']:
        fig.savefig(os.path.join(FIG_PATH, f'rg_timeseries.{ext}'), dpi=150, bbox_inches='tight')
    plt.close()

    # Distributions
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for pi, (prop, label) in enumerate([('rg', 'Rg (nm)'), ('ree', 'Ree (nm)')]):
        ax = axes[pi]
        positions, labels, colors, all_vals = [], [], [], []
        for idx, set_key in enumerate(active_sets):
            vals = all_data[set_key][prop]
            if vals:
                all_vals.append([np.mean(v) for v in vals])
                positions.append(idx); labels.append(SETS[set_key]['label'])
                colors.append(SETS[set_key]['color'])
        if all_vals:
            vp = ax.violinplot(all_vals, positions=positions, showmeans=True, showextrema=True)
            for i, body in enumerate(vp['bodies']):
                body.set_facecolor(colors[i]); body.set_alpha(0.6)
            for i, (pos, eq) in enumerate(zip(positions, all_vals)):
                jitter = np.random.default_rng(42).normal(0, 0.05, len(eq))
                ax.scatter(pos+jitter, eq, c=colors[i], s=20, alpha=0.7, edgecolors='black', linewidths=0.5)
            ax.set_xticks(positions); ax.set_xticklabels(labels, fontsize=9)
            ax.set_ylabel(label); ax.set_title(f'{prop.upper()} Distribution')
    fig.tight_layout()
    for ext in ['png', 'svg']:
        fig.savefig(os.path.join(FIG_PATH, f'conf_prop_distributions.{ext}'), dpi=150, bbox_inches='tight')
    plt.close()

    # Convergence
    fig, axes = plt.subplots(1, n_sets, figsize=(5*n_sets, 4), squeeze=False)
    for idx, set_key in enumerate(active_sets):
        ax = axes[0, idx]; info = SETS[set_key]
        rg_list = all_data[set_key]['rg']
        if not rg_list:
            ax.set_title(info['label']); continue
        min_frames = min(len(rg) for rg in rg_list)
        rg_arr = np.array([rg[:min_frames] for rg in rg_list])
        n_rep, n_frames = rg_arr.shape
        cum_mean = np.cumsum(rg_arr, axis=1) / np.arange(1, n_frames+1)
        time_ns = np.arange(n_frames)*0.01 + SKIP_FRAMES*0.01
        for rep in range(n_rep):
            ax.plot(time_ns, cum_mean[rep], lw=0.3, alpha=0.4, color=info['color'])
        ax.plot(time_ns, np.mean(cum_mean, axis=0), lw=2, color='black', label='Ensemble')
        ax.set_xlabel('Time (ns)'); ax.set_ylabel('Cumulative Rg (nm)')
        ax.set_title(info['label']); ax.legend(fontsize=8)
    fig.suptitle('Rg Convergence', fontsize=14, y=1.02); fig.tight_layout()
    for ext in ['png', 'svg']:
        fig.savefig(os.path.join(FIG_PATH, f'rg_convergence.{ext}'), dpi=150, bbox_inches='tight')
    plt.close()

    # Secondary version: block-averaging standard error (Flyvbjerg & Petersen
    # 1989). The panel above is a cumulative running mean, which always looks
    # "converged" at long times regardless of autocorrelation -- it doesn't
    # actually test anything. Block SE should plateau once block size exceeds
    # the Rg autocorrelation time; a plateau is the real convergence signal,
    # and its height is the correct (correlation-corrected) standard error.
    fig, axes = plt.subplots(1, n_sets, figsize=(5*n_sets, 4), squeeze=False)
    for idx, set_key in enumerate(active_sets):
        ax = axes[0, idx]; info = SETS[set_key]
        rg_list = all_data[set_key]['rg']
        if not rg_list:
            ax.set_title(info['label']); continue
        min_frames = min(len(rg) for rg in rg_list)
        max_pow = int(np.floor(np.log2(max(min_frames // 2, 1))))
        block_sizes = 2 ** np.arange(0, max_pow + 1)
        block_ns = block_sizes * 0.01
        curves = []
        for rg in rg_list:
            x = np.asarray(rg[:min_frames])
            se = np.full(len(block_sizes), np.nan)
            for bi, b in enumerate(block_sizes):
                nb = min_frames // b
                if nb < 2:
                    continue
                block_means = x[:nb*b].reshape(nb, b).mean(axis=1)
                se[bi] = np.std(block_means, ddof=1) / np.sqrt(nb)
            curves.append(se)
        curves = np.array(curves)
        for c in curves:
            ax.plot(block_ns, c, lw=0.5, alpha=0.3, color=info['color'])
        ax.plot(block_ns, np.nanmean(curves, axis=0), lw=2, color='black', label='Mean')
        ax.set_xscale('log')
        ax.set_xlabel('Block size (ns)'); ax.set_ylabel('Block SE of Rg (nm)')
        ax.set_title(info['label']); ax.legend(fontsize=8)
    fig.suptitle('Rg Convergence — Block-Averaging Standard Error (Flyvbjerg & Petersen)',
                 fontsize=13, y=1.02)
    fig.tight_layout()
    for ext in ['png', 'svg']:
        fig.savefig(os.path.join(FIG_PATH, f'rg_convergence_block.{ext}'), dpi=150, bbox_inches='tight')
    plt.close()

    print(f"  Saved: figures/rg_timeseries, conf_prop_distributions, rg_convergence(_block)")


# ============================================================
# 2. Distance Maps
# ============================================================

def run_dmap_analysis(active_sets):
    print("\n" + "=" * 70)
    print("DISTANCE MAP ANALYSIS — parallel")
    print("=" * 70)

    # Cache-aware: separate sets into cached vs needs-compute
    needs_compute = [
        sk for sk in active_sets
        if FORCE_RECOMPUTE or
        cached_data_path(sk, '{set}_interdomain_com_dist.npy') is None
    ]
    cached = [sk for sk in active_sets if sk not in needs_compute]
    if cached:
        print(f"  Loading cached data for {len(cached)} sets: "
              f"{', '.join(cached)}")
    if not needs_compute:
        # All cached — just rebuild figures from saved data
        by_set = {}
        for sk in active_sets:
            p_mean = cached_data_path(sk, '{set}_interdomain_com_dist.npy')
            p_std = cached_data_path(sk, '{set}_interdomain_com_dist_std.npy')
            if p_mean and p_std:
                # Reconstruct minimal reps list (one fake "rep" with the
                # already-aggregated mean/std)
                mean_d = np.load(p_mean)
                std_d = np.load(p_std)
                # Try to recover dnames from CONSTRUCT_DOMAINS or fragment metadata
                if sk.startswith('frag_'):
                    units = CONSTRUCT_UNITS.get(sk, [])
                    dnames = [u.upper() for u in units]
                else:
                    domains = get_construct_domains_for_set(sk)
                    dnames = list(domains.keys())
                # Pad dnames to match matrix size if mismatch
                n = mean_d.shape[0]
                if len(dnames) != n:
                    dnames = [f'D{i}' for i in range(n)]
                by_set[sk] = [{'mean': mean_d, 'std': std_d,
                               'dnames': dnames}]
    else:
        print(f"  Computing for {len(needs_compute)} sets: "
              f"{', '.join(needs_compute)}")
        jobs = replicate_jobs(needs_compute)
        results = parallel_map(_worker_dmap, jobs, "dmap")
        by_set = group_by_set(results)
        # Also load any cached
        for sk in cached:
            p_mean = cached_data_path(sk, '{set}_interdomain_com_dist.npy')
            p_std = cached_data_path(sk, '{set}_interdomain_com_dist_std.npy')
            if p_mean and p_std:
                mean_d = np.load(p_mean)
                std_d = np.load(p_std)
                if sk.startswith('frag_'):
                    units = CONSTRUCT_UNITS.get(sk, [])
                    dnames = [u.upper() for u in units]
                else:
                    domains = get_construct_domains_for_set(sk)
                    dnames = list(domains.keys())
                if len(dnames) != mean_d.shape[0]:
                    dnames = [f'D{i}' for i in range(mean_d.shape[0])]
                by_set[sk] = [{'mean': mean_d, 'std': std_d,
                               'dnames': dnames}]

    for set_key in active_sets:
        info = SETS[set_key]
        reps = by_set.get(set_key, [])
        if not reps:
            continue
        dnames = reps[0]['dnames']
        n_dom = len(dnames)
        # If only one "rep" entry from cache, use it as-is; else average
        if len(reps) == 1 and 'mean' in reps[0]:
            mean_d = reps[0]['mean']
            std_d = reps[0]['std']
        else:
            mean_d = np.mean([r['mean'] for r in reps], axis=0)
            std_d = np.mean([r['std'] for r in reps], axis=0)
            np.save(os.path.join(DATA_PATH,
                                  f'{set_key}_interdomain_com_dist.npy'),
                    mean_d)
            np.save(os.path.join(DATA_PATH,
                                  f'{set_key}_interdomain_com_dist_std.npy'),
                    std_d)

        print(f"  {info['label']}: {len(reps)} reps, {n_dom} domains")

        fig, ax = plt.subplots(figsize=(8, 7))
        im = ax.imshow(mean_d, cmap='viridis', origin='lower')
        ax.set_xticks(range(n_dom)); ax.set_xticklabels(dnames, rotation=45, ha='right', fontsize=8)
        ax.set_yticks(range(n_dom)); ax.set_yticklabels(dnames, fontsize=8)
        plt.colorbar(im, ax=ax, label='COM distance (nm)', shrink=0.8)
        for i in range(n_dom):
            for j in range(n_dom):
                ax.text(j, i, f'{mean_d[i,j]:.1f}', ha='center', va='center', fontsize=6,
                        color='white' if mean_d[i,j] > np.median(mean_d) else 'black')
        ax.set_title(f'{info["label"]} — Interdomain COM Distance')
        fig.tight_layout()
        for ext in ['png', 'svg']:
            fig.savefig(os.path.join(FIG_PATH, f'{set_key}_interdomain_dmap.{ext}'), dpi=150, bbox_inches='tight')
        plt.close()
    print(f"  Saved: figures/*_interdomain_dmap")


# ============================================================
# 3. Contact Maps
# ============================================================

def run_cmap_analysis(active_sets):
    print("\n" + "=" * 70)
    print("CONTACT MAP ANALYSIS — parallel")
    print("=" * 70)

    # Cache-aware: split into cached vs needs-compute
    needs_compute = [
        sk for sk in active_sets
        if FORCE_RECOMPUTE or
        cached_data_path(sk, '{set}_interdomain_contact_freq.npy') is None
    ]
    cached_sets = [sk for sk in active_sets if sk not in needs_compute]
    if cached_sets:
        print(f"  Cached: {', '.join(cached_sets)}")
    if not needs_compute:
        print(f"  All data cached (--force to recompute). Rebuilding figures.")
        for sk in active_sets:
            p = cached_data_path(sk, '{set}_interdomain_contact_freq.npy')
            if not p:
                continue
            freq = np.load(p)
            info = SETS[sk]
            n_dom = freq.shape[0]
            if sk.startswith('frag_'):
                units = CONSTRUCT_UNITS.get(sk, [])
                dnames = [u.upper() for u in units][:n_dom]
            else:
                domains = get_construct_domains_for_set(sk)
                dnames = list(domains.keys())[:n_dom]
            if len(dnames) < n_dom:
                dnames = [f'D{i}' for i in range(n_dom)]
            _plot_cmap(sk, info, freq, dnames)
        return

    print(f"  Computing for: {', '.join(needs_compute)}")
    # Use seed=1 only (cmap_traj is expensive)
    jobs = [(sk, 1, sa) for sk in needs_compute for sa in SAMPLES
            if os.path.isfile(get_sim_paths(sk, 1, sa)[0]) and
               os.path.isfile(get_sim_paths(sk, 1, sa)[1])]

    results = parallel_map(_worker_cmap, jobs, "cmap")
    by_set = group_by_set(results)

    for set_key in active_sets:
        info = SETS[set_key]
        reps = by_set.get(set_key, [])
        if not reps:
            continue
        dnames = reps[0]['dnames']
        n_dom = len(dnames)
        freq = np.mean([r['freq'] for r in reps], axis=0)

        np.save(os.path.join(DATA_PATH, f'{set_key}_interdomain_contact_freq.npy'), freq)
        print(f"  {info['label']}: {len(reps)} reps, {n_dom} domains")
        _plot_cmap(set_key, info, freq, dnames)
    print(f"  Saved: figures/*_interdomain_cmap")


def _plot_cmap(set_key, info, freq, dnames):
    """Plot one cmap heatmap from frequency matrix + domain names."""
    n_dom = freq.shape[0]
    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(freq, cmap='hot_r', origin='lower', vmin=0)
    ax.set_xticks(range(n_dom))
    ax.set_xticklabels(dnames, rotation=45, ha='right', fontsize=8)
    ax.set_yticks(range(n_dom))
    ax.set_yticklabels(dnames, fontsize=8)
    plt.colorbar(im, ax=ax, label='Mean contact frequency', shrink=0.8)
    for i in range(n_dom):
        for j in range(n_dom):
            ax.text(j, i, f'{freq[i,j]:.3f}', ha='center', va='center',
                    fontsize=6,
                    color='white' if freq[i, j] > 0.3 * np.max(freq) else 'black')
    ax.set_title(f'{info["label"]} — Interdomain Contact Frequency')
    fig.tight_layout()
    for ext in ['png', 'svg']:
        fig.savefig(os.path.join(FIG_PATH,
                                  f'{set_key}_interdomain_cmap.{ext}'),
                    dpi=150, bbox_inches='tight')
    plt.close()


# ============================================================
# 4. FNC
# ============================================================

def run_fnc_analysis(active_sets):
    print("\n" + "=" * 70)
    print("FRACTION OF NATIVE CONTACTS — parallel")
    print("=" * 70)

    # Cache-aware
    needs_compute = [
        sk for sk in active_sets
        if FORCE_RECOMPUTE or
        cached_data_path(sk, '{set}_fnc.npz') is None
    ]
    if not needs_compute:
        print(f"  All FNC data cached for active sets (--force to "
              f"recompute). Skipping.")
        return None
    if len(needs_compute) < len(active_sets):
        print(f"  Cached: {set(active_sets) - set(needs_compute)}, "
              f"computing: {needs_compute}")
    active_sets = needs_compute

    jobs = replicate_jobs(active_sets)
    results = parallel_map(_worker_fnc, jobs, "fnc")
    by_set = group_by_set(results)

    for set_key in active_sets:
        info = SETS[set_key]
        reps = by_set.get(set_key, [])
        if not reps:
            continue
        sysname = SETS[set_key]['sysname']
        domain_ranges = CONSTRUCT_DOMAINS[set_key][sysname]

        fnc_means = {}
        for ds, de in domain_ranges:
            key = f'{ds}_{de}'
            vals = [r['fnc'][key] for r in reps if key in r['fnc']]
            if vals:
                fnc_means[f'dom_{key}'] = np.array(vals)
                print(f"  {info['label']} {ds}-{de}: FNC={np.mean(vals):.3f}+/-{np.std(vals):.3f}")

        np.savez(os.path.join(DATA_PATH, f'{set_key}_fnc.npz'), **fnc_means)

    # Plot (violin + points per restrained domain)
    rng = np.random.default_rng(42)
    fig, axes = plt.subplots(1, len(active_sets), figsize=(5*len(active_sets), 5), squeeze=False)
    for idx, set_key in enumerate(active_sets):
        ax = axes[0, idx]; info = SETS[set_key]; color = info['color']
        fnc_file = os.path.join(DATA_PATH, f'{set_key}_fnc.npz')
        if not os.path.isfile(fnc_file):
            ax.set_title(info['label']); continue
        data = np.load(fnc_file)
        dnames = list(data.keys())
        short = [k.replace('dom_', '').replace('_', '-') for k in dnames]
        for di, k in enumerate(dnames):
            vals = data[k]
            if len(vals) >= 2:
                parts = ax.violinplot([vals], positions=[di], widths=0.6,
                                      showmeans=False, showextrema=False)
                for body in parts['bodies']:
                    body.set_facecolor(color); body.set_alpha(0.35); body.set_edgecolor(color)
                m, s = np.mean(vals), np.std(vals)
                ax.plot([di-0.15, di+0.15], [m, m], color=color, lw=2, zorder=4)
                ax.plot([di, di], [m-s, m+s], color=color, lw=1.5, zorder=4)
                ax.plot([di-0.08, di+0.08], [m-s, m-s], color=color, lw=1.5, zorder=4)
                ax.plot([di-0.08, di+0.08], [m+s, m+s], color=color, lw=1.5, zorder=4)
            jitter = rng.normal(0, 0.05, len(vals))
            ax.scatter(di+jitter, vals, c=color, s=18, alpha=0.7, zorder=5,
                       edgecolors='black', linewidths=0.4)
        ax.set_xticks(range(len(dnames))); ax.set_xticklabels(short, fontsize=7, rotation=45, ha='right')
        ax.set_ylabel('FNC'); ax.set_ylim(0, 1.1); ax.set_title(info['label'])
    fig.suptitle('Fraction of Native Contacts per Restrained Domain', fontsize=13, y=1.02)
    fig.tight_layout()
    for ext in ['png', 'svg']:
        fig.savefig(os.path.join(FIG_PATH, f'fnc_per_domain.{ext}'), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: figures/fnc_per_domain")


# ============================================================
# 5. Energy Decomposition (not parallelized — 1 rep per set)
# ============================================================

# Canonical domain order (all 11 FL domains) used for all energy plots
ALL_DOMAIN_NAMES = list(FL_DOMAINS.keys())  # RRM1..ART
N_ALL_DOMAINS = len(ALL_DOMAIN_NAMES)

def _compute_energy_matrix(u, domains, sig, lam, qmap, kappa):
    """Compute interdomain AH and Yukawa energy matrices."""
    dnames = list(domains.keys())
    n_dom = len(dnames)
    dom_ags = {d: u.select_atoms(f'resid {s}:{e}') for d, (s, e) in domains.items()}

    n_total = len(u.trajectory)
    sample_frames = np.linspace(SKIP_FRAMES, n_total - 1,
                                min(10, n_total - SKIP_FRAMES), dtype=int)

    ah_matrix = np.zeros((n_dom, n_dom))
    yu_matrix = np.zeros((n_dom, n_dom))

    for frame_idx in sample_frames:
        u.trajectory[frame_idx]
        for i, di in enumerate(dnames):
            for j, dj in enumerate(dnames):
                if j <= i:
                    continue
                idx_i = dom_ags[di].indices
                idx_j = dom_ags[dj].indices
                dmap = calc_dmap(dom_ags[di], dom_ags[dj])
                u_ah, u_yu = calc_energy(dmap, sig[np.ix_(idx_i, idx_j)],
                                         lam[np.ix_(idx_i, idx_j)],
                                         CUTOFF_LJ, EPS_LJ,
                                         qmap[np.ix_(idx_i, idx_j)],
                                         kappa, rc_yu=CUTOFF_YU)
                ah_matrix[i, j] += np.sum(u_ah)
                yu_matrix[i, j] += np.sum(u_yu)

    ah_matrix /= len(sample_frames)
    yu_matrix /= len(sample_frames)
    ah_matrix += ah_matrix.T
    yu_matrix += yu_matrix.T
    return ah_matrix, yu_matrix, dnames


def _embed_in_canonical(mat, dnames):
    """Embed a construct-specific domain matrix into the canonical 11x11 grid.
    Returns NaN for missing domains (rendered as gray in plots)."""
    full = np.full((N_ALL_DOMAINS, N_ALL_DOMAINS), np.nan)
    idx_map = {d: ALL_DOMAIN_NAMES.index(d) for d in dnames if d in ALL_DOMAIN_NAMES}
    for i, di in enumerate(dnames):
        if di not in idx_map:
            continue
        for j, dj in enumerate(dnames):
            if dj not in idx_map:
                continue
            full[idx_map[di], idx_map[dj]] = mat[i, j]
    return full


# Fixed energy scale for all constructs: ±12 kJ/mol encompasses the
# full range (strongest attraction MD1-MD2 ~ -11 kJ/mol, residual
# boundary repulsion FL KHb-KH8-WWE ~ +12 kJ/mol).
ENERGY_VMIN = -12.0
ENERGY_VMAX = 12.0


def _plot_energy_matrix(ax, mat, title, vmin=ENERGY_VMIN, vmax=ENERGY_VMAX):
    """Plot an 11x11 energy heatmap with gray for NaN (missing domains)."""
    n = mat.shape[0]

    # Gray background for NaN
    gray_bg = np.ones((*mat.shape, 4)) * 0.85
    ax.imshow(gray_bg, origin='lower', aspect='equal')

    # Diverging colormap
    cmap = plt.cm.RdBu_r.copy()
    cmap.set_bad(color=(0.85, 0.85, 0.85, 1.0))
    im = ax.imshow(mat, cmap=cmap, origin='lower', vmin=vmin, vmax=vmax, aspect='equal')

    ax.set_xticks(range(n))
    ax.set_xticklabels(ALL_DOMAIN_NAMES[:n], rotation=45, ha='right', fontsize=7)
    ax.set_yticks(range(n))
    ax.set_yticklabels(ALL_DOMAIN_NAMES[:n], fontsize=7)

    # Annotate non-NaN off-diagonal cells (lower triangle only)
    for i in range(n):
        for j in range(n):
            v = mat[i, j]
            if np.isnan(v) or i == j or j > i:
                continue
            if abs(v) < 0.01:
                continue
            if abs(v) >= 10:
                txt = f'{v:.1f}'
            else:
                txt = f'{v:.2f}'
            norm_v = np.clip((v - vmin) / (vmax - vmin), 0, 1)
            tc = 'white' if norm_v < 0.25 or norm_v > 0.75 else 'black'
            ax.text(j, i, txt, ha='center', va='center', fontsize=5, color=tc)

    plt.colorbar(im, ax=ax, label='kJ/mol', shrink=0.8)
    ax.set_title(title, fontsize=11)
    return im


def run_energy_analysis(active_sets):
    print("\n" + "=" * 70)
    print("ENERGY DECOMPOSITION (calc_energy, trimmed domain boundaries)")
    print("=" * 70)

    # Cache-aware: only compute sets whose energy.npz is missing
    needs_compute = [
        sk for sk in active_sets
        if FORCE_RECOMPUTE or
        cached_data_path(sk, '{set}_energy.npz') is None
    ]
    cached_sets = [sk for sk in active_sets if sk not in needs_compute]
    if cached_sets:
        print(f"  Cached: {', '.join(cached_sets)}")
    if not needs_compute:
        print(f"  All energy data cached (--force to recompute). Skipping.")
        return None
    active_sets = needs_compute

    residues = pd.read_csv(RESIDUES_FILE).set_index('three')

    # Store all canonical-grid matrices for difference maps
    all_ah = {}
    all_yu = {}
    all_total = {}

    for set_key in active_sets:
        info = SETS[set_key]
        domains = get_construct_domains_for_set(set_key)
        dnames = list(domains.keys())
        n_dom = len(dnames)
        print(f"\n--- {info['label']} ({n_dom} domains) ---")

        pdb, dcd = get_sim_paths(set_key, 1, 0)
        u = load_universe(pdb, dcd)
        if u is None:
            print(f"  SKIP: no data"); continue

        sig, lam, qmap, kappa = build_force_maps(u, residues)
        ah_matrix, yu_matrix, dnames = _compute_energy_matrix(
            u, domains, sig, lam, qmap, kappa)
        total = ah_matrix + yu_matrix

        # Save raw (construct-sized) matrices
        np.savez(os.path.join(DATA_PATH, f'{set_key}_energy.npz'),
                 ah=ah_matrix, yu=yu_matrix, total=total, domains=dnames)
        print(f"  Saved: data/{set_key}_energy.npz")

        # Print summary
        for i, di in enumerate(dnames):
            for j, dj in enumerate(dnames):
                if j <= i:
                    continue
                if abs(total[i, j]) > 0.05:
                    print(f"    {di:>8s}-{dj:<8s}: AH={ah_matrix[i,j]:7.2f}  "
                          f"YU={yu_matrix[i,j]:7.2f}  TOT={total[i,j]:7.2f} kJ/mol")

        # Embed into canonical 11x11 grid
        ah_full = _embed_in_canonical(ah_matrix, dnames)
        yu_full = _embed_in_canonical(yu_matrix, dnames)
        tot_full = _embed_in_canonical(total, dnames)
        all_ah[set_key] = ah_full
        all_yu[set_key] = yu_full
        all_total[set_key] = tot_full

        # --- Per-construct 3-panel figure (fixed ±1 kJ/mol scale) ---
        fig, axes = plt.subplots(1, 3, figsize=(20, 6))
        _plot_energy_matrix(axes[0], ah_full, 'Ashbaugh-Hatch (hydrophobic)')
        _plot_energy_matrix(axes[1], yu_full, 'Yukawa (electrostatic)')
        _plot_energy_matrix(axes[2], tot_full, 'Total')
        fig.suptitle(f'{info["label"]} — Interdomain Energy (kJ/mol)\n'
                     f'Blue = attractive, Red = repulsive, Gray = domain not in construct',
                     fontsize=13, y=1.04)
        fig.tight_layout()
        for ext in ['png', 'svg']:
            fig.savefig(os.path.join(FIG_PATH, f'{set_key}_energy.{ext}'),
                        dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved: figures/{set_key}_energy")

    # --- Difference maps: each construct vs FL ---
    if 'fl' in all_total:
        fl_tot = all_total['fl']
        diff_sets = [s for s in active_sets if s != 'fl' and s in all_total]

        if diff_sets:
            n_diff = len(diff_sets)
            fig, axes = plt.subplots(1, n_diff, figsize=(7 * n_diff, 6), squeeze=False)
            for idx, set_key in enumerate(diff_sets):
                ax = axes[0, idx]
                diff = all_total[set_key] - fl_tot
                # NaN where either construct or FL is missing
                _plot_energy_matrix(ax, diff,
                                    f'{SETS[set_key]["label"]} minus FL',
)
            fig.suptitle('Energy Difference vs Full-Length (kJ/mol)\n'
                         'Blue = more attractive in construct, '
                         'Red = more repulsive in construct, Gray = missing',
                         fontsize=13, y=1.06)
            fig.tight_layout()
            for ext in ['png', 'svg']:
                fig.savefig(os.path.join(FIG_PATH, f'energy_diff_vs_fl.{ext}'),
                            dpi=150, bbox_inches='tight')
            plt.close()
            print(f"  Saved: figures/energy_diff_vs_fl")

    # --- All-vs-all pairwise difference maps ---
    computed_sets = [s for s in active_sets if s in all_total]
    n_sets = len(computed_sets)
    if n_sets >= 2:
        fig, axes = plt.subplots(n_sets, n_sets, figsize=(5 * n_sets, 5 * n_sets),
                                 squeeze=False)
        for ri, si in enumerate(computed_sets):
            for ci, sj in enumerate(computed_sets):
                ax = axes[ri][ci]
                if ri == ci:
                    # Diagonal: show the construct's own total energy
                    _plot_energy_matrix(ax, all_total[si],
                                        SETS[si]['label'].split('(')[0].strip(),
    )
                elif ri > ci:
                    # Lower triangle: row minus col
                    diff = all_total[si] - all_total[sj]
                    _plot_energy_matrix(ax, diff,
                                        f'{SETS[si]["label"].split("(")[0].strip()} '
                                        f'- {SETS[sj]["label"].split("(")[0].strip()}',
    )
                else:
                    # Upper triangle: hide
                    ax.set_visible(False)
        fig.suptitle('Interdomain Energy (Total): All vs All (kJ/mol)\n'
                     'Diagonal = absolute, Lower triangle = row minus column\n'
                     'Blue = attractive, Red = repulsive, Gray = missing domain',
                     fontsize=14, y=1.01)
        fig.tight_layout()
        for ext in ['png', 'svg']:
            fig.savefig(os.path.join(FIG_PATH, f'energy_all_vs_all.{ext}'),
                        dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved: figures/energy_all_vs_all")

    # --- All-vs-all: Yukawa (electrostatic) only ---
    if n_sets >= 2:
        yu_vmin, yu_vmax = -1.0, 1.0
        fig, axes = plt.subplots(n_sets, n_sets, figsize=(5 * n_sets, 5 * n_sets),
                                 squeeze=False)
        for ri, si in enumerate(computed_sets):
            for ci, sj in enumerate(computed_sets):
                ax = axes[ri][ci]
                if ri == ci:
                    _plot_energy_matrix(ax, all_yu[si],
                                        SETS[si]['label'].split('(')[0].strip(),
                                        vmin=yu_vmin, vmax=yu_vmax)
                elif ri > ci:
                    diff = all_yu[si] - all_yu[sj]
                    _plot_energy_matrix(ax, diff,
                                        f'{SETS[si]["label"].split("(")[0].strip()} '
                                        f'- {SETS[sj]["label"].split("(")[0].strip()}',
                                        vmin=yu_vmin, vmax=yu_vmax)
                else:
                    ax.set_visible(False)
        fig.suptitle('Yukawa (Electrostatic) Energy: All vs All (kJ/mol)\n'
                     'Diagonal = absolute, Lower triangle = row minus column\n'
                     'Blue = attractive (opposite charges), Red = repulsive (like charges)',
                     fontsize=14, y=1.01)
        fig.tight_layout()
        for ext in ['png', 'svg']:
            fig.savefig(os.path.join(FIG_PATH, f'energy_yukawa_all_vs_all.{ext}'),
                        dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved: figures/energy_yukawa_all_vs_all")

    # --- All-vs-all: Ashbaugh-Hatch (hydrophobic) only ---
    if n_sets >= 2:
        fig, axes = plt.subplots(n_sets, n_sets, figsize=(5 * n_sets, 5 * n_sets),
                                 squeeze=False)
        for ri, si in enumerate(computed_sets):
            for ci, sj in enumerate(computed_sets):
                ax = axes[ri][ci]
                if ri == ci:
                    _plot_energy_matrix(ax, all_ah[si],
                                        SETS[si]['label'].split('(')[0].strip())
                elif ri > ci:
                    diff = all_ah[si] - all_ah[sj]
                    _plot_energy_matrix(ax, diff,
                                        f'{SETS[si]["label"].split("(")[0].strip()} '
                                        f'- {SETS[sj]["label"].split("(")[0].strip()}')
                else:
                    ax.set_visible(False)
        fig.suptitle('Ashbaugh-Hatch (Hydrophobic) Energy: All vs All (kJ/mol)\n'
                     'Diagonal = absolute, Lower triangle = row minus column\n'
                     'Blue = attractive (hydrophobic packing), Red = repulsive (steric)',
                     fontsize=14, y=1.01)
        fig.tight_layout()
        for ext in ['png', 'svg']:
            fig.savefig(os.path.join(FIG_PATH, f'energy_ah_all_vs_all.{ext}'),
                        dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved: figures/energy_ah_all_vs_all")


# ============================================================
# 6. WCN
# ============================================================

def run_wcn_analysis(active_sets):
    print("\n" + "=" * 70)
    print("WEIGHTED COORDINATION NUMBER — parallel")
    print("=" * 70)

    # Per-set caching: each set saved to data/{set}_wcn.npz
    needs_compute = [
        sk for sk in active_sets
        if FORCE_RECOMPUTE or
        cached_data_path(sk, '{set}_wcn.npz') is None
    ]
    cached_sets = [sk for sk in active_sets if sk not in needs_compute]
    if cached_sets:
        print(f"  Cached: {', '.join(cached_sets)}")

    # Compute fresh data for sets that need it
    all_wcn = {}
    if needs_compute:
        print(f"  Computing: {', '.join(needs_compute)}")
        jobs = replicate_jobs(needs_compute)
        results = parallel_map(_worker_wcn, jobs, "wcn")
        by_set = group_by_set(results)

        for set_key in needs_compute:
            info = SETS[set_key]
            reps = by_set.get(set_key, [])
            available = CONSTRUCT_SITES[set_key]
            site_vals = {s: [r['wcn'][s] for r in reps
                              if not np.isnan(r['wcn'].get(s, np.nan))]
                         for s in available}
            all_wcn[set_key] = site_vals
            for s in available:
                v = site_vals[s]
                if v:
                    print(f"  {info['label']} {s}: WCN={np.mean(v):.2f}"
                          f"+/-{np.std(v):.2f}")

            # Save per-set cache
            save_dict = {f'{s}_wcn': np.array(v)
                         for s, v in site_vals.items() if v}
            np.savez(os.path.join(DATA_PATH, f'{set_key}_wcn.npz'),
                     **save_dict)

    # Load cached data for cross-set figure
    for sk in cached_sets:
        p = cached_data_path(sk, '{set}_wcn.npz')
        if not p:
            continue
        cached = np.load(p, allow_pickle=True)
        site_vals = {}
        for key in cached.files:
            if key.endswith('_wcn'):
                site = key[:-4]
                site_vals[site] = cached[key].tolist()
        all_wcn[sk] = site_vals

    # Also write the aggregated (legacy) global file for backward compat
    legacy = {}
    for sk, sw in all_wcn.items():
        for s, v in sw.items():
            if v:
                legacy[f'{sk}_{s}_wcn'] = np.array(v)
    merge_savez(os.path.join(DATA_PATH, 'wcn_active_sites.npz'), legacy, active_sets)

    _plot_violins(active_sets, all_wcn, 'WCN',
                   'Weighted Coordination Number at Active Sites',
                   'wcn_active_sites')
    print(f"  Saved: figures/wcn_active_sites")


# ============================================================
# 7. Active Site Analysis
# ============================================================

def run_active_site_analysis(active_sets):
    print("\n" + "=" * 70)
    print("ACTIVE SITE ANALYSIS — parallel")
    print("=" * 70)

    cache_file = os.path.join(DATA_PATH, 'active_site_stats.npz')
    if os.path.isfile(cache_file) and not FORCE_RECOMPUTE:
        print(f"  Cached: {cache_file} (--force to recompute)")
        return None

    jobs = replicate_jobs(active_sets)
    results = parallel_map(_worker_active_sites, jobs, "active_sites")
    by_set = group_by_set(results)

    all_results = {}
    for set_key in active_sets:
        info = SETS[set_key]
        reps = by_set.get(set_key, [])
        available = CONSTRUCT_SITES[set_key]

        radial = {s: [r['radial'][s] for r in reps if s in r['radial']] for s in available}
        rmsf = {s: [r['rmsf'][s] for r in reps if s in r['rmsf']] for s in available}
        inter_dists = defaultdict(list)
        for r in reps:
            for pair, d in r['inter_dists'].items():
                inter_dists[pair].append(d)

        all_results[set_key] = {'radial': radial, 'rmsf': rmsf,
                                'inter_dists': dict(inter_dists), 'available': available}

        print(f"  {info['label']}: {len(reps)} reps")
        for s in available:
            rp = radial[s]; rf = rmsf[s]
            print(f"    {s}: Radial={np.mean(rp):.3f}+/-{np.std(rp):.3f}, "
                  f"RMSF={np.mean(rf):.3f}+/-{np.std(rf):.3f} nm" if rp else f"    {s}: N/A")

    # Save
    save_dict = {}
    for sk, res in all_results.items():
        for s in res['available']:
            if res['radial'][s]:
                save_dict[f'{sk}_{s}_radial'] = np.array(res['radial'][s])
            if res['rmsf'][s]:
                save_dict[f'{sk}_{s}_rmsf'] = np.array(res['rmsf'][s])
        for pair, vals in res['inter_dists'].items():
            save_dict[f'{sk}_dist_{pair}'] = np.array(vals)
    merge_savez(os.path.join(DATA_PATH, 'active_site_stats.npz'), save_dict, active_sets)

    # RMSF bars
    rmsf_data = {sk: all_results[sk]['rmsf'] for sk in active_sets}
    _plot_violins(active_sets, rmsf_data, 'RMSF (nm)', 'Active Site Catalytic RMSF',
                       'active_site_rmsf')

    # Inter-site distances (violin + points per construct)
    all_pairs = sorted(set(p for r in all_results.values() for p in r['inter_dists']))
    if all_pairs:
        n_pairs = min(len(all_pairs), 6)
        fig, axes = plt.subplots(2, 3, figsize=(15, 8))
        rng = np.random.default_rng(42)
        for pi, pair in enumerate(all_pairs[:n_pairs]):
            ax = axes.flat[pi]
            sd, sl, sc = [], [], []
            for sk in active_sets:
                if pair in all_results[sk]['inter_dists']:
                    sd.append(all_results[sk]['inter_dists'][pair])
                    sl.append(SETS[sk]['label'].split('(')[0].strip())
                    sc.append(SETS[sk]['color'])
            if sd:
                positions = np.arange(len(sd))
                for i, (vals, c, lab) in enumerate(zip(sd, sc, sl)):
                    vals_arr = np.array(vals)
                    if len(vals_arr) >= 2:
                        parts = ax.violinplot([vals_arr], positions=[i], widths=0.6,
                                              showmeans=False, showextrema=False)
                        for body in parts['bodies']:
                            body.set_facecolor(c); body.set_alpha(0.35); body.set_edgecolor(c)
                        m, s = np.mean(vals_arr), np.std(vals_arr)
                        ax.plot([i-0.15, i+0.15], [m, m], color=c, lw=2, zorder=4)
                        ax.plot([i, i], [m-s, m+s], color=c, lw=1.5, zorder=4)
                        ax.plot([i-0.08, i+0.08], [m-s, m-s], color=c, lw=1.5, zorder=4)
                        ax.plot([i-0.08, i+0.08], [m+s, m+s], color=c, lw=1.5, zorder=4)
                    jitter = rng.normal(0, 0.05, len(vals_arr))
                    ax.scatter(i+jitter, vals_arr, c=c, s=15, alpha=0.7,
                               edgecolors='black', linewidths=0.3, zorder=5)
                ax.set_xticks(positions)
                ax.set_xticklabels(sl, fontsize=7, rotation=30, ha='right')
            ax.set_title(pair, fontsize=11); ax.set_ylabel('Distance (nm)')
        for pi in range(n_pairs, 6):
            axes.flat[pi].set_visible(False)
        fig.suptitle('Inter-Active-Site Distances', fontsize=13, y=1.01); fig.tight_layout()
        for ext in ['png', 'svg']:
            fig.savefig(os.path.join(FIG_PATH, f'intersite_distances.{ext}'), dpi=150, bbox_inches='tight')
        plt.close()

        # Secondary version: every pair as one distinctly-colored violin in a
        # single axes sharing one y-axis, so magnitudes compare directly
        # instead of across 6 independently-scaled subplots.
        pair_colors = plt.cm.tab10(np.linspace(0, 1, len(all_pairs[:n_pairs])))
        fig, ax = plt.subplots(figsize=(1.4*n_pairs*max(len(active_sets), 1)+2, 6))
        pos = 0
        xt_pos, xt_lab = [], []
        for pi, pair in enumerate(all_pairs[:n_pairs]):
            pc = pair_colors[pi]
            group_start = pos
            for sk in active_sets:
                vals = all_results[sk]['inter_dists'].get(pair)
                if not vals:
                    continue
                vals_arr = np.array(vals)
                if len(vals_arr) >= 2:
                    parts = ax.violinplot([vals_arr], positions=[pos], widths=0.7,
                                          showmeans=False, showextrema=False)
                    for body in parts['bodies']:
                        body.set_facecolor(pc); body.set_alpha(0.5); body.set_edgecolor(pc)
                    m, s = np.mean(vals_arr), np.std(vals_arr)
                    ax.plot([pos-0.15, pos+0.15], [m, m], color='black', lw=2, zorder=4)
                    ax.plot([pos, pos], [m-s, m+s], color='black', lw=1.2, zorder=4)
                jitter = rng.normal(0, 0.05, len(vals_arr))
                ax.scatter(pos+jitter, vals_arr, c=[pc], s=12, alpha=0.6,
                           edgecolors='black', linewidths=0.3, zorder=5)
                lab = pair if len(active_sets) == 1 else f"{pair} ({SETS[sk]['label'].split('(')[0].strip()})"
                xt_pos.append(pos); xt_lab.append(lab)
                pos += 1
            if pos > group_start:
                pos += 0.5  # gap between pair groups
        ax.set_xticks(xt_pos); ax.set_xticklabels(xt_lab, fontsize=8, rotation=30, ha='right')
        ax.set_ylabel('Distance (nm)')
        ax.set_title('Inter-Active-Site Distances (all pairs, shared y-axis)', fontsize=13)
        fig.tight_layout()
        for ext in ['png', 'svg']:
            fig.savefig(os.path.join(FIG_PATH, f'intersite_distances_v2.{ext}'), dpi=150, bbox_inches='tight')
        plt.close()

    # Radial position violin + points
    radial_data = {sk: all_results[sk]['radial'] for sk in active_sets}
    _plot_violins(active_sets, radial_data, 'Dist-to-COM / Rg',
                  'Active Site Radial Position (>1 = surface)',
                  'active_site_radial', hline=1.0)

    # Radial heatmap (secondary summary)
    fig, ax = plt.subplots(figsize=(8, 5))
    mat = np.full((len(active_sets), len(SITE_NAMES)), np.nan)
    for si, sk in enumerate(active_sets):
        for sj, sn in enumerate(SITE_NAMES):
            v = all_results[sk]['radial'].get(sn, [])
            if v:
                mat[si, sj] = np.mean(v)
    im = ax.imshow(mat, cmap='RdYlBu_r', aspect='auto', vmin=0.5, vmax=1.5)
    ax.set_yticks(range(len(active_sets)))
    ax.set_yticklabels([SETS[s]['label'] for s in active_sets])
    ax.set_xticks(range(len(SITE_NAMES))); ax.set_xticklabels(SITE_NAMES)
    for i in range(len(active_sets)):
        for j in range(len(SITE_NAMES)):
            v = mat[i,j]; ax.text(j, i, f'{v:.2f}' if not np.isnan(v) else '—',
                                   ha='center', va='center', fontsize=10)
    plt.colorbar(im, ax=ax, label='Dist-to-COM / Rg', shrink=0.8)
    ax.set_title('Active Site Radial Position (>1 = surface)')
    fig.tight_layout()
    for ext in ['png', 'svg']:
        fig.savefig(os.path.join(FIG_PATH, f'active_site_summary.{ext}'), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: figures/active_site_rmsf, active_site_radial, intersite_distances(_v2), active_site_summary")

    return all_results


# ============================================================
# 8. Accessibility (SAA)
# ============================================================

def run_accessibility_analysis(active_sets):
    print("\n" + "=" * 70)
    print("STERIC ACCESSIBILITY (SAA) — parallel")
    print("=" * 70)
    print(f"  Probe={PROBE_RADIUS}nm, MaxDist={MAX_DIST}nm, Rays={N_RAYS}, SeqSep=+/-{SEQ_SEP}")

    cache_file = os.path.join(DATA_PATH, 'accessibility_stats.npz')
    if os.path.isfile(cache_file) and not FORCE_RECOMPUTE:
        print(f"  Cached: {cache_file} (--force to recompute)")
        return None

    jobs = replicate_jobs(active_sets)
    results = parallel_map(_worker_accessibility, jobs, "accessibility")
    by_set = group_by_set(results)

    all_stats = {}
    for set_key in active_sets:
        info = SETS[set_key]
        reps = by_set.get(set_key, [])
        available = CONSTRUCT_SITES[set_key]

        saa = {s: [r['saa'][s] for r in reps if not np.isnan(r['saa'].get(s, np.nan))]
               for s in available}
        cone = {s: [r['cone'][s] for r in reps if not np.isnan(r['cone'].get(s, np.nan))]
                for s in available}
        shell = {s: [r['shell'][s] for r in reps if not np.isnan(r['shell'].get(s, np.nan))]
                 for s in available}
        all_stats[set_key] = {'saa': saa, 'cone': cone, 'shell': shell}

        print(f"  {info['label']}: {len(reps)} reps")
        for s in available:
            if saa[s]:
                print(f"    {s}: SAA={np.mean(saa[s]):.3f}+/-{np.std(saa[s]):.3f}  "
                      f"Cone={np.mean(cone[s]):.1f}+/-{np.std(cone[s]):.1f}deg  "
                      f"Shell={np.mean(shell[s]):.4f}+/-{np.std(shell[s]):.4f}")

    # Save
    save_dict = {}
    for sk, st in all_stats.items():
        for sn in CONSTRUCT_SITES[sk]:
            for metric in ['saa', 'cone', 'shell']:
                v = st[metric][sn]
                if v:
                    save_dict[f'{sk}_{sn}_{metric}'] = np.array(v)
    merge_savez(os.path.join(DATA_PATH, 'accessibility_stats.npz'), save_dict, active_sets)

    # SAA bars
    _plot_violins(active_sets, {sk: all_stats[sk]['saa'] for sk in active_sets},
                       'SAA (fraction open)', 'Solid Angle Accessibility', 'accessibility_saa',
                       hline=0.5)

    # Cone bars
    _plot_violins(active_sets, {sk: all_stats[sk]['cone'] for sk in active_sets},
                       'Max cone angle (deg)', 'Maximum Approach Cone Angle', 'accessibility_cone',
                       hline=90)

    # Summary heatmap
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    for pi, (metric, label, cm, vr) in enumerate([
        ('saa', 'SAA', 'RdYlGn', (0, 0.8)),
        ('cone', 'Cone (deg)', 'RdYlGn', (30, 100)),
        ('shell', 'Shell density', 'RdYlGn_r', (0, 0.05)),
    ]):
        ax = axes[pi]
        mat = np.full((len(active_sets), len(SITE_NAMES)), np.nan)
        for si, sk in enumerate(active_sets):
            for sj, sn in enumerate(SITE_NAMES):
                v = all_stats[sk][metric].get(sn, [])
                if v:
                    mat[si, sj] = np.mean(v)
        im = ax.imshow(mat, cmap=cm, aspect='auto', vmin=vr[0], vmax=vr[1])
        ax.set_yticks(range(len(active_sets)))
        ax.set_yticklabels([SETS[s]['label'] for s in active_sets], fontsize=9)
        ax.set_xticks(range(len(SITE_NAMES))); ax.set_xticklabels(SITE_NAMES, fontsize=11)
        for i in range(len(active_sets)):
            for j in range(len(SITE_NAMES)):
                v = mat[i,j]
                ax.text(j, i, f'{v:.3f}' if not np.isnan(v) else '—', ha='center', va='center', fontsize=9)
        plt.colorbar(im, ax=ax, label=label, shrink=0.8); ax.set_title(label)
    fig.suptitle('Steric Accessibility Summary', fontsize=14, y=1.02); fig.tight_layout()
    for ext in ['png', 'svg']:
        fig.savefig(os.path.join(FIG_PATH, f'accessibility_summary.{ext}'), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: figures/accessibility_saa, accessibility_cone, accessibility_summary")

    return all_stats


# ============================================================
# Shared plotting helper
# ============================================================

def _plot_violins(active_sets, data_dict, ylabel, title, fname, hline=None):
    """Violin + jitter points + mean/SD annotation per site per construct.
    data_dict[set_key][site_name] = list of per-replicate values."""
    n_sites = len(SITE_NAMES)
    n_sets = len(active_sets)
    width = 0.7 / n_sets  # violin half-width spacing
    fig, ax = plt.subplots(figsize=(max(10, 2.5*n_sites), 5.5))
    x = np.arange(n_sites)
    rng = np.random.default_rng(42)

    for si, sk in enumerate(active_sets):
        sd = data_dict.get(sk, {})
        color = SETS[sk]['color']
        offset = (si - n_sets/2 + 0.5) * width
        positions_for_set = x + offset

        for sj, sn in enumerate(SITE_NAMES):
            vals = sd.get(sn, [])
            if not vals or len(vals) < 2:
                # Single point or missing: just plot the point
                if vals:
                    ax.scatter(positions_for_set[sj], vals[0], c=color, s=30, zorder=5,
                               edgecolors='black', linewidths=0.5)
                continue

            vals_arr = np.array(vals)
            pos = positions_for_set[sj]

            # Violin
            parts = ax.violinplot([vals_arr], positions=[pos], widths=width*0.85,
                                   showmeans=False, showextrema=False)
            for body in parts['bodies']:
                body.set_facecolor(color)
                body.set_alpha(0.35)
                body.set_edgecolor(color)

            # Mean + SD bar
            m, s = np.mean(vals_arr), np.std(vals_arr)
            ax.plot([pos - width*0.2, pos + width*0.2], [m, m],
                    color=color, lw=2, zorder=4)
            ax.plot([pos, pos], [m - s, m + s], color=color, lw=1.5, zorder=4)
            ax.plot([pos - width*0.1, pos + width*0.1], [m - s, m - s],
                    color=color, lw=1.5, zorder=4)
            ax.plot([pos - width*0.1, pos + width*0.1], [m + s, m + s],
                    color=color, lw=1.5, zorder=4)

            # Jitter points
            jitter = rng.normal(0, width*0.12, len(vals_arr))
            ax.scatter(pos + jitter, vals_arr, c=color, s=18, alpha=0.7, zorder=5,
                       edgecolors='black', linewidths=0.4)

        # Legend entry (invisible scatter for label)
        ax.scatter([], [], c=color, s=40, label=SETS[sk]['label'], edgecolors='black', linewidths=0.5)

    ax.set_xticks(x); ax.set_xticklabels(SITE_NAMES, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=11); ax.set_title(title, fontsize=13)
    if hline is not None:
        ax.axhline(hline, color='gray', ls=':', lw=1, label=f'{hline}')
    ax.legend(fontsize=9, loc='best'); fig.tight_layout()
    for ext in ['png', 'svg']:
        fig.savefig(os.path.join(FIG_PATH, f'{fname}.{ext}'), dpi=150, bbox_inches='tight')
    plt.close()


# ============================================================
# Main
# ============================================================

def main():
    parser = ArgumentParser(description='PARP14 trajectory analysis (parallelized)')
    parser.add_argument('--workers', type=int, default=0,
                        help='Number of parallel workers (0 = auto-detect, 1 = sequential)')
    parser.add_argument('--conf-prop', action='store_true')
    parser.add_argument('--dmap', action='store_true')
    parser.add_argument('--cmap', action='store_true')
    parser.add_argument('--fnc', action='store_true')
    parser.add_argument('--energy', action='store_true')
    parser.add_argument('--wcn', action='store_true')
    parser.add_argument('--active-sites', action='store_true')
    parser.add_argument('--accessibility', action='store_true')
    parser.add_argument('--set', nargs='+', default=None,
                        help='Sets to analyze. Can be: fl/md/core/mka/norrm/'
                             'noart/md3art/fl_optimized, fragment name '
                             "(e.g. md1l1_md2), or 'all' for everything")
    parser.add_argument('--include-fragments', action='store_true',
                        help='Auto-discover and include all fragments/')
    parser.add_argument('--sim-folder', nargs='+', default=None, metavar='PATH',
                        help='One or more simulation-set folders to analyze '
                             '(each containing seed-*_sample-*/ replicates). '
                             "Domain units are read from the folder's metadata.json, "
                             'or pass --units. Use this for new simulations '
                             'without editing the script.')
    parser.add_argument('--units', nargs='+', default=None, metavar='UNIT',
                        help='FL domain units in the --sim-folder construct '
                             '(e.g. md1l1 md2 md3). Needed only when the folder '
                             'has no metadata.json.')
    parser.add_argument('--force', '--force-recompute', dest='force',
                        action='store_true',
                        help='Force re-computation of cached data files. '
                             'Default: skip data computation when .npz/.npy '
                             'files already exist; just rebuild figures.')
    args = parser.parse_args()

    # Expose force flag globally for run_* functions to check
    global FORCE_RECOMPUTE
    FORCE_RECOMPUTE = args.force

    # Auto-discover fragments only if --include-fragments is set.
    # `--set all` without --include-fragments runs only the 8 NAMED sets;
    # `--set all --include-fragments` runs named sets + all fragments.
    # Specific `--set <fragment_name>` works without --include-fragments.
    # Register any folders passed via --sim-folder (before the worker pool is
    # created, so forked workers inherit the registration).
    folder_keys = []
    if args.sim_folder:
        for folder in args.sim_folder:
            k = register_sim_folder(folder, units=args.units)
            folder_keys.append(k)
            print(f"  Registered sim folder: {folder} -> set '{k}' "
                  f"(units: {', '.join(CONSTRUCT_UNITS[k])})")

    if args.include_fragments:
        discovered = discover_fragments()
        print(f"  Discovered {len(discovered)} fragment simulations")
    elif args.set:
        # User may have named specific fragments — register those (no full
        # auto-discovery scan).
        for s in args.set:
            if s in SETS:
                continue
            frag_dir = os.path.join(CWD, 'fragments', s)
            meta_file = os.path.join(frag_dir, 'metadata.json')
            if os.path.isdir(frag_dir) and os.path.isfile(meta_file):
                with open(meta_file) as f:
                    register_fragment(s, json.load(f))

    global N_WORKERS
    if args.workers == 0:
        N_WORKERS = min(os.cpu_count() or 4, 25)
    else:
        N_WORKERS = args.workers

    os.makedirs(DATA_PATH, exist_ok=True)
    os.makedirs(FIG_PATH, exist_ok=True)

    any_flag = (args.conf_prop or args.dmap or args.cmap or args.fnc or
                args.energy or args.wcn or args.active_sites or args.accessibility)
    run_all = not any_flag

    # Resolve active_sets: handle 'all', bare fragment names, normal sets
    if args.set:
        if 'all' in args.set:
            # 'all' = everything currently registered in SETS.
            # If --include-fragments was passed, that includes fragments.
            # Otherwise, only the 8 named sets.
            active_sets = list(SETS.keys())
        else:
            active_sets = []
            for s in args.set:
                if s in SETS:
                    active_sets.append(s)
                elif f'frag_{s}' in SETS:
                    active_sets.append(f'frag_{s}')
                else:
                    print(f"  WARN: unknown set '{s}', skipping")
        # --sim-folder always adds its folders on top of --set selection
        active_sets += [k for k in folder_keys if k not in active_sets]
    elif folder_keys:
        # Only --sim-folder given: analyze just those folders
        active_sets = folder_keys
    else:
        # Nothing specified: default = all named/registered sets
        active_sets = list(SETS.keys())

    print("=" * 70)
    print("PARP14 Trajectory Analysis (parallelized)")
    print("=" * 70)
    print(f"  Sets: {', '.join(active_sets)}")
    print(f"  Workers: {N_WORKERS}")
    print(f"  Equilibration: {SKIP_FRAMES} frames ({SKIP_FRAMES*0.01:.1f} ns)")
    print(f"  Output: data/ + figures/")

    for sk in active_sets:
        detect_n_frames(sk)

    if run_all or args.conf_prop:
        run_conf_prop(active_sets)
    if run_all or args.dmap:
        run_dmap_analysis(active_sets)
    if run_all or args.cmap:
        run_cmap_analysis(active_sets)
    if run_all or args.fnc:
        run_fnc_analysis(active_sets)
    if run_all or args.energy:
        run_energy_analysis(active_sets)
    if run_all or args.wcn:
        run_wcn_analysis(active_sets)
    if run_all or args.active_sites:
        run_active_site_analysis(active_sets)
    if run_all or args.accessibility:
        run_accessibility_analysis(active_sets)

    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print(f"  Data: {DATA_PATH}")
    print(f"  Figures: {FIG_PATH}")
    print("=" * 70)


if __name__ == '__main__':
    main()
