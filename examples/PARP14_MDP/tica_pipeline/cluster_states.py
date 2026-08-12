#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cluster CALVADOS CG-MD frames into conformational states via TICA + MSM +
PCCA+ (or geometric K-means as a baseline).

Pipeline (with --reduce tica --pcca):
  1. Extract features per frame (--features mode; see below)
  2. Standardize (z-score per feature)
  3. ITS sweep across multiple lags (PyEMMA-style validation)
  4. TICA at chosen --tica-lag
  5. K-means into N microstates (default 200; --n-microstates)
  6. Build MSM at --msm-lag
  7. PCCA+ coarse-grain to --k macrostates
  8. Optional Chapman-Kolmogorov validation (--ck-test)

Featurization (--features <mode>):
  com           inter-domain COM-COM distances (default, simple baseline)
  ca            pairwise CA-CA distances within each domain at stride
                --ca-stride (intermediate detail)
  linker_ca     positions (xyz, COM-centered) of unrestrained linker CAs.
                Captures the SLOW dofs in restrained CG-MD.
  orient        inter-domain rigid-body rotation matrices (9 elements per
                pair). Captures hinge rotations between domains.
  interface_ca  CA-CA distances between edge residues of adjacent
                domains. Captures hinge-contact making/breaking.
  segments      pairwise COM-COM distances between consecutive
                --ca-stride-residue chain windows. Coarse global view.
  torsion       CA pseudo-dihedral angles within each domain (mostly
                frozen by restraints; usually scores low).

Domain selection (--domains):
  Default = all 11 grouped units (RRM1, RRM2, RRM3, KH1-KH6, KH7a,
            MD1L1, MD2, MD3, KHb-KH8, WWE, ART). Auto-filters to those
            present in the set's construct.

Outputs (per run, in figures/05_clustering/<date>/<set>_<feat>_<reduce>/):
    acf.png              feature autocorrelation diagnostic
    its.png              MSM ITS validation
    tica_spectrum.png    TICA eigenvalues + implied timescales
    silhouette.png       k-sweep cluster validity
    pca_clusters.png     2D PCA/TICA scatter colored by state
    distance_profile.png per-state mean features (heatmap)
    population.png       state populations
    replicate_state.png  per-replicate state distribution
    microstates_landscape.png  (if --pcca) microstate centers on landscape
    pcca_decomposition.png     (if --pcca) microstate→macrostate grouping
    ck_test.png          (if --ck-test) Chapman-Kolmogorov validation

Plus representative_frames/<date>/<run_tag>/state_{1..k}.pdb

Usage examples:
    python cluster_states.py --set fl_optimized
    python cluster_states.py --set md --reduce tica --pcca --ck-test
    python cluster_states.py --set fl_optimized --features linker_ca \\
                              --reduce tica --pcca --k 5
    python cluster_states.py --set fragment_name --reduce tica
"""

import os
# Cap BLAS threads BEFORE numpy/sklearn imports — many-core machines (>128
# cores) trigger OpenBLAS "too many memory regions" segfault otherwise.
os.environ.setdefault('OPENBLAS_NUM_THREADS', '8')
os.environ.setdefault('OMP_NUM_THREADS', '8')
os.environ.setdefault('MKL_NUM_THREADS', '8')
os.environ.setdefault('NUMEXPR_NUM_THREADS', '8')

import json
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

# ============================================================
# Configuration
# ============================================================

CWD = Path(__file__).resolve().parent.parent
from datetime import date as _date
DATA_PATH = CWD / 'data'
FIG_ROOT = CWD / 'figures'
REP_PATH = CWD / 'representative_frames'
for p in (DATA_PATH, FIG_ROOT, REP_PATH):
    p.mkdir(exist_ok=True)

# FIG_PATH is computed PER-RUN inside main() based on set + featurization
# so different featurizations don't overwrite each other:
#   figures/05_clustering/2026-05-18/<set>_<feat>_<reduce>/
FIG_PATH = FIG_ROOT  # placeholder; overwritten in main()

DOMAINS_FL = {
    'RRM1':    (6, 88),
    'RRM2':    (146, 224), 'RRM3': (225, 314),
    'KH1-KH6': (315, 737), 'KH7a': (738, 789),
    'MD1L1':   (790, 1004), 'MD2': (1005, 1193), 'MD3': (1207, 1388),
    'KHb-KH8': (1389, 1533), 'WWE': (1534, 1602), 'ART': (1603, 1801),
}

# 11 grouped domain units used to define constructs (lowercase FL ranges)
DOMAIN_UNITS = {
    'rrm1':    (1, 145),    'rrm2':    (146, 224), 'rrm3':    (225, 314),
    'kh1-kh6': (315, 737),  'kh7a':    (738, 789),
    'md1l1':   (790, 1004), 'md2':     (1004, 1193), 'md3':   (1207, 1388),
    'khb-kh8': (1389, 1533), 'wwe':    (1534, 1602), 'art':   (1603, 1801),
}

# Which units are in each named construct
CONSTRUCT_UNITS = {
    'fl':           list(DOMAIN_UNITS.keys()),
    'fl_optimized': list(DOMAIN_UNITS.keys()),
    'md':           ['md1l1', 'md2', 'md3'],
    'core':         ['kh7a', 'md1l1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art'],
    'mka':          ['md1l1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art'],
    'norrm':        ['kh1-kh6', 'kh7a', 'md1l1', 'md2', 'md3', 'khb-kh8',
                     'wwe', 'art'],
    'noart':        ['kh1-kh6', 'kh7a', 'md1l1', 'md2', 'md3', 'khb-kh8',
                     'wwe'],
    'md3art':       ['md3', 'khb-kh8', 'wwe', 'art'],
}

# Contiguous constructs: built from a single uninterrupted FL span (NO inter-domain
# residues deleted), so FL->local is one offset (local = FL - (fl_lo - 1)). Unlike the
# unit-concatenated constructs above, these keep inter-domain linkers (e.g. md_full
# spans FL 790-1388 incl. the 1194-1206 MD2-MD3 linker the old `md` set dropped).
CONTIGUOUS_CONSTRUCTS = {
    'md_full': (790, 1388),
}

# Map UI domain names (uppercase from --domains arg) → unit key
DOMAIN_TO_UNIT = {
    'RRM1': 'rrm1', 'RRM2': 'rrm2', 'RRM3': 'rrm3',
    'KH1-KH6': 'kh1-kh6', 'KH7a': 'kh7a',
    'MD1L1': 'md1l1', 'MD2': 'md2', 'MD3': 'md3',
    'KHb-KH8': 'khb-kh8', 'WWE': 'wwe', 'ART': 'art',
}

# The full 11-unit default (matches --domains' argparse default). Any --domains
# subset that differs from this gets a fingerprint suffix on its cache file /
# output dir so it can never collide with another subset run at the same
# features_mode on the same set (see _domains_suffix below).
_FULL_DOMAINS = ['RRM1', 'RRM2', 'RRM3', 'KH1-KH6', 'KH7a',
                 'MD1L1', 'MD2', 'MD3', 'KHb-KH8', 'WWE', 'ART']


def _domains_suffix(domains):
    if set(domains) == set(_FULL_DOMAINS):
        return ''
    import hashlib
    h = hashlib.sha1('_'.join(sorted(domains)).encode()).hexdigest()[:6]
    return f'-d{h}'

# Directories registered via --sim-folder: set_key -> absolute folder Path.
# Populated in main() before the worker pool is created so forked workers inherit it.
EXTERNAL_DIRS = {}


def register_sim_folder(path, units=None):
    """Register an arbitrary simulation folder (from --sim-folder) as a set so it
    can be clustered without editing this script. Domain units are read from the
    folder's metadata.json, the `units` arg, or (if the basename is a known set)
    that set. Returns the set_key under which it is registered.

    Also registers with sim_registry (reg.SETS) so reg.get_active_sites_for_set()
    works for --features site_orient (active-site pocket lookups).
    """
    folder = Path(path).resolve()
    if not folder.is_dir():
        raise FileNotFoundError(f"sim folder not found: {folder}")
    set_key = reg.register_external_folder(str(folder), units=units)
    CONSTRUCT_UNITS[set_key] = reg.SETS[set_key]['units']
    EXTERNAL_DIRS[set_key] = folder
    return set_key


SEEDS = range(1, 6)
SAMPLES = range(0, 5)
SKIP_FRAMES = 50

# Cache: set_key -> replicate dirs discovered by _flat_replicate_dirs() below.
_FLAT_DIRS_CACHE = {}


def _flat_replicate_dirs(folder):
    """Every immediate subdirectory of `folder` containing a .dcd file.

    Fallback for --sim-folder layouts that aren't a seed-{1-5}_sample-{0-4}
    grid -- e.g. a handful of independent long runs named arbitrarily
    (fl_go's state-*_tica_seed-*_sample-*_fr*/ dirs)."""
    dirs = []
    for entry in sorted(os.listdir(folder)):
        d = folder / entry
        if d.is_dir() and list(d.glob('*.dcd')):
            dirs.append(d)
    return dirs


def resolve_sim_dir(set_key, seed, sample):
    """Replicate directory for any set_key, falling back to flat discovery
    for --sim-folder sets that don't follow the seed-N_sample-M grid."""
    if set_key in EXTERNAL_DIRS:
        base = EXTERNAL_DIRS[set_key]
    elif set_key.startswith('frag_'):
        base = CWD / 'fragments' / set_key[5:]
    else:
        base = CWD / set_key
    sim_dir = base / f'seed-{seed}_sample-{sample}'
    if set_key in EXTERNAL_DIRS and not sim_dir.is_dir():
        if set_key not in _FLAT_DIRS_CACHE:
            _FLAT_DIRS_CACHE[set_key] = _flat_replicate_dirs(base)
        flat = _FLAT_DIRS_CACHE[set_key]
        idx = (seed - 1) * len(SAMPLES) + sample
        if idx < len(flat):
            sim_dir = flat[idx]
    return sim_dir


def resolve_topology_pdb(sim_dir):
    """Topology PDB for one replicate dir, preferring top.pdb but falling
    back to restart.pdb / checkpoint.pdb (long fresh-start runs like fl_go
    never write a top.pdb)."""
    for name in ('top.pdb', 'restart.pdb', 'checkpoint.pdb'):
        p = sim_dir / name
        if p.is_file():
            return p
    return sim_dir / 'top.pdb'


def get_construct_domain_ranges(set_key, requested_domains):
    """Return {domain_name: (start, end)} in CONSTRUCT-LOCAL numbering for
    the requested domains that exist in this set.

    Filters out domains absent from the construct.
    For fl/fl_optimized, returns FL numbering directly.
    For fragments (frag_NAME), reads metadata.json.
    """
    if set_key in ('fl', 'fl_optimized'):
        return {d: DOMAINS_FL[d] for d in requested_domains
                if d in DOMAINS_FL}

    # Fragment: read metadata.json
    if set_key.startswith('frag_'):
        frag_name = set_key[5:]
        meta_file = CWD / 'fragments' / frag_name / 'metadata.json'
        if not meta_file.exists():
            return {}
        with open(meta_file) as f:
            meta = json.load(f)
        # metadata has 'domain_ranges_construct' parallel to 'restraint_labels'
        labels = meta.get('restraint_labels', [])
        ranges = meta.get('domain_ranges_construct', [])
        # Map restraint labels (MD1L1, MD2, ...) back to construct positions
        # Note: restraint labels may include sub-domains like KH1, KH2 etc
        # — for simplicity we use the requested_domains as-is and look them up
        label_to_range = dict(zip(labels, ranges))
        result = {}
        for d in requested_domains:
            # Try the exact uppercase label first
            if d in label_to_range:
                result[d] = tuple(label_to_range[d])
                continue
            # Special case: MD1L1 vs MD1
            if d == 'MD1L1' and 'MD1L1' in label_to_range:
                result[d] = tuple(label_to_range['MD1L1'])
                continue
            # Fall through: try to compute from FL mapping below
        if result:
            return result
        # else fall through to FL-block mapping using fragment's units
        units = meta.get('units', [])
        return _map_fl_to_construct(units, requested_domains)

    # Contiguous construct (md_full): single FL offset, linkers retained.
    if set_key in CONTIGUOUS_CONSTRUCTS:
        fl_lo, fl_hi = CONTIGUOUS_CONSTRUCTS[set_key]
        result = {}
        for d in requested_domains:
            unit = DOMAIN_TO_UNIT.get(d)
            if unit is None or unit not in DOMAIN_UNITS:
                continue
            us, ue = DOMAIN_UNITS[unit]
            s, e = max(us, fl_lo), min(ue, fl_hi)
            if s > e:
                continue
            result[d] = (s - (fl_lo - 1), e - (fl_lo - 1))
        return result

    # Named construct (md, core, mka, norrm, noart, md3art)
    units = CONSTRUCT_UNITS.get(set_key)
    if not units:
        return {}
    return _map_fl_to_construct(units, requested_domains)


def _map_fl_to_construct(units, requested_domains):
    """Build FL→construct residue map and return construct-local ranges for
    requested domains that have any overlap with the construct."""
    # Compute continuous FL blocks
    ranges = sorted([DOMAIN_UNITS[u] for u in units if u in DOMAIN_UNITS])
    if not ranges:
        return {}
    merged = [list(ranges[0])]
    for s, e in ranges[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])

    # Build FL→construct mapping
    segments = []
    pos = 1
    for fl_s, fl_e in merged:
        segments.append((fl_s, fl_e, pos - fl_s))
        pos += (fl_e - fl_s + 1)

    def fl_to_c(fl):
        for fl_s, fl_e, off in segments:
            if fl_s <= fl <= fl_e:
                return fl + off
        return None

    # FL coverage as set
    fl_set = set()
    for fl_s, fl_e in merged:
        fl_set.update(range(fl_s, fl_e + 1))

    result = {}
    for d in requested_domains:
        if d not in DOMAINS_FL:
            continue
        fl_s, fl_e = DOMAINS_FL[d]
        # Clip to construct coverage
        while fl_s <= fl_e and fl_s not in fl_set:
            fl_s += 1
        while fl_e >= fl_s and fl_e not in fl_set:
            fl_e -= 1
        if fl_s > fl_e:
            continue
        c_s = fl_to_c(fl_s)
        c_e = fl_to_c(fl_e)
        if c_s is not None and c_e is not None and c_e - c_s + 1 >= 5:
            result[d] = (c_s, c_e)
    return result


# ============================================================
# Feature extraction
# ============================================================

def extract_features_one_replicate(args):
    """Compute inter-domain features for one replicate.

    Returns dict with:
      'features':   (n_frames, n_features) array
      'metadata':   per-frame (seed, sample, frame_idx) tuples
      'feature_names': list of feature labels
    """
    (set_key, seed, sample, domains, add_rg,
     features_mode, ca_stride) = args
    import MDAnalysis as mda

    sim_dir = resolve_sim_dir(set_key, seed, sample)
    pdb = resolve_topology_pdb(sim_dir)
    # CALVADOS uses sysname as DCD filename — try a few possibilities
    candidate_dcds = [
        sim_dir / 'parp14.dcd',
        sim_dir / 'parp14_macrodomains.dcd',
        sim_dir / 'parp14_core.dcd',
        sim_dir / 'parp14_mka.dcd',
        sim_dir / 'parp14_norrm.dcd',
        sim_dir / 'parp14_noart.dcd',
        sim_dir / 'parp14_md3art.dcd',
    ]
    dcd = next((p for p in candidate_dcds if p.is_file()), None)
    if dcd is None:
        # Fall back: any .dcd in the directory
        dcds = list(sim_dir.glob('*.dcd'))
        dcd = dcds[0] if dcds else None
    if not pdb.is_file() or dcd is None:
        return None

    # Map requested FL domain names to construct-local residue ranges
    construct_ranges = get_construct_domain_ranges(set_key, domains)
    if not construct_ranges:
        return None

    u = mda.Universe(str(pdb), str(dcd))
    ags = {}
    for dn, (ds, de) in construct_ranges.items():
        ag = u.select_atoms(f'resid {ds}:{de}')
        if len(ag) > 0:
            ags[dn] = ag
    if not ags:
        return None

    # Build pair list of domains
    domain_list = [d for d in domains if d in ags]
    pairs = [(domain_list[i], domain_list[j])
             for i in range(len(domain_list))
             for j in range(i + 1, len(domain_list))]

    rows = []
    metadata = []
    feature_names = None

    if features_mode == 'ca':
        # Per-domain CA selections, strided
        ca_sel = {}
        for d in domain_list:
            ag = ags[d]
            # AG is already a residue-by-residue CA list
            idx_sub = list(range(0, len(ag), ca_stride))
            ca_sel[d] = (ag, idx_sub)

        # Build feature names (compute once)
        names = []
        for a, b in pairs:
            for ia in ca_sel[a][1]:
                ra = ca_sel[a][0].residues[ia].resid
                for ib in ca_sel[b][1]:
                    rb = ca_sel[b][0].residues[ib].resid
                    names.append(f'ca_{a}_{ra}__{b}_{rb}')
        if add_rg:
            names += [f'rg_{d}' for d in domain_list]
        feature_names = names

        for frame_idx, ts in enumerate(u.trajectory[SKIP_FRAMES:]):
            row = []
            # Cache positions per domain (in nm)
            pos_cache = {d: ca_sel[d][0].positions / 10.0
                         for d in domain_list}
            for a, b in pairs:
                pa = pos_cache[a][ca_sel[a][1]]
                pb = pos_cache[b][ca_sel[b][1]]
                diff = pa[:, None, :] - pb[None, :, :]
                dists = np.sqrt(np.sum(diff ** 2, axis=2))
                row.extend(dists.ravel().tolist())
            if add_rg:
                for d in domain_list:
                    p = pos_cache[d]
                    com = p.mean(axis=0)
                    rg = float(np.sqrt(np.mean(
                        np.sum((p - com) ** 2, axis=1))))
                    row.append(rg)
            rows.append(row)
            metadata.append((seed, sample, frame_idx + SKIP_FRAMES))
    elif features_mode == 'torsion':
        # CA pseudo-dihedral angles between 4 consecutive CAs.
        # For N residues in the domain, we get N-3 dihedrals.
        # Each dihedral is encoded as (sin, cos) to avoid wrapping.
        # We compute torsions within each domain only.
        names = []
        for d in domain_list:
            n_res = len(ags[d])
            for i in range(n_res - 3):
                names.append(f'tor_{d}_{i+1}_sin')
                names.append(f'tor_{d}_{i+1}_cos')
        if add_rg:
            names += [f'rg_{d}' for d in domain_list]
        feature_names = names

        for frame_idx, ts in enumerate(u.trajectory[SKIP_FRAMES:]):
            row = []
            for d in domain_list:
                pos = ags[d].positions / 10.0  # nm
                # CA-CA-CA-CA dihedrals
                p1 = pos[:-3]
                p2 = pos[1:-2]
                p3 = pos[2:-1]
                p4 = pos[3:]
                b1 = p2 - p1
                b2 = p3 - p2
                b3 = p4 - p3
                n1 = np.cross(b1, b2)
                n2 = np.cross(b2, b3)
                b2n = b2 / (np.linalg.norm(b2, axis=1, keepdims=True) + 1e-12)
                m = np.cross(n1, b2n)
                x = np.einsum('ij,ij->i', n1, n2)
                y = np.einsum('ij,ij->i', m, n2)
                phi = np.arctan2(y, x)
                # Interleave sin, cos
                for v in phi:
                    row.append(float(np.sin(v)))
                    row.append(float(np.cos(v)))
            if add_rg:
                for d in domain_list:
                    p = ags[d].positions / 10.0
                    com = p.mean(axis=0)
                    rg = float(np.sqrt(np.mean(
                        np.sum((p - com) ** 2, axis=1))))
                    row.append(rg)
            rows.append(row)
            metadata.append((seed, sample, frame_idx + SKIP_FRAMES))

    elif features_mode == 'linker_ca':
        # Positions of all CA atoms that are NOT inside any restraint
        # domain — the actual unrestrained slow degrees of freedom.
        # Encoded as CA positions relative to chain COM (translation-
        # invariant after standardization).
        all_ca = u.select_atoms('name CA')
        n_total = len(all_ca)

        # Find the ACTUAL restraint ranges from the construct's
        # domains.yaml. construct_ranges from get_construct_domain_ranges
        # uses full domain extents which may cover the entire chain (no
        # linkers). The yaml has the trimmed restraint cores.
        import yaml as _yaml
        if set_key.startswith('frag_'):
            yaml_path = CWD / 'fragments' / set_key[5:] / 'input' / 'domains.yaml'
        else:
            yaml_path = CWD / set_key / 'input' / 'domains.yaml'
        restrained_resids = set()
        if yaml_path.exists():
            with open(yaml_path) as f:
                ydata = _yaml.safe_load(f)
            for ranges in ydata.values():
                for ds, de in ranges:
                    for r in range(ds, de + 1):
                        restrained_resids.add(r)
        else:
            # Fallback: use full domain extents (will produce no linkers
            # for chains where extents cover everything)
            for dn, (ds, de) in construct_ranges.items():
                for r in range(ds, de + 1):
                    restrained_resids.add(r)
        linker_resids = sorted(r for r in range(1, n_total + 1)
                                if r not in restrained_resids)
        if not linker_resids:
            return None  # nothing unrestrained — skip

        # Sub-sample linker residues at user-specified stride
        linker_idx = linker_resids[::max(1, ca_stride)]
        names = [f'linker_r{r}_{axis}'
                 for r in linker_idx for axis in 'xyz']
        if add_rg:
            names += [f'rg_{d}' for d in domain_list]
        feature_names = names

        for frame_idx, ts in enumerate(u.trajectory[SKIP_FRAMES:]):
            row = []
            pos = all_ca.positions / 10.0  # nm
            chain_com = pos.mean(axis=0)
            for r in linker_idx:
                p = pos[r - 1] - chain_com  # relative to chain COM
                row.extend([float(p[0]), float(p[1]), float(p[2])])
            if add_rg:
                for d in domain_list:
                    p = ags[d].positions / 10.0
                    com = p.mean(axis=0)
                    rg = float(np.sqrt(np.mean(
                        np.sum((p - com) ** 2, axis=1))))
                    row.append(rg)
            rows.append(row)
            metadata.append((seed, sample, frame_idx + SKIP_FRAMES))

    elif features_mode == 'orient':
        # Inter-domain relative orientations.
        # For each domain, compute principal axes via SVD on centered
        # CA positions. For each pair of domains, the relative rotation
        # matrix R = R_b @ R_a^T (9 elements per pair). Captures hinge
        # motions between rigid domains.
        names = []
        for a, b in pairs:
            for i in range(3):
                for j in range(3):
                    names.append(f'rot_{a}_{b}_{i}{j}')
        if add_rg:
            names += [f'rg_{d}' for d in domain_list]
        feature_names = names

        for frame_idx, ts in enumerate(u.trajectory[SKIP_FRAMES:]):
            row = []
            # Per-domain principal axes
            R_per = {}
            for d in domain_list:
                p = ags[d].positions / 10.0  # nm
                centered = p - p.mean(axis=0)
                # SVD: U * S * Vt where Vt rows are principal axes
                _, _, Vt = np.linalg.svd(centered, full_matrices=False)
                # Ensure right-handed: det(Vt) = +1
                if np.linalg.det(Vt) < 0:
                    Vt[2] = -Vt[2]
                R_per[d] = Vt  # (3, 3) rotation matrix domain → lab
            for a, b in pairs:
                R = R_per[b] @ R_per[a].T  # rotation a → b
                row.extend(R.flatten().tolist())
            if add_rg:
                for d in domain_list:
                    p = ags[d].positions / 10.0
                    com = p.mean(axis=0)
                    rg = float(np.sqrt(np.mean(
                        np.sum((p - com) ** 2, axis=1))))
                    row.append(rg)
            rows.append(row)
            metadata.append((seed, sample, frame_idx + SKIP_FRAMES))

    elif features_mode == 'pose':
        # Robust rigid-body RELATIVE POSE per domain pair -- captures BOTH where
        # B sits relative to A AND how B is rotated relative to A (i.e. what you
        # see by eye). Distances (ca/com) are rotation-invariant and miss this;
        # the old 'orient' used sign-ambiguous SVD principal axes (hence noisy).
        # Here each domain's body frame comes from a Kabsch fit to a SHARED
        # reference (seed-1_sample-0), so frames are sign-consistent across all
        # replicates. Per pair: relative position (3) + relative rotation matrix
        # (9, continuous, no wrap/sign pitfalls) = 12 features.
        ref_pdb = resolve_topology_pdb(resolve_sim_dir(set_key, 1, 0))
        uref = mda.Universe(str(ref_pdb))
        ref_ca = {}
        for dn, (ds, de) in construct_ranges.items():
            ag = uref.select_atoms(f'resid {ds}:{de}')
            if len(ag) > 0:
                ref_ca[dn] = ag.positions.copy() / 10.0
        def _kabsch(P, Q):  # rotation R mapping centered P onto centered Q
            H = P.T @ Q
            U, S, Vt = np.linalg.svd(H)
            d = np.sign(np.linalg.det(Vt.T @ U.T))
            return Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
        names = []
        for a, b in pairs:
            names += [f'pos_{a}_{b}_{c}' for c in 'xyz']
            names += [f'rot_{a}_{b}_{i}{j}' for i in range(3) for j in range(3)]
        feature_names = names
        for frame_idx, ts in enumerate(u.trajectory[SKIP_FRAMES:]):
            row = []
            Rd = {}; com = {}
            for d in domain_list:
                P = ags[d].positions / 10.0
                Q = ref_ca[d]
                R = _kabsch(P - P.mean(0), Q - Q.mean(0))  # current frame -> reference
                Rd[d] = R; com[d] = P.mean(0)
            for a, b in pairs:
                v = Rd[a] @ (com[b] - com[a])     # B's COM in A's body frame (nm)
                Rrel = Rd[a] @ Rd[b].T            # relative orientation
                row += v.tolist() + Rrel.flatten().tolist()
            rows.append(row)
            metadata.append((seed, sample, frame_idx + SKIP_FRAMES))

    elif features_mode == 'inter_exposed':
        # Inter-domain CA-CA distances between EXPOSED residues only,
        # where exposure is determined from frame 1.
        # CG SASA proxy: count CA neighbors within EXPOSURE_CUTOFF nm
        # (ignoring i±2 sequence neighbors). Residues in the bottom
        # EXPOSURE_FRACTION of neighbor counts are considered exposed.
        EXPOSURE_CUTOFF_NM = 1.2  # CA-CA neighbor cutoff
        EXPOSURE_FRACTION = 0.30  # bottom 30% neighbor count = exposed

        all_ca = u.select_atoms('name CA')
        n_total = len(all_ca)

        # Frame 1 (index 0) — determine exposure
        u.trajectory[0]
        pos0 = all_ca.positions / 10.0  # nm
        diff0 = pos0[:, None, :] - pos0[None, :, :]
        d_all0 = np.sqrt(np.sum(diff0 ** 2, axis=2))
        seq_mask = (
            np.abs(np.arange(n_total)[:, None] -
                   np.arange(n_total)[None, :]) > 2)
        n_neighbors = np.sum((d_all0 < EXPOSURE_CUTOFF_NM) & seq_mask,
                              axis=1)
        threshold = np.percentile(n_neighbors, EXPOSURE_FRACTION * 100)
        exposed_global = n_neighbors <= threshold  # bool array, len n_total

        # Map each domain's residues → indices in all_ca, then filter
        # to exposed and subsample by ca_stride
        exposed_per_domain = {}
        for d in domain_list:
            ag = ags[d]
            domain_resids = sorted(set(int(r) for r in ag.residues.resids))
            # CALVADOS top.pdb numbers residues 1..N matching all_ca order
            idx_in_all = [r - 1 for r in domain_resids
                          if 0 <= r - 1 < n_total]
            exposed_idx = [i for i in idx_in_all if exposed_global[i]]
            # Subsample to keep feature count manageable
            exposed_per_domain[d] = exposed_idx[::max(1, ca_stride)]

        # Build feature list — only inter-domain pairs
        pair_indices = []
        names = []
        for a, b in pairs:
            for ia in exposed_per_domain[a]:
                for ib in exposed_per_domain[b]:
                    pair_indices.append((ia, ib))
                    names.append(f'exposed_{a}_r{ia+1}__{b}_r{ib+1}')

        if add_rg:
            names += [f'rg_{d}' for d in domain_list]
        feature_names = names

        if not pair_indices:
            return None  # no exposed cross-domain pairs

        print(f"    [{seed}/{sample}] exposed residues per domain: "
              + ", ".join(f"{d}={len(exposed_per_domain[d])}"
                          for d in domain_list)
              + f" → {len(pair_indices)} cross-domain pairs", flush=True)

        # Rewind to start frame, then iterate from SKIP_FRAMES
        for frame_idx, ts in enumerate(u.trajectory[SKIP_FRAMES:]):
            row = []
            pos = all_ca.positions / 10.0
            for ia, ib in pair_indices:
                row.append(float(np.linalg.norm(pos[ia] - pos[ib])))
            if add_rg:
                for d in domain_list:
                    p = ags[d].positions / 10.0
                    com = p.mean(axis=0)
                    rg = float(np.sqrt(np.mean(
                        np.sum((p - com) ** 2, axis=1))))
                    row.append(rg)
            rows.append(row)
            metadata.append((seed, sample, frame_idx + SKIP_FRAMES))

    elif features_mode == 'interface_ca':
        # CA-CA distances between BOUNDARY residues of adjacent domains
        # (first/last ~10 residues of each domain to its neighbors).
        # These are at the hinges and move slowly with inter-domain
        # rearrangement.
        names = []
        n_edge = max(1, ca_stride // 5)  # use ca-stride to control edge size
        # For each domain, take its first n_edge and last n_edge CAs
        edges = {}
        for d in domain_list:
            ag = ags[d]
            n = len(ag)
            edges[d] = (list(range(min(n_edge, n))),
                         list(range(max(0, n - n_edge), n)))

        for a, b in pairs:
            _, edge_a = edges[a]  # end residues of a
            edge_b, _ = edges[b]  # start residues of b
            for ia in edge_a:
                for ib in edge_b:
                    names.append(f'iface_{a}_{ia}__{b}_{ib}')
        if add_rg:
            names += [f'rg_{d}' for d in domain_list]
        feature_names = names

        for frame_idx, ts in enumerate(u.trajectory[SKIP_FRAMES:]):
            row = []
            pos_cache = {d: ags[d].positions / 10.0 for d in domain_list}
            for a, b in pairs:
                _, edge_a = edges[a]
                edge_b, _ = edges[b]
                pa = pos_cache[a][edge_a]
                pb = pos_cache[b][edge_b]
                diff = pa[:, None, :] - pb[None, :, :]
                d_ij = np.sqrt(np.sum(diff ** 2, axis=2))
                row.extend(d_ij.ravel().tolist())
            if add_rg:
                for d in domain_list:
                    p = pos_cache[d]
                    com = p.mean(axis=0)
                    rg = float(np.sqrt(np.mean(
                        np.sum((p - com) ** 2, axis=1))))
                    row.append(rg)
            rows.append(row)
            metadata.append((seed, sample, frame_idx + SKIP_FRAMES))

    elif features_mode == 'pose_iface':
        # COMBINED rigid-body POSE (relative position + orientation per pair) AND
        # INTERFACE contact distances (domain-edge CA-CA). Captures global
        # arrangement *and* the contact register in one vector. TICA is invariant
        # to per-feature linear scaling (generalized eigenproblem with C(0)), so
        # the two blocks are concatenated directly. Same feature definition for any
        # construct that contains the requested domains -> comparable components
        # across construct sizes.
        ref_pdb = resolve_topology_pdb(resolve_sim_dir(set_key, 1, 0))
        uref = mda.Universe(str(ref_pdb))
        ref_ca = {}
        for dn, (ds, de) in construct_ranges.items():
            ag = uref.select_atoms(f'resid {ds}:{de}')
            if len(ag) > 0:
                ref_ca[dn] = ag.positions.copy() / 10.0
        def _kabsch(P, Q):
            H = P.T @ Q
            U, S, Vt = np.linalg.svd(H)
            d = np.sign(np.linalg.det(Vt.T @ U.T))
            return Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
        n_edge = max(1, ca_stride // 5)
        edges = {}
        for d in domain_list:
            n = len(ags[d])
            edges[d] = (list(range(min(n_edge, n))), list(range(max(0, n - n_edge), n)))
        names = []
        for a, b in pairs:                       # pose block
            names += [f'pos_{a}_{b}_{c}' for c in 'xyz']
            names += [f'rot_{a}_{b}_{i}{j}' for i in range(3) for j in range(3)]
        for a, b in pairs:                       # interface block
            _, edge_a = edges[a]; edge_b, _ = edges[b]
            for ia in edge_a:
                for ib in edge_b:
                    names.append(f'iface_{a}_{ia}__{b}_{ib}')
        feature_names = names

        for frame_idx, ts in enumerate(u.trajectory[SKIP_FRAMES:]):
            pos_cache = {d: ags[d].positions / 10.0 for d in domain_list}
            Rd = {}; com = {}
            for d in domain_list:
                P = pos_cache[d]; Q = ref_ca[d]
                Rd[d] = _kabsch(P - P.mean(0), Q - Q.mean(0)); com[d] = P.mean(0)
            row = []
            for a, b in pairs:                   # pose
                v = Rd[a] @ (com[b] - com[a])
                Rrel = Rd[a] @ Rd[b].T
                row += v.tolist() + Rrel.flatten().tolist()
            for a, b in pairs:                   # interface
                _, edge_a = edges[a]; edge_b, _ = edges[b]
                pa = pos_cache[a][edge_a]; pb = pos_cache[b][edge_b]
                diff = pa[:, None, :] - pb[None, :, :]
                row += np.sqrt(np.sum(diff ** 2, axis=2)).ravel().tolist()
            rows.append(row)
            metadata.append((seed, sample, frame_idx + SKIP_FRAMES))

    elif features_mode == 'segments':
        # Divide chain into consecutive ca_stride-length segments,
        # compute COM of each, then pairwise COM-COM distances.
        # Uses the WHOLE chain (not just domains) — captures global
        # arrangement at coarse resolution.
        all_ca = u.select_atoms('name CA')
        n_total = len(all_ca)
        seg_size = ca_stride  # 50 by default
        n_segs = max(2, n_total // seg_size)
        # Segment residue ranges (0-indexed inclusive start, exclusive end)
        seg_ranges = []
        for i in range(n_segs):
            s = i * seg_size
            e = (i + 1) * seg_size if i < n_segs - 1 else n_total
            seg_ranges.append((s, e))

        seg_pairs = [(i, j) for i in range(n_segs)
                     for j in range(i + 1, n_segs)]
        names = [f'seg{i+1}_seg{j+1}' for i, j in seg_pairs]
        if add_rg:
            names += [f'rg_seg{i+1}' for i in range(n_segs)]
        feature_names = names

        for frame_idx, ts in enumerate(u.trajectory[SKIP_FRAMES:]):
            row = []
            pos = all_ca.positions / 10.0  # nm
            coms = np.array([pos[s:e].mean(axis=0)
                              for s, e in seg_ranges])
            for i, j in seg_pairs:
                row.append(float(np.linalg.norm(coms[i] - coms[j])))
            if add_rg:
                for s, e in seg_ranges:
                    p = pos[s:e]
                    com = p.mean(axis=0)
                    rg = float(np.sqrt(np.mean(
                        np.sum((p - com) ** 2, axis=1))))
                    row.append(rg)
            rows.append(row)
            metadata.append((seed, sample, frame_idx + SKIP_FRAMES))

    elif features_mode == 'site_orient':
        # COM-COM distances (translation, rotation-blind) + a rotationally-aware
        # complement: each active-site domain's "site vector" = unit(pocket COM -
        # domain COM). The cosine angle between that vector and the direction to
        # another domain's COM is invariant to the whole complex's overall
        # tumbling (a global rotation rotates both vectors together, leaving their
        # relative angle unchanged) but changes if the site domain reorients
        # in place -- i.e. captures "is this active site facing toward or away
        # from its neighbor", far cheaper than a full 3x3 relative-rotation
        # matrix per domain pair (cf. 'pose').
        active_sites = reg.get_active_sites_for_set(set_key)
        site_domains = [d for d in domain_list
                        if DOMAIN_TO_UNIT.get(d) in reg.UNIT_TO_SITE
                        and reg.UNIT_TO_SITE[DOMAIN_TO_UNIT[d]] in active_sites]
        pocket_ags = {}
        for d in site_domains:
            sname = reg.UNIT_TO_SITE[DOMAIN_TO_UNIT[d]]
            pocket_resids = active_sites[sname]['pocket']
            sel = ' or '.join(f'resid {r}' for r in pocket_resids)
            ag = u.select_atoms(sel)
            if len(ag) > 0:
                pocket_ags[d] = ag
        site_domains = [d for d in site_domains if d in pocket_ags]

        names = [f'dist_{a}_{b}' for a, b in pairs]
        for d in site_domains:
            for other in domain_list:
                if other != d:
                    names.append(f'face_{d}_to_{other}')
        for i in range(len(site_domains)):
            for j in range(i + 1, len(site_domains)):
                names.append(f'facedot_{site_domains[i]}_{site_domains[j]}')
        if add_rg:
            names += [f'rg_{d}' for d in domain_list]
        feature_names = names

        for frame_idx, ts in enumerate(u.trajectory[SKIP_FRAMES:]):
            row = []
            coms = {d: ags[d].center_of_geometry() / 10.0 for d in domain_list}
            for a, b in pairs:
                row.append(float(np.linalg.norm(coms[a] - coms[b])))
            site_vecs = {}
            for d in site_domains:
                pocket_com = pocket_ags[d].center_of_geometry() / 10.0
                v = pocket_com - coms[d]
                n = np.linalg.norm(v)
                site_vecs[d] = v / n if n > 1e-9 else v
            for d in site_domains:
                for other in domain_list:
                    if other == d:
                        continue
                    dirv = coms[other] - coms[d]
                    dn = np.linalg.norm(dirv)
                    row.append(float(np.dot(site_vecs[d], dirv / dn)) if dn > 1e-9 else 0.0)
            for i in range(len(site_domains)):
                for j in range(i + 1, len(site_domains)):
                    row.append(float(np.dot(site_vecs[site_domains[i]], site_vecs[site_domains[j]])))
            if add_rg:
                for d in domain_list:
                    p = ags[d].positions / 10.0
                    com = p.mean(axis=0)
                    rg = float(np.sqrt(np.mean(np.sum((p - com) ** 2, axis=1))))
                    row.append(rg)
            rows.append(row)
            metadata.append((seed, sample, frame_idx + SKIP_FRAMES))

    else:  # COM-COM (default)
        feature_names = [f'dist_{a}_{b}' for a, b in pairs]
        if add_rg:
            feature_names += [f'rg_{d}' for d in domain_list]
        for frame_idx, ts in enumerate(u.trajectory[SKIP_FRAMES:]):
            row = []
            for a, b in pairs:
                com_a = ags[a].center_of_geometry() / 10.0
                com_b = ags[b].center_of_geometry() / 10.0
                row.append(float(np.linalg.norm(com_a - com_b)))
            if add_rg:
                for d in domain_list:
                    positions = ags[d].positions / 10.0
                    com = positions.mean(axis=0)
                    rg = float(np.sqrt(np.mean(
                        np.sum((positions - com) ** 2, axis=1))))
                    row.append(rg)
            rows.append(row)
            metadata.append((seed, sample, frame_idx + SKIP_FRAMES))

    return {
        'set_key': set_key,
        'features': np.array(rows),
        'metadata': metadata,
        'feature_names': feature_names,
    }


def collect_features(set_key, domains, add_rg, workers,
                     features_mode='com', ca_stride=10):
    """Run extraction across all replicates and concatenate."""
    jobs = []
    for seed in SEEDS:
        for sample in SAMPLES:
            jobs.append((set_key, seed, sample, tuple(domains), add_rg,
                         features_mode, ca_stride))

    print(f"  Processing {len(jobs)} replicates...")
    all_feats = []
    all_meta = []
    feature_names = None

    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(extract_features_one_replicate, j): j
                   for j in jobs}
        for fut in as_completed(futures):
            res = fut.result()
            if res is None:
                continue
            all_feats.append(res['features'])
            all_meta.extend(res['metadata'])
            feature_names = res['feature_names']
            print(f"    seed-{res['metadata'][0][0]}_sample-"
                  f"{res['metadata'][0][1]}: {len(res['features'])} frames")

    if not all_feats:
        return None

    # Record replicate boundaries (start/end indices in the concatenated array)
    rep_boundaries = []
    pos = 0
    for arr in all_feats:
        rep_boundaries.append((pos, pos + len(arr)))
        pos += len(arr)

    return {
        'features': np.concatenate(all_feats, axis=0),
        'metadata': all_meta,
        'feature_names': feature_names,
        'rep_boundaries': rep_boundaries,
    }


# ============================================================
# TICA (time-lagged Independent Component Analysis)
# ============================================================

def compute_tica(X, rep_boundaries, lag=10, n_components=None):
    """Time-lagged Independent Component Analysis on standardized data.

    Builds the instantaneous covariance C(0) and the symmetrized
    time-lagged covariance C(tau). Solves the generalized eigenvalue
    problem C(tau) @ v = lambda * C(0) @ v. Eigenvectors corresponding
    to the largest eigenvalues (slowest motions) are the TICA
    components.

    Time-lag pairs are NOT taken across replicate boundaries.

    Parameters
    ----------
    X : (n_frames, n_features) array, already mean-centered or standardized
    rep_boundaries : list of (start, end) tuples per replicate
    lag : lag time in frames
    n_components : if None, returns all components

    Returns
    -------
    Y : (n_frames, n_kept) projected data sorted by slowness (slowest first)
    eigvals : eigenvalues (corresponding to autocorrelation at lag tau)
    eigvecs : eigenvector matrix
    """
    from scipy.linalg import eigh

    n_frames, n_features = X.shape

    # Build instantaneous and time-lagged covariances by summing across reps
    c00 = np.zeros((n_features, n_features))
    ctt = np.zeros((n_features, n_features))
    n_pairs = 0

    for start, end in rep_boundaries:
        if end - start <= lag:
            continue
        X0 = X[start:end - lag]
        Xt = X[start + lag:end]
        c00 += X0.T @ X0 + Xt.T @ Xt  # both endpoints
        ctt += X0.T @ Xt
        n_pairs += len(X0)

    if n_pairs == 0:
        raise RuntimeError(f"lag={lag} is too large; no usable pairs")

    c00 /= 2 * n_pairs
    ctt /= n_pairs
    # Symmetrize the time-lagged covariance
    ctt = 0.5 * (ctt + ctt.T)

    # Regularize C(0) slightly to ensure positive definiteness
    c00 += 1e-6 * np.eye(n_features)

    # Generalized eigenvalue problem
    eigvals, eigvecs = eigh(ctt, c00)
    # Sort by descending eigenvalue (slowest motions first)
    order = np.argsort(-eigvals)
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]

    if n_components is not None:
        eigvecs = eigvecs[:, :n_components]
        eigvals = eigvals[:n_components]

    # Project
    Y = X @ eigvecs
    return Y, eigvals, eigvecs


def implied_timescales(eigvals, lag_frames, frame_dt_ns=0.01):
    """Convert TICA eigenvalues to implied timescales (ns).

    t_i = -lag / log(|lambda_i|).
    """
    lag_ns = lag_frames * frame_dt_ns
    ts = []
    for lam in eigvals:
        if abs(lam) >= 1.0 or abs(lam) <= 1e-12:
            ts.append(np.nan)
        else:
            ts.append(-lag_ns / np.log(abs(lam)))
    return np.array(ts)


def its_sweep(X, rep_boundaries, lags, n_components=10, frame_dt_ns=0.01):
    """Run TICA at multiple lag times and return implied timescales matrix.

    Returns
    -------
    ts_matrix : (n_lags, n_components) implied timescales in ns
    eig_matrix : (n_lags, n_components) eigenvalues
    """
    n_features = X.shape[1]
    n_comp = min(n_components, n_features)
    ts_matrix = np.full((len(lags), n_comp), np.nan)
    eig_matrix = np.full((len(lags), n_comp), np.nan)
    for i, lag in enumerate(lags):
        try:
            _, eigvals, _ = compute_tica(X, rep_boundaries, lag=lag,
                                          n_components=n_comp)
            eig_matrix[i, :len(eigvals)] = eigvals
            ts = implied_timescales(eigvals, lag, frame_dt_ns)
            ts_matrix[i, :len(ts)] = ts
        except RuntimeError as e:
            print(f"  lag={lag}: skipped ({e})")
    return ts_matrix, eig_matrix


def vamp2_score(X, rep_boundaries, lag, n_components=10):
    """VAMP-2 score: sum of squared singular values of the Koopman operator
    estimate at lag tau. Higher = the featurization captures more slow
    dynamics. Used to compare featurizations: pick the one with the highest
    VAMP-2 at a chosen lag.

    Reference: Wu & Noé (2020), J. Nonlinear Sci 30, 23-66.
    """
    from scipy.linalg import sqrtm

    n_frames, n_features = X.shape
    c00 = np.zeros((n_features, n_features))
    ctt = np.zeros((n_features, n_features))
    c0t = np.zeros((n_features, n_features))
    n_pairs = 0

    for s, e in rep_boundaries:
        if e - s <= lag:
            continue
        X0 = X[s:e - lag]
        Xt = X[s + lag:e]
        c00 += X0.T @ X0
        ctt += Xt.T @ Xt
        c0t += X0.T @ Xt
        n_pairs += len(X0)

    if n_pairs == 0:
        return np.nan

    c00 /= n_pairs
    ctt /= n_pairs
    c0t /= n_pairs

    # Regularize for numerical stability
    reg = 1e-6 * np.eye(n_features)
    c00 += reg
    ctt += reg

    try:
        c00_inv_sqrt = np.real(sqrtm(np.linalg.inv(c00)))
        ctt_inv_sqrt = np.real(sqrtm(np.linalg.inv(ctt)))
    except np.linalg.LinAlgError:
        return np.nan

    K = c00_inv_sqrt @ c0t @ ctt_inv_sqrt
    sing_vals = np.linalg.svd(K, compute_uv=False)
    # VAMP-2 = sum of top n squared singular values
    n_keep = min(n_components, len(sing_vals))
    return float(np.sum(sing_vals[:n_keep] ** 2))


def vamp_sweep(X, rep_boundaries, lags, n_components=10):
    """Compute VAMP-2 score across multiple lag times."""
    return np.array([vamp2_score(X, rep_boundaries, lag, n_components)
                     for lag in lags])


def vamp2_score_cv(X, rep_boundaries, lag, n_components=10,
                   n_folds=5, seed=42):
    """Cross-validated VAMP-2 score (k-fold over replicates).

    For each fold: hold out 1/k replicates as test, fit TICA on rest, then
    project test data using the train eigenvectors and compute the
    test-set autocorrelation² summed over modes.

    Returns array of shape (n_folds,) with one CV score per fold.
    """
    n_reps = len(rep_boundaries)
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n_reps)
    fold_size = max(1, n_reps // n_folds)

    def build(idx_list):
        parts, bnd = [], []
        pos = 0
        for r in idx_list:
            s, e = rep_boundaries[r]
            parts.append(X[s:e])
            bnd.append((pos, pos + (e - s)))
            pos += (e - s)
        if not parts:
            return None, []
        return np.concatenate(parts, axis=0), bnd

    n_feat = X.shape[1]
    n_comp = min(n_components, n_feat)

    scores = np.full(n_folds, np.nan)
    for f in range(n_folds):
        if f == n_folds - 1:
            test_idx = perm[f * fold_size:]
        else:
            test_idx = perm[f * fold_size:(f + 1) * fold_size]
        train_idx = np.setdiff1d(perm, test_idx)

        X_tr, bnd_tr = build(train_idx)
        X_te, bnd_te = build(test_idx)
        if X_tr is None or X_te is None:
            continue

        try:
            # Fit TICA on train — get right eigenvectors
            _, eigvals_tr, eigvecs_tr = compute_tica(
                X_tr, bnd_tr, lag=lag, n_components=n_comp)
        except (RuntimeError, np.linalg.LinAlgError):
            continue

        # Project TEST data using TRAIN eigenvectors
        X_te_proj = X_te @ eigvecs_tr
        # Center test projections per-replicate to remove offsets
        X_te_proj = X_te_proj - X_te_proj.mean(axis=0)
        # Variance per mode (for normalization)
        var_te = X_te_proj.var(axis=0) + 1e-12

        # Compute test-set autocorrelation² in each TICA mode
        score = 0.0
        for c in range(n_comp):
            xc = X_te_proj[:, c]
            num = 0.0
            den = 0
            for s, e in bnd_te:
                if e - s <= lag:
                    continue
                x0 = xc[s:e - lag]
                xt = xc[s + lag:e]
                num += np.sum(x0 * xt)
                den += len(x0)
            if den > 0:
                rho = num / den / var_te[c]
                score += rho ** 2
        scores[f] = score
    return scores


def vamp_sweep_cv(X, rep_boundaries, lags, n_components=10,
                   n_folds=5, seed=42):
    """Cross-validated VAMP-2 sweep across lags."""
    cv_mat = np.full((len(lags), n_folds), np.nan)
    for i, lag in enumerate(lags):
        cv_mat[i] = vamp2_score_cv(X, rep_boundaries, lag,
                                    n_components=n_components,
                                    n_folds=n_folds, seed=seed)
    with np.errstate(invalid='ignore'):
        med = np.nanmedian(cv_mat, axis=1)
        lo = np.nanpercentile(cv_mat, 2.5, axis=1)
        hi = np.nanpercentile(cv_mat, 97.5, axis=1)
    return med, lo, hi, cv_mat


def vamp_sweep_bootstrap(X, rep_boundaries, lags, n_components=10,
                          n_bootstrap=20, seed=42):
    """Block-bootstrap VAMP-2 across lags for 95% CI."""
    n_reps = len(rep_boundaries)
    rng = np.random.default_rng(seed)
    vamp_all = np.full((len(lags), n_bootstrap), np.nan)

    for b in range(n_bootstrap):
        sampled = rng.choice(n_reps, n_reps, replace=True)
        parts, new_bnd = [], []
        pos = 0
        for r in sampled:
            s, e = rep_boundaries[r]
            parts.append(X[s:e])
            new_bnd.append((pos, pos + (e - s)))
            pos += (e - s)
        X_b = np.concatenate(parts, axis=0)
        for i, lag in enumerate(lags):
            try:
                vamp_all[i, b] = vamp2_score(X_b, new_bnd, lag, n_components)
            except Exception:
                pass

    with np.errstate(invalid='ignore'):
        median = np.nanmedian(vamp_all, axis=1)
        lower = np.nanpercentile(vamp_all, 2.5, axis=1)
        upper = np.nanpercentile(vamp_all, 97.5, axis=1)
    return median, lower, upper


def its_sweep_bootstrap(X, rep_boundaries, lags, n_components=10,
                         n_bootstrap=30, frame_dt_ns=0.01, seed=42):
    """Bootstrap ITS sweep with 95% confidence intervals.

    Resamples replicates with replacement (block bootstrap) — preserves
    temporal correlation within each replicate while sampling uncertainty
    in which replicates we have.

    Returns
    -------
    ts_median : (n_lags, n_components) median timescales (ns)
    ts_lower  : (n_lags, n_components) 2.5%
    ts_upper  : (n_lags, n_components) 97.5%
    """
    n_reps = len(rep_boundaries)
    n_features = X.shape[1]
    n_comp = min(n_components, n_features)
    rng = np.random.default_rng(seed)

    # Store all bootstrap timescales: (n_lags, n_components, n_bootstrap)
    ts_all = np.full((len(lags), n_comp, n_bootstrap), np.nan)

    for b in range(n_bootstrap):
        # Block bootstrap: sample replicates with replacement
        sampled = rng.choice(n_reps, n_reps, replace=True)
        # Build X subset and new boundaries
        parts = []
        new_bnd = []
        pos = 0
        for r in sampled:
            s, e = rep_boundaries[r]
            length = e - s
            parts.append(X[s:e])
            new_bnd.append((pos, pos + length))
            pos += length
        X_b = np.concatenate(parts, axis=0)

        for i, lag in enumerate(lags):
            try:
                _, eigvals, _ = compute_tica(X_b, new_bnd, lag=lag,
                                              n_components=n_comp)
                ts = implied_timescales(eigvals, lag, frame_dt_ns)
                ts_all[i, :len(ts), b] = ts
            except (RuntimeError, np.linalg.LinAlgError):
                pass

    # Aggregate over bootstrap samples
    with np.errstate(invalid='ignore'):
        ts_median = np.nanmedian(ts_all, axis=2)
        ts_lower = np.nanpercentile(ts_all, 2.5, axis=2)
        ts_upper = np.nanpercentile(ts_all, 97.5, axis=2)

    return ts_median, ts_lower, ts_upper


# ============================================================
# Clustering
# ============================================================

def silhouette_sweep(X, k_range, sample_size=5000):
    """Compute silhouette score for each k in k_range."""
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score

    # Subsample for speed
    rng = np.random.default_rng(42)
    if len(X) > sample_size:
        idx = rng.choice(len(X), sample_size, replace=False)
        X_sub = X[idx]
    else:
        X_sub = X

    scores = []
    inertias = []
    for k in k_range:
        km = KMeans(n_clusters=k, n_init=10, random_state=42)
        labels = km.fit_predict(X_sub)
        score = silhouette_score(X_sub, labels) if k > 1 else 0.0
        scores.append(score)
        inertias.append(km.inertia_)
        print(f"    k={k}: silhouette={score:.4f}, inertia={km.inertia_:.1f}")
    return np.array(scores), np.array(inertias)


def cluster_kmeans(X, k, seed=42):
    """K-means clustering. Returns (labels, centroids)."""
    from sklearn.cluster import KMeans
    km = KMeans(n_clusters=k, n_init=20, random_state=seed)
    labels = km.fit_predict(X)
    return labels, km.cluster_centers_, km


def cluster_regspace(X, target_n, max_iter=20, seed=42, use_2d_only=False):
    """Regular-space clustering — uniform spatial coverage.

    Adds frames as new cluster centers if they're further than dmin from
    all existing centers. Binary-searches dmin to hit ~target_n centers.

    Parameters
    ----------
    X : (N, D) feature matrix
    target_n : approximate number of microstates desired
    max_iter : binary search iterations for dmin
    use_2d_only : if True, place centers using only the top 2 dims of X
        (gives EVEN coverage on the 2D IC1-IC2 plane specifically).
        Frames still get assigned in full D-dim space.

    Returns
    -------
    labels : (N,) integer cluster IDs
    centers : (n_centers, D) center coordinates
    """
    rng = np.random.default_rng(seed)

    # Subsample for speed
    n_frames = len(X)
    if n_frames > 50000:
        sub_idx = rng.choice(n_frames, 50000, replace=False)
    else:
        sub_idx = np.arange(n_frames)
    X_sub = X[sub_idx]

    # Choose placement space
    X_place = X_sub[:, :2] if use_2d_only else X_sub

    # Binary search for dmin to get ~target_n centers
    dmin_lo, dmin_hi = 0.0, float(np.max(np.std(X_place, axis=0)) * 5)
    best_centers = None
    best_n = -1
    for it in range(max_iter):
        dmin = 0.5 * (dmin_lo + dmin_hi)
        centers_idx = [0]
        for i in range(1, len(X_place)):
            d = np.min(np.linalg.norm(X_place[centers_idx] - X_place[i],
                                       axis=1))
            if d > dmin:
                centers_idx.append(i)
                if len(centers_idx) > target_n * 3:
                    break
        n_centers = len(centers_idx)
        if abs(n_centers - target_n) < abs(best_n - target_n):
            best_centers = centers_idx
            best_n = n_centers
        if n_centers > target_n:
            dmin_lo = dmin
        else:
            dmin_hi = dmin
        if abs(n_centers - target_n) <= max(2, target_n // 20):
            break

    # Use full-D coordinates of selected center frames as actual centers
    centers_full = X_sub[best_centers]

    # Assign ALL frames (not just subsample) to nearest center in full-D
    labels = np.zeros(n_frames, dtype=int)
    chunk_size = 10000
    for s in range(0, n_frames, chunk_size):
        e = min(s + chunk_size, n_frames)
        # (chunk, n_centers) distances
        d_chunk = np.linalg.norm(
            X[s:e, None, :] - centers_full[None, :, :], axis=2)
        labels[s:e] = np.argmin(d_chunk, axis=1)

    return labels, centers_full, None


def find_centroid_frames(X, labels, k):
    """For each cluster, find the frame closest to the cluster centroid.

    Returns dict {cluster_id: frame_index_in_X}.
    """
    centroids = []
    for c in range(k):
        mask = labels == c
        if mask.sum() == 0:
            centroids.append(None)
            continue
        cluster_pts = X[mask]
        cluster_mean = cluster_pts.mean(axis=0)
        # Find the frame closest to the mean
        dists = np.linalg.norm(cluster_pts - cluster_mean, axis=1)
        local_idx = np.argmin(dists)
        # Map back to global index
        global_idx = np.where(mask)[0][local_idx]
        centroids.append(int(global_idx))
    return centroids


def find_spread_frames(X2d, k, min_pct=20.0):
    """Pick k frames EVENLY SPREAD across the 2D landscape via farthest-point
    sampling (max-min), restricted to reasonably populated regions so we don't
    select sparse outliers. Use instead of density-centroid selection when the
    macrostates overlap / nest (e.g. a tiny PCCA+ state inside a big basin).

    Returns (frame_indices[k], labels) where labels assigns every frame to its
    nearest selected representative (Voronoi tiling -> non-overlapping states).
    """
    # 2D density floor: keep frames in cells above the min_pct percentile of
    # occupied-cell density, so representatives sit in populated landscape.
    H, xe, ye = np.histogram2d(X2d[:, 0], X2d[:, 1], bins=50)
    xi = np.clip(np.digitize(X2d[:, 0], xe) - 1, 0, H.shape[0] - 1)
    yi = np.clip(np.digitize(X2d[:, 1], ye) - 1, 0, H.shape[1] - 1)
    dens = H[xi, yi]
    occ = dens[dens > 0]
    floor = np.percentile(occ, min_pct) if occ.size else 0
    cand = np.where(dens >= floor)[0]
    if cand.size < k:
        cand = np.arange(len(X2d))
    P = X2d[cand]
    # seed FPS at the densest candidate (most representative of the main basin)
    start = int(cand[np.argmax(dens[cand])])
    chosen = [start]
    dmin = np.linalg.norm(P - X2d[start], axis=1)
    while len(chosen) < k:
        nxt = int(cand[np.argmax(dmin)])
        chosen.append(nxt)
        dmin = np.minimum(dmin, np.linalg.norm(P - X2d[nxt], axis=1))
    # Voronoi reassignment of ALL frames to nearest representative
    reps = X2d[chosen]
    d_all = np.linalg.norm(X2d[:, None, :] - reps[None, :, :], axis=2)
    labels = np.argmin(d_all, axis=1)
    return chosen, labels


# ============================================================
# Extract representative frame as PDB
# ============================================================

def extract_pdb_for_frame(set_key, seed, sample, frame_idx, out_pdb):
    """Extract a single frame as PDB from a replicate trajectory."""
    import MDAnalysis as mda

    sim_dir = resolve_sim_dir(set_key, seed, sample)
    pdb = resolve_topology_pdb(sim_dir)
    dcds = list(sim_dir.glob('*.dcd'))
    if not pdb.is_file() or not dcds:
        return False
    u = mda.Universe(str(pdb), str(dcds[0]))
    if frame_idx >= len(u.trajectory):
        return False
    u.trajectory[frame_idx]
    with mda.Writer(str(out_pdb), n_atoms=u.atoms.n_atoms) as W:
        W.write(u.atoms)
    return True


# ============================================================
# Plots
# ============================================================

def plot_feature_autocorrelation(X, rep_boundaries, feature_names,
                                  max_lag, frame_dt_ns, outpath, units='ps',
                                  n_show=20):
    """Per-feature time-autocorrelation function (ACF) for each feature.

    For each feature, compute ACF(tau) = <(x_t - mean)(x_{t+tau} - mean)> /
    var(x), averaged across replicates. Plot all features on one axis to
    see which ones have slow dynamics (slow decay) vs fast/noisy (rapid
    decay).

    A feature with no slow signal will have ACF that drops to ~0 within
    a few frames. A feature with metastable dynamics will plateau or
    decay slowly.
    """
    n_frames, n_features = X.shape
    lags = np.arange(0, max_lag + 1)
    acf_matrix = np.zeros((n_features, len(lags)))

    # Center features
    X_centered = X - X.mean(axis=0)
    var = X.var(axis=0)
    var[var < 1e-12] = 1e-12

    for k, lag in enumerate(lags):
        # Sum across replicates, weighted by usable pairs
        num = np.zeros(n_features)
        den = 0
        for s, e in rep_boundaries:
            length = e - s - lag
            if length <= 0:
                continue
            x0 = X_centered[s:e - lag]
            xt = X_centered[s + lag:e]
            num += np.sum(x0 * xt, axis=0)
            den += length
        if den > 0:
            acf_matrix[:, k] = num / den / var

    # Plot
    if units == 'ps':
        lag_axis = lags * frame_dt_ns * 1000  # ps
        unit_str = 'ps'
    else:
        lag_axis = lags * frame_dt_ns
        unit_str = 'ns'

    fig, ax = plt.subplots(figsize=(11, 6))
    colors = plt.cm.tab20(np.linspace(0, 1, max(n_features, 1)))
    # Show only n_show features to keep plot readable. For CA mode with
    # 2000+ features, pick the slowest (most informative): rank by mean
    # ACF over the second half of the lag range.
    show = min(n_features, n_show)
    if n_features > n_show:
        # Rank by mean ACF over second half (slow features have higher mean)
        half = max_lag // 2
        score = acf_matrix[:, half:max_lag + 1].mean(axis=1)
        order = np.argsort(-score)
        feature_ranks = order[:show]
    else:
        feature_ranks = list(range(show))
    for k, i in enumerate(feature_ranks):
        ax.plot(lag_axis, acf_matrix[i], color=colors[k % len(colors)],
                linewidth=1.5, alpha=0.8,
                label=feature_names[i][:30])

    ax.axhline(0, color='black', linewidth=0.5, alpha=0.5)
    ax.axhline(0.1, color='red', linewidth=0.8, linestyle='--', alpha=0.6,
               label='ACF=0.1 (~3τ)')

    ax.set_xlabel(f'Lag / {unit_str}', fontsize=12)
    ax.set_ylabel('Autocorrelation C(τ)', fontsize=12)
    ax.set_title('Per-Feature Time Autocorrelation\n'
                 'Slow features stay high; noise drops to 0 fast',
                 fontsize=13)
    ax.legend(fontsize=6, loc='upper right', ncol=2, framealpha=0.9)
    ax.grid(alpha=0.3)
    ax.set_ylim(-0.2, 1.05)

    plt.tight_layout()
    fig.savefig(f'{outpath}.png', dpi=200, bbox_inches='tight')
    fig.savefig(f'{outpath}.svg', bbox_inches='tight')
    plt.close()

    # Also save the matrix
    return acf_matrix


def plot_its(lags, ts_matrix, frame_dt_ns, outpath, n_show=5, units='ps',
             ts_lower=None, ts_upper=None, x_scale='linear', title=None):
    """PyEMMA-style MSM Implied Timescales (ITS) validation plot.

    Parameters
    ----------
    lags : list of int (lag times in frames)
    ts_matrix : (n_lags, n_components) implied timescales in ns (median or
        point estimate). Plotted as the central line.
    ts_lower, ts_upper : optional 95% CI bands (same shape as ts_matrix).
    frame_dt_ns : ns per frame (CALVADOS at wfreq=1000 → 0.01)
    units : 'ps' (PyEMMA convention) or 'ns'
    x_scale : 'linear' (start at 0) or 'log'
    """
    fig, ax = plt.subplots(figsize=(9, 6))

    # Convert to requested units
    if units == 'ps':
        scale = 1000.0
        unit_str = 'ps'
    else:
        scale = 1.0
        unit_str = 'ns'

    lags_arr = np.array(lags)
    lag_axis = lags_arr * frame_dt_ns * scale
    ts_scaled = ts_matrix * scale
    if ts_lower is not None:
        ts_lower_scaled = ts_lower * scale
        ts_upper_scaled = ts_upper * scale

    colors = plt.cm.viridis(np.linspace(0, 0.85, n_show))
    n_modes = min(n_show, ts_matrix.shape[1])

    for m in range(n_modes):
        valid = ~np.isnan(ts_scaled[:, m])
        if valid.any():
            # CI band
            if ts_lower is not None:
                ax.fill_between(lag_axis[valid],
                                ts_lower_scaled[valid, m],
                                ts_upper_scaled[valid, m],
                                color=colors[m], alpha=0.22)
            ax.plot(lag_axis[valid], ts_scaled[valid, m], 'o-',
                    color=colors[m], linewidth=2, markersize=8,
                    label=f'IC{m+1}')

    # Trust region: timescale must be ≥ lag to be resolved
    if x_scale == 'linear':
        x_extend = np.linspace(0, lag_axis.max() * 1.05, 50)
        ax.fill_between(x_extend, 0, x_extend, color='gray', alpha=0.15,
                        label='Below lag (unresolved)')
        ax.plot(x_extend, x_extend, 'k--', linewidth=1, alpha=0.5)
    else:
        ax.fill_between(lag_axis, 0, lag_axis, color='gray', alpha=0.15,
                        label='Below lag (unresolved)')
        ax.plot(lag_axis, lag_axis, 'k--', linewidth=1, alpha=0.5)

    ax.set_xlabel(f'Lag time / {unit_str}', fontsize=12)
    ax.set_ylabel(f'Timescale / {unit_str}', fontsize=12)
    if title is None:
        title = ('MSM Implied Timescales (95% CI shaded)' if ts_lower is not None
                 else 'MSM Implied Timescales')
    ax.set_title(title + '\nLines should plateau when lag is large enough',
                 fontsize=12)

    if x_scale == 'log':
        ax.set_xscale('log')
    else:
        ax.set_xlim(0, lag_axis.max() * 1.05)
        # Explicit tick marks at the actual lag values so each one is
        # labeled (rather than auto-ticks that often skip or duplicate)
        ax.set_xticks(lag_axis)
        # Tick labels: show in ns if values are large, else ps
        if lag_axis.max() >= 1000:  # ≥ 1 ns
            ax.set_xticklabels([f'{v/1000:.1f}' for v in lag_axis])
            ax.set_xlabel('Lag time (ns)', fontsize=12)
        else:
            ax.set_xticklabels([f'{v:.0f}' for v in lag_axis])
    ax.set_yscale('log')

    ax.legend(fontsize=10, loc='best')
    ax.grid(alpha=0.3, which='both')

    plt.tight_layout()
    fig.savefig(f'{outpath}.png', dpi=200, bbox_inches='tight')
    fig.savefig(f'{outpath}.svg', bbox_inches='tight')
    plt.close()


def plot_tica_spectrum(eigvals, timescales, lag_frames, outpath):
    """TICA eigenvalues and implied timescales — diagnostic plot."""
    n = len(eigvals)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

    idx = np.arange(1, n + 1)
    ax1.bar(idx, eigvals, color='#3498db', alpha=0.7,
            edgecolor='black', linewidth=0.5)
    ax1.axhline(0, color='black', linewidth=0.5)
    ax1.set_xlabel('TICA component', fontsize=11)
    ax1.set_ylabel('Eigenvalue (autocorrelation at lag)', fontsize=11)
    ax1.set_title(f'TICA spectrum (lag={lag_frames} frames = '
                  f'{lag_frames * 0.01:.2f} ns)', fontsize=12)
    ax1.set_xticks(idx)
    ax1.grid(alpha=0.3, axis='y')

    # Implied timescales
    valid = ~np.isnan(timescales)
    ax2.bar(idx[valid], timescales[valid], color='#9b59b6', alpha=0.7,
            edgecolor='black', linewidth=0.5)
    ax2.axhline(lag_frames * 0.01, color='red', linestyle='--', alpha=0.6,
                label=f'Lag ({lag_frames * 0.01:.2f} ns)')
    ax2.set_xlabel('TICA component', fontsize=11)
    ax2.set_ylabel('Implied timescale (ns)', fontsize=11)
    ax2.set_title('Implied timescales (slowest motions)', fontsize=12)
    ax2.set_xticks(idx)
    ax2.legend(fontsize=9)
    ax2.grid(alpha=0.3, axis='y')

    plt.tight_layout()
    fig.savefig(f'{outpath}.png', dpi=200, bbox_inches='tight')
    plt.close()


def plot_silhouette(k_range, sil, inertias, outpath):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
    ax1.plot(k_range, sil, 'o-', color='#9b59b6', linewidth=2,
             markersize=8)
    ax1.set_xlabel('Number of clusters k', fontsize=11)
    ax1.set_ylabel('Silhouette score', fontsize=11)
    ax1.set_title('Cluster validity (higher = better)', fontsize=12)
    ax1.set_xticks(list(k_range))
    ax1.grid(alpha=0.3)

    ax2.plot(k_range, inertias, 's-', color='#e67e22', linewidth=2,
             markersize=8)
    ax2.set_xlabel('Number of clusters k', fontsize=11)
    ax2.set_ylabel('Inertia (within-cluster SSE)', fontsize=11)
    ax2.set_title('Elbow plot (look for kink)', fontsize=12)
    ax2.set_xticks(list(k_range))
    ax2.grid(alpha=0.3)

    plt.tight_layout()
    fig.savefig(f'{outpath}.png', dpi=200, bbox_inches='tight')
    plt.close()


def plot_vamp(lags, vamp_median, vamp_lower, vamp_upper, frame_dt_ns,
              outpath, x_scale='linear', label=None):
    """VAMP-2 score vs lag plot with 95% CI band.

    Higher VAMP-2 = the featurization captures more slow dynamics.
    Used for objective featurization comparison.
    """
    fig, ax = plt.subplots(figsize=(9, 6))
    lag_axis_ns = np.array(lags) * frame_dt_ns

    if vamp_lower is not None and vamp_upper is not None:
        ax.fill_between(lag_axis_ns, vamp_lower, vamp_upper,
                        color='#3498db', alpha=0.22)
    ax.plot(lag_axis_ns, vamp_median, 'o-', color='#2980b9',
            linewidth=2, markersize=9, label=label or 'VAMP-2')

    ax.set_xlabel('Lag time (ns)', fontsize=12)
    ax.set_ylabel('VAMP-2 score', fontsize=12)
    ax.set_title('VAMP-2 featurization quality\n'
                 '(higher = better slow-dynamics capture)', fontsize=12)
    if x_scale == 'linear':
        ax.set_xlim(0, lag_axis_ns.max() * 1.05)
        ax.set_xticks(lag_axis_ns)
        ax.set_xticklabels([f'{v:.1f}' for v in lag_axis_ns])
    else:
        ax.set_xscale('log')
    ax.grid(alpha=0.3)
    ax.legend(fontsize=10)

    plt.tight_layout()
    fig.savefig(f'{outpath}.png', dpi=200, bbox_inches='tight')
    fig.savefig(f'{outpath}.svg', bbox_inches='tight')
    plt.close()


def plot_ck_test(k_arr, T_obs, T_pred, base_lag_frames, frame_dt_ns,
                 outpath, T_obs_lower=None, T_obs_upper=None):
    """Standard Chapman-Kolmogorov test plot.

    Grid of (n_macro x n_macro) subplots. Each shows P(i→j) vs k*tau.
    - Solid black line: T(tau)^k Markov prediction
    - Red dots with CI: observed T(k*tau)

    Good Markov model: dots overlap the line. Mismatch = MSM doesn't
    capture the dynamics at this lag.
    """
    n_macro = T_obs.shape[1]
    fig, axes = plt.subplots(n_macro, n_macro,
                              figsize=(3 * n_macro, 2.5 * n_macro),
                              squeeze=False, sharex=True)

    lag_axis_ns = k_arr * base_lag_frames * frame_dt_ns

    for i in range(n_macro):
        for j in range(n_macro):
            ax = axes[i, j]
            # Markov prediction
            ax.plot(lag_axis_ns, T_pred[:, i, j], '-', color='black',
                    linewidth=2, label='Markov prediction' if (i == 0 and j == 0) else None)
            # Observed with CI
            if T_obs_lower is not None:
                ax.fill_between(lag_axis_ns, T_obs_lower[:, i, j],
                                T_obs_upper[:, i, j],
                                color='#e74c3c', alpha=0.25)
            ax.plot(lag_axis_ns, T_obs[:, i, j], 'o', color='#c0392b',
                    markersize=8, markeredgecolor='black', markeredgewidth=0.5,
                    label='Estimate' if (i == 0 and j == 0) else None)

            ax.set_title(f'{i+1} → {j+1}', fontsize=10)
            ax.set_ylim(-0.02, 1.05)
            if i == n_macro - 1:
                ax.set_xlabel('Lag (ns)', fontsize=9)
            if j == 0:
                ax.set_ylabel('P', fontsize=9)
            ax.grid(alpha=0.3)

    # Shared legend
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right', fontsize=10, ncol=2)

    base_lag_ns = base_lag_frames * frame_dt_ns
    fig.suptitle(f'Chapman-Kolmogorov test '
                 f'(base τ = {base_lag_ns:.2f} ns, '
                 f'{n_macro} macrostates)\n'
                 'Dots = observed, line = Markov prediction',
                 fontsize=13)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(f'{outpath}.png', dpi=200, bbox_inches='tight')
    fig.savefig(f'{outpath}.svg', bbox_inches='tight')
    plt.close()


def plot_microstate_centers(X_2d, centers_2d, populations, outpath,
                              axis_label='IC', method_name='TICA'):
    """Show microstate cluster centers on the 2D landscape, sized by
    population. Background = hexbin density of frames."""
    fig, ax = plt.subplots(figsize=(10, 8))

    # Background density
    hb = ax.hexbin(X_2d[:, 0], X_2d[:, 1], gridsize=80, bins='log',
                    cmap='Greys', mincnt=1, edgecolors='none')
    plt.colorbar(hb, ax=ax, label='log(frame count)', shrink=0.6)

    # Microstate centers, sized by population, colored by population
    pop_norm = populations / populations.max()
    sizes = 20 + 200 * pop_norm
    sc = ax.scatter(centers_2d[:, 0], centers_2d[:, 1], s=sizes,
                     c=populations, cmap='plasma',
                     edgecolors='black', linewidths=0.3, alpha=0.85,
                     zorder=10)
    cbar = plt.colorbar(sc, ax=ax, label='Microstate population (frames)',
                         shrink=0.6)

    ax.set_xlabel(f'{axis_label}1', fontsize=12)
    ax.set_ylabel(f'{axis_label}2', fontsize=12)
    ax.set_title(f'{method_name} 2D landscape with {len(centers_2d)} '
                 f'microstate centers\n'
                 f'(dot size + color = population per microstate)',
                 fontsize=13)

    plt.tight_layout()
    fig.savefig(f'{outpath}.png', dpi=200, bbox_inches='tight')
    fig.savefig(f'{outpath}.svg', bbox_inches='tight')
    plt.close()


def plot_microstate_to_macro(X_2d, micro_labels, macro_labels,
                              centers_2d, macro_assign, outpath,
                              axis_label='IC', method_name='TICA'):
    """Show microstate centers colored by macrostate assignment.

    Demonstrates the PCCA+ coarse-graining: which microstates were
    grouped together into each macrostate.
    """
    from matplotlib.patches import Ellipse

    fig, ax = plt.subplots(figsize=(10, 8))
    n_macro = macro_assign.max() + 1
    colors = plt.cm.tab10(np.linspace(0, 1, n_macro))

    # Background: frames colored by their macrostate (faint)
    for m in range(n_macro):
        mask = macro_labels == m
        if mask.sum() == 0:
            continue
        ax.scatter(X_2d[mask, 0], X_2d[mask, 1], s=1, alpha=0.15,
                   c=[colors[m]], edgecolors='none')

    # Microstate centers, colored by macrostate
    for m in range(n_macro):
        mask = macro_assign == m
        if mask.sum() == 0:
            continue
        ax.scatter(centers_2d[mask, 0], centers_2d[mask, 1],
                   s=90, c=[colors[m]], edgecolors='black',
                   linewidths=0.6, alpha=0.95, zorder=10,
                   label=f'Macrostate {m+1} (n_micro={mask.sum()}, '
                         f'frames={int((macro_labels == m).sum())})')

        # 1.5-sigma ellipse around microstate centers of this macrostate
        if mask.sum() >= 3:
            pts = centers_2d[mask]
            center = pts.mean(axis=0)
            cov = np.cov(pts[:, :2].T)
            try:
                eigvals, eigvecs = np.linalg.eigh(cov)
                angle = np.degrees(np.arctan2(eigvecs[1, 1], eigvecs[0, 1]))
                width, height = 2 * 1.5 * np.sqrt(np.abs(eigvals))
                ax.add_patch(Ellipse(center, width, height, angle=angle,
                                       edgecolor=colors[m], facecolor='none',
                                       linewidth=2.5, alpha=0.85, zorder=11))
            except Exception:
                pass

    ax.set_xlabel(f'{axis_label}1', fontsize=12)
    ax.set_ylabel(f'{axis_label}2', fontsize=12)
    ax.set_title(f'PCCA+ macrostate decomposition\n'
                 f'(dots = microstate centers; same color = same macrostate)',
                 fontsize=13)
    ax.legend(fontsize=9, loc='best')

    plt.tight_layout()
    fig.savefig(f'{outpath}.png', dpi=200, bbox_inches='tight')
    fig.savefig(f'{outpath}.svg', bbox_inches='tight')
    plt.close()


def ck_test(micro_labels, macro_assign, rep_boundaries, base_lag,
            n_lag_multiples=5):
    """Chapman-Kolmogorov test for MSM validation.

    For an MSM at base_lag tau, the prediction is T(k*tau) = T(tau)^k.
    Compute T(tau), T(2*tau), ... T(n*tau) from the actual data and
    compare to T(tau)^k.

    Returns (k_array, T_obs, T_pred) where:
      k_array : [1, 2, ..., n_lag_multiples]
      T_obs[k-1, i, j] = empirical transition probability at lag k*tau
      T_pred[k-1, i, j] = (T_obs[0])^k = Markov prediction
    """
    macro_assign = np.asarray(macro_assign)
    n_macro = macro_assign.max() + 1

    # Map every frame's microstate to its macrostate
    macro_traj = macro_assign[micro_labels]

    T_obs = np.zeros((n_lag_multiples, n_macro, n_macro))
    for k in range(1, n_lag_multiples + 1):
        klag = k * base_lag
        C = np.zeros((n_macro, n_macro))
        for s, e in rep_boundaries:
            if e - s <= klag:
                continue
            a = macro_traj[s:e - klag]
            b = macro_traj[s + klag:e]
            for ai, bj in zip(a, b):
                C[ai, bj] += 1
        # Row-normalize (transition probability)
        row_sums = C.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1
        T_obs[k - 1] = C / row_sums

    # Predicted: T(tau)^k
    T_pred = np.zeros((n_lag_multiples, n_macro, n_macro))
    T1 = T_obs[0]
    T_pred[0] = T1
    for k in range(2, n_lag_multiples + 1):
        T_pred[k - 1] = np.linalg.matrix_power(T1, k)

    return np.arange(1, n_lag_multiples + 1), T_obs, T_pred


def ck_test_bootstrap(micro_labels, macro_assign, rep_boundaries, base_lag,
                      n_lag_multiples=5, n_bootstrap=20, seed=42):
    """Bootstrap CK test for 95% CI on the observed transition probabilities."""
    n_reps = len(rep_boundaries)
    rng = np.random.default_rng(seed)
    n_macro = max(macro_assign) + 1

    obs_all = np.full((n_lag_multiples, n_macro, n_macro, n_bootstrap), np.nan)

    for b in range(n_bootstrap):
        sampled = rng.choice(n_reps, n_reps, replace=True)
        # Build bootstrap micro-trajectory and new boundaries
        parts = []
        new_bnd = []
        pos = 0
        for r in sampled:
            s, e = rep_boundaries[r]
            parts.append(micro_labels[s:e])
            new_bnd.append((pos, pos + (e - s)))
            pos += (e - s)
        boot_traj = np.concatenate(parts)

        try:
            _, T_b, _ = ck_test(boot_traj, macro_assign, new_bnd,
                                 base_lag, n_lag_multiples)
            obs_all[..., b] = T_b
        except Exception:
            pass

    with np.errstate(invalid='ignore'):
        obs_lower = np.nanpercentile(obs_all, 2.5, axis=-1)
        obs_upper = np.nanpercentile(obs_all, 97.5, axis=-1)
    return obs_lower, obs_upper


def build_msm_transition_matrix(micro_labels, rep_boundaries, lag, n_states):
    """Build column-stochastic MSM transition matrix at lag tau from
    microstate trajectories. Skips pairs that cross replicate boundaries.

    Uses reversible MLE: T = 0.5 * (C + C.T) / row_sum.
    """
    C = np.zeros((n_states, n_states))
    for s, e in rep_boundaries:
        if e - s <= lag:
            continue
        a = micro_labels[s:e - lag]
        b = micro_labels[s + lag:e]
        for i, j in zip(a, b):
            C[i, j] += 1

    # Symmetrize (reversible MLE — equilibrium condition)
    C_sym = 0.5 * (C + C.T)
    row_sums = C_sym.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    T = C_sym / row_sums

    # Stationary distribution (proportional to row sums)
    pi = C_sym.sum(axis=1)
    if pi.sum() > 0:
        pi = pi / pi.sum()
    return T, pi, C


def pcca_plus(T, n_macro):
    """PCCA+ coarse-graining via spectral analysis of the transition matrix.

    1. Compute right eigenvectors of T
    2. Take top n_macro (slowest modes)
    3. Apply a simple sign-based assignment (PCCA-style)

    This is a simplified PCCA+ (true PCCA+ uses sign-structure optimization
    via the Schur decomposition). Sufficient for k_macro=3-10.
    """
    # Right eigendecomposition
    eigvals, eigvecs = np.linalg.eig(T.T)  # left eigenvectors of T
    # Sort by real part descending
    order = np.argsort(-eigvals.real)
    eigvals = eigvals[order].real
    eigvecs = eigvecs[:, order].real

    # First eigenvector is stationary (constant); use next (n_macro - 1)
    # for clustering microstates into n_macro groups
    if n_macro <= 1:
        return np.zeros(T.shape[0], dtype=int), eigvals

    # Use sign-structure of top n_macro-1 right eigenvectors of T
    # (T has same eigenvalues as T.T; right eigvec of T = left eigvec of T.T)
    eigvals_T, right_evecs = np.linalg.eig(T)
    order = np.argsort(-eigvals_T.real)
    right_evecs = right_evecs[:, order].real

    # Skip the stationary eigenvector (column 0), use cols 1..n_macro-1
    # as features for K-means clustering of microstates
    if n_macro == 2:
        assign = (right_evecs[:, 1] > np.median(right_evecs[:, 1])).astype(int)
    else:
        from sklearn.cluster import KMeans
        features = right_evecs[:, 1:n_macro]
        km = KMeans(n_clusters=n_macro, n_init=10, random_state=42)
        assign = km.fit_predict(features)
    return assign, eigvals_T[order].real


def _replot_landscapes_from_cache(set_key, args):
    """Re-make landscape_with_states.png from cached data without
    re-running TICA / MSM / clustering."""
    import glob

    # Find matching cache files (set_key may match multiple feature configs)
    suffix_glob = set_key.replace('frag_', '')
    pattern = str(DATA_PATH / f'landscape_{suffix_glob}*.npz')
    matches = sorted(glob.glob(pattern))
    if not matches:
        print(f"  No cache files match {pattern}")
        print(f"  Run cluster_states.py without --replot-landscapes first")
        return 1

    print(f"Found {len(matches)} cached landscape(s) for {set_key}:")
    for m in matches:
        print(f"  {m}")

    for cache_file in matches:
        d = np.load(cache_file, allow_pickle=True)
        X_pca_2d = d['X_pca_2d']
        labels = d['labels']
        centroid_2d = d['centroid_2d']
        k_val = int(d['k'])
        axis_label = str(d['axis_label'])
        method_name = str(d['method_name'])

        # Reconstruct the run_tag from the cache filename
        # Filename: landscape_<suffix>.npz  where suffix is e.g. 'md_ca25_tica'
        cache_name = Path(cache_file).stem.replace('landscape_', '')
        out_dir = FIG_ROOT / '05_clustering' / _date.today().isoformat() / cache_name
        out_dir.mkdir(parents=True, exist_ok=True)
        out = out_dir / 'landscape_with_states'
        plot_density_with_state_boundaries(
            X_pca_2d, labels, k_val, centroid_2d,
            outpath=out, axis_label=axis_label,
            method_name=method_name)
        print(f"  Saved: {out}.png")
    return 0


def plot_density_with_state_boundaries(X_2d, labels, k, centroid_2d,
                                         outpath, axis_label='IC',
                                         method_name='TICA'):
    """TICA/PCA landscape colored by density with outlines around each
    macrostate's region and a small circle at each state's centroid.

    - Background: log-density hexbin of all frames (viridis)
    - Each macrostate gets a smoothed contour outline drawn around the
      frames assigned to it
    - Small circle dot marks the centroid frame of each state
    """
    from matplotlib.colors import LogNorm
    from scipy.ndimage import gaussian_filter

    fig, ax = plt.subplots(figsize=(10, 8))

    # Density background
    hb = ax.hexbin(X_2d[:, 0], X_2d[:, 1], gridsize=100, bins='log',
                    cmap='viridis', mincnt=1, edgecolors='none')
    cbar = plt.colorbar(hb, ax=ax, label='log(frame count)', shrink=0.6)

    # Build a 2D grid covering the data range
    x_min, x_max = X_2d[:, 0].min(), X_2d[:, 0].max()
    y_min, y_max = X_2d[:, 1].min(), X_2d[:, 1].max()
    pad = 0.05 * max(x_max - x_min, y_max - y_min)
    x_min -= pad; x_max += pad; y_min -= pad; y_max += pad
    grid_n = 100
    xg = np.linspace(x_min, x_max, grid_n)
    yg = np.linspace(y_min, y_max, grid_n)
    XG, YG = np.meshgrid(xg, yg)

    # For each macrostate, build a density grid and draw the contour
    state_colors = plt.cm.tab10(np.linspace(0, 1, k))
    for c in range(k):
        mask = labels == c
        if mask.sum() < 10:
            continue
        pts = X_2d[mask]
        # 2D histogram of this state's points
        H, xedges, yedges = np.histogram2d(
            pts[:, 0], pts[:, 1], bins=[xg, yg])
        H = gaussian_filter(H, sigma=2)
        # Outline at a low-density iso-contour to enclose the state
        if H.max() > 0:
            level = H.max() * 0.10
            ax.contour(
                0.5 * (xedges[:-1] + xedges[1:]),
                0.5 * (yedges[:-1] + yedges[1:]),
                H.T, levels=[level],
                colors=[state_colors[c]], linewidths=2.5, alpha=0.9)
            # Add label at population centroid
            ax.plot(centroid_2d[c, 0], centroid_2d[c, 1],
                    marker='o', markersize=10, color=state_colors[c],
                    markeredgecolor='white', markeredgewidth=1.5,
                    zorder=20)
            # Text label slightly offset
            ax.annotate(
                f'S{c+1}\n({mask.sum()}, {mask.sum()/len(labels)*100:.1f}%)',
                xy=(centroid_2d[c, 0], centroid_2d[c, 1]),
                xytext=(8, 8), textcoords='offset points',
                fontsize=10, fontweight='bold',
                color=state_colors[c],
                bbox=dict(boxstyle='round,pad=0.3',
                          facecolor='white', alpha=0.85,
                          edgecolor=state_colors[c], linewidth=1))

    ax.set_xlabel(f'{axis_label}1', fontsize=12)
    ax.set_ylabel(f'{axis_label}2', fontsize=12)
    ax.set_title(f'{method_name} landscape with {k} state boundaries\n'
                 f'(density background, contours = state regions, '
                 f'dots = centroids)', fontsize=12)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    plt.tight_layout()
    fig.savefig(f'{outpath}.png', dpi=200, bbox_inches='tight')
    fig.savefig(f'{outpath}.svg', bbox_inches='tight')
    plt.close()


def plot_density_with_clusters(X_2d, labels, k, outpath, axis_label='PC',
                                method_name='PCA'):
    """2D density (hexbin) with K-means cluster overlay (contours + centroids).

    For each cluster:
      - Draw a 1.5-sigma ellipse around its center
      - Mark the centroid with a star
    Background = log-density hexbin of all frames.
    """
    from matplotlib.patches import Ellipse

    fig, ax = plt.subplots(figsize=(10, 8))

    # Background hexbin density (log scale)
    hb = ax.hexbin(X_2d[:, 0], X_2d[:, 1], gridsize=80, bins='log',
                    cmap='Greys', mincnt=1, edgecolors='none')
    plt.colorbar(hb, ax=ax, label='log(frame count)', shrink=0.6)

    # Cluster ellipses and centroids
    colors = plt.cm.tab10(np.linspace(0, 1, k))
    for c in range(k):
        mask = labels == c
        if mask.sum() < 5:
            continue
        pts = X_2d[mask]
        center = pts.mean(axis=0)

        # Covariance for 1.5-sigma ellipse
        cov = np.cov(pts[:, :2].T)
        eigvals, eigvecs = np.linalg.eigh(cov)
        angle = np.degrees(np.arctan2(eigvecs[1, 1], eigvecs[0, 1]))
        width, height = 2 * 1.5 * np.sqrt(eigvals)

        ax.add_patch(Ellipse(center, width, height, angle=angle,
                              edgecolor=colors[c], facecolor='none',
                              linewidth=2.5, alpha=0.85, zorder=10))
        ax.plot(*center, marker='*', markersize=22, color=colors[c],
                markeredgecolor='black', markeredgewidth=1.5, zorder=11,
                label=f'State {c+1} (n={mask.sum()}, '
                      f'{mask.sum()/len(labels)*100:.1f}%)')

    ax.set_xlabel(f'{axis_label}1', fontsize=12)
    ax.set_ylabel(f'{axis_label}2', fontsize=12)
    ax.set_title(f'{method_name} 2D density with K-means cluster overlay\n'
                 f'(stars = cluster centers, ellipses = 1.5σ)', fontsize=13)
    ax.legend(fontsize=9, loc='best')

    plt.tight_layout()
    fig.savefig(f'{outpath}.png', dpi=200, bbox_inches='tight')
    fig.savefig(f'{outpath}.svg', bbox_inches='tight')
    plt.close()


def plot_pca_clusters(X_pca, labels, k, centroid_indices, outpath,
                      axis_label='PC'):
    """2D scatter colored by cluster (works for PCA or TICA projections)."""
    fig, ax = plt.subplots(figsize=(10, 8))
    colors = plt.cm.tab10(np.linspace(0, 1, k))

    # Plot each cluster
    for c in range(k):
        mask = labels == c
        ax.scatter(X_pca[mask, 0], X_pca[mask, 1], s=2, alpha=0.3,
                   c=[colors[c]], label=f'State {c+1} (n={mask.sum()})',
                   edgecolors='none')

    # Mark centroid frames
    for c, idx in enumerate(centroid_indices):
        if idx is not None:
            ax.plot(X_pca[idx, 0], X_pca[idx, 1], marker='*',
                    markersize=20, color=colors[c],
                    markeredgecolor='black', markeredgewidth=1.5, zorder=10)

    ax.set_xlabel(f'{axis_label}1', fontsize=12)
    ax.set_ylabel(f'{axis_label}2', fontsize=12)
    method = 'TICA (slow modes)' if axis_label == 'IC' else 'PCA (variance)'
    ax.set_title(f'Clustered conformational states in 2D {method}\n'
                 '(stars = centroid frames per state)', fontsize=13)
    ax.legend(fontsize=10, loc='best')

    plt.tight_layout()
    fig.savefig(f'{outpath}.png', dpi=200, bbox_inches='tight')
    plt.close()


def plot_state_distance_profile(X_raw, labels, feature_names, k, outpath):
    """Mean feature values per state (rows = states, cols = features).

    Annotated heatmap.
    """
    means = np.zeros((k, X_raw.shape[1]))
    stds = np.zeros((k, X_raw.shape[1]))
    for c in range(k):
        mask = labels == c
        if mask.sum() > 0:
            means[c] = X_raw[mask].mean(axis=0)
            stds[c] = X_raw[mask].std(axis=0)

    n_feats = len(feature_names)

    # Cap figure width to avoid matplotlib's 2^16 pixel limit.
    # For >100 features, switch to a fixed-width heatmap with no labels.
    if n_feats > 100:
        # No per-feature tick labels — just imshow as a heatmap strip
        fig_w = 18
        annotate = False
    else:
        # Width scales with feature count, with cap
        fig_w = min(40, max(8, n_feats * 0.8))
        annotate = True

    fig, ax = plt.subplots(figsize=(fig_w, max(4, k * 0.6)))
    im = ax.imshow(means, cmap='RdYlBu_r', aspect='auto')
    if annotate:
        ax.set_xticks(range(n_feats))
        ax.set_xticklabels(feature_names, rotation=45, ha='right', fontsize=9)
    else:
        # Show every Nth label to keep legible
        stride = max(1, n_feats // 30)
        ticks = list(range(0, n_feats, stride))
        ax.set_xticks(ticks)
        ax.set_xticklabels([feature_names[i] for i in ticks],
                            rotation=90, fontsize=6)
    ax.set_yticks(range(k))
    ax.set_yticklabels([f'State {c+1}' for c in range(k)], fontsize=10)
    plt.colorbar(im, ax=ax, label='Mean distance (nm)', shrink=0.8)

    # Annotate only when feature count is small
    if annotate:
        for i in range(k):
            for j in range(n_feats):
                val = means[i, j]
                color = ('white'
                         if val < means.min() + 0.3 * (means.max() - means.min())
                         else 'black')
                ax.text(j, i, f'{val:.2f}\n±{stds[i,j]:.2f}',
                        ha='center', va='center', fontsize=7, color=color)

    ax.set_title(f'Mean feature values per cluster state '
                 f'(n_features={n_feats})', fontsize=13)
    plt.tight_layout()
    fig.savefig(f'{outpath}.png', dpi=200, bbox_inches='tight')
    plt.close()


def plot_population(labels, k, outpath):
    """Bar chart of cluster populations + percentages."""
    counts = np.bincount(labels, minlength=k)
    fracs = counts / counts.sum()

    fig, ax = plt.subplots(figsize=(8, 5))
    colors = plt.cm.tab10(np.linspace(0, 1, k))
    bars = ax.bar(range(1, k + 1), counts, color=colors, edgecolor='black',
                  linewidth=0.5)
    for bar, frac in zip(bars, fracs):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                f'{frac*100:.1f}%', ha='center', va='bottom', fontsize=10)
    ax.set_xlabel('State', fontsize=12)
    ax.set_ylabel('Number of frames', fontsize=12)
    ax.set_title(f'Cluster populations (total {counts.sum()} frames)',
                 fontsize=13)
    ax.set_xticks(range(1, k + 1))
    plt.tight_layout()
    fig.savefig(f'{outpath}.png', dpi=200, bbox_inches='tight')
    plt.close()


def plot_replicate_state_distribution(metadata, labels, k, outpath):
    """Heatmap showing how each replicate's frames distribute across states."""
    rep_state = {}  # (seed, sample) -> [count per state]
    for (seed, sample, _), c in zip(metadata, labels):
        key = (seed, sample)
        if key not in rep_state:
            rep_state[key] = np.zeros(k, dtype=int)
        rep_state[key][c] += 1

    rep_keys = sorted(rep_state.keys())
    mat = np.array([rep_state[k] for k in rep_keys])
    mat_norm = mat / mat.sum(axis=1, keepdims=True)

    fig, ax = plt.subplots(figsize=(8, max(6, len(rep_keys) * 0.25)))
    im = ax.imshow(mat_norm, cmap='YlGnBu', aspect='auto', vmin=0, vmax=1)
    ax.set_xticks(range(k))
    ax.set_xticklabels([f'State {c+1}' for c in range(k)], fontsize=10)
    ax.set_yticks(range(len(rep_keys)))
    ax.set_yticklabels([f'seed-{s}_sample-{p}' for s, p in rep_keys],
                       fontsize=7)
    plt.colorbar(im, ax=ax, label='Fraction of frames', shrink=0.8)

    # Annotate
    for i in range(len(rep_keys)):
        for j in range(k):
            v = mat_norm[i, j]
            if v > 0.05:
                ax.text(j, i, f'{v*100:.0f}', ha='center', va='center',
                        fontsize=6,
                        color='white' if v > 0.5 else 'black')

    ax.set_title('Replicate-by-state population (% of frames)', fontsize=12)
    plt.tight_layout()
    fig.savefig(f'{outpath}.png', dpi=200, bbox_inches='tight')
    plt.close()


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--set', required=False, default=None,
                        help='Set name (e.g. fl_optimized, fl, or fragment '
                             'name like md1l1_md2_md3)')
    parser.add_argument('--sim-folder', default=None, metavar='PATH',
                        help='Cluster an arbitrary simulation folder (dir with '
                             'seed-{1-5}_sample-{0-4}/ subdirs each having '
                             'top.pdb + <sysname>.dcd). Domain composition comes '
                             'from its metadata.json, --units, or a known '
                             'basename. Use instead of --set for new sims.')
    parser.add_argument('--units', nargs='+', default=None, metavar='UNIT',
                        help='FL domain units present in the --sim-folder '
                             'construct (e.g. md1l1 md2 md3). Overrides '
                             'metadata.json. Ignored without --sim-folder.')
    parser.add_argument('--k', type=int, default=5,
                        help='Number of clusters (default 5)')
    parser.add_argument('--rep-mode', choices=['density', 'spread'], default='density',
                        help="How to pick the k representative states. 'density' "
                             "(default): frame nearest each macrostate mean (can "
                             "overlap/nest). 'spread': farthest-point sampling across "
                             "the populated 2D landscape + Voronoi reassignment, so "
                             "states are evenly distributed and non-overlapping.")
    parser.add_argument('--domains', nargs='+',
                        default=['RRM1', 'RRM2', 'RRM3',
                                  'KH1-KH6', 'KH7a',
                                  'MD1L1', 'MD2', 'MD3',
                                  'KHb-KH8', 'WWE', 'ART'],
                        help='Domains for inter-domain features (default '
                             '= all 11 grouped units; auto-filtered to '
                             'those present in this set)')
    parser.add_argument('--add-rg', action='store_true',
                        help='Include per-domain Rg as additional features')
    parser.add_argument('--pca', type=int, default=None,
                        help='Reduce to N components before clustering '
                             '(default for PCA: keep 95%% variance; for TICA: '
                             'min(n_features, 10))')
    parser.add_argument('--reduce', choices=['pca', 'tica'], default='pca',
                        help='Dimensionality reduction: pca (variance) or '
                             'tica (slow conformational motions)')
    parser.add_argument('--tica-lag', type=int, default=10,
                        help='TICA lag time in frames (default 10 = 0.1 ns '
                             'at wfreq=1000 / dt=0.01 ps). Larger = slower '
                             'motions emphasized.')
    parser.add_argument('--its-lags', nargs='+', type=int, default=None,
                        help='List of TICA lag times (frames) to compute an '
                             'implied-timescales (ITS) plot. e.g. '
                             '"--its-lags 5 10 25 50". Default: no ITS plot.')
    parser.add_argument('--features',
                        choices=['com', 'ca', 'torsion', 'segments',
                                  'linker_ca', 'orient', 'interface_ca',
                                  'inter_exposed', 'pose', 'pose_iface',
                                  'site_orient'],
                        default='com',
                        help='Feature type: com (COM-COM domain distances, '
                             'default) or ca (pairwise CA-CA distances '
                             'between strided CA atoms). site_orient: COM-COM '
                             'distances + active-site-facing direction per '
                             'active-site domain (rotationally aware, far '
                             'lower-dimensional than pose).')
    parser.add_argument('--ca-stride', type=int, default=10,
                        help='When --features ca, take every Nth CA per '
                             'domain (default 10)')
    parser.add_argument('--acf-max-lag', type=int, default=500,
                        help='Max lag (frames) for per-feature ACF '
                             'diagnostic plot (default 500 = 5 ns)')
    parser.add_argument('--n-microstates', type=int, default=None,
                        help='Microstate discretization: N cluster centers '
                             '(default: 200 if --reduce tica, else skip). '
                             'Microstates feed into PCCA+ to produce final '
                             'macrostates of size --k.')
    parser.add_argument('--clustering', choices=['kmeans', 'regspace'],
                        default='kmeans',
                        help='Microstate clustering method: kmeans (default, '
                             'density-aware; more centers in dense regions) '
                             'or regspace (uniform spatial coverage; one '
                             'center per IC-space region regardless of '
                             'density).')
    parser.add_argument('--regspace-2d-only', action='store_true',
                        help='With --clustering regspace, place centers '
                             'using only the top 2 ICs (gives EVEN coverage '
                             'on the IC1-IC2 plane). Otherwise uses full '
                             'TICA space for placement.')
    parser.add_argument('--pcca', action='store_true',
                        help='Use PCCA+ to coarse-grain microstates into '
                             '--k macrostates (kinetically meaningful). '
                             'Requires --n-microstates and MSM build.')
    parser.add_argument('--msm-lag', type=int, default=None,
                        help='Lag time (frames) for MSM transition matrix. '
                             'Default: same as --tica-lag.')
    parser.add_argument('--ck-test', action='store_true',
                        help='Run Chapman-Kolmogorov MSM validation test '
                             'after PCCA+. Requires --pcca (or default '
                             'TICA path with microstates).')
    parser.add_argument('--ck-nlags', type=int, default=5,
                        help='Number of lag multiples in CK test '
                             '(plot will show k=1..N). Default 5.')
    parser.add_argument('--ck-bootstrap', type=int, default=20,
                        help='Bootstrap resamples for CK test 95%% CI '
                             '(default 20). Set 0 to skip.')
    parser.add_argument('--replot-landscapes', action='store_true',
                        help='Skip everything else and just re-make the '
                             'landscape_with_states plot from cached data '
                             'at data/landscape_*.npz. Useful for changing '
                             'plot style without recomputing TICA + MSM.')
    parser.add_argument('--screen', action='store_true',
                        help='Screen mode: produce diagnostic plots only '
                             '(ACF + 2D density of PCA1/PC2 and TIC1/TIC2 '
                             'with cluster overlay). Skips silhouette '
                             'sweep, PDB extraction, full plots. Use to '
                             'choose featurization strategy fast.')
    parser.add_argument('--n-acf-features', type=int, default=20,
                        help='Number of features to show on ACF plot '
                             '(default 20). Helps when --features ca '
                             'creates 2000+ features.')
    parser.add_argument('--no-silhouette', action='store_true',
                        help='Skip silhouette k-sweep (faster)')
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--skip-pdb', action='store_true',
                        help='Skip extracting representative PDBs')
    args = parser.parse_args()

    if bool(args.set) == bool(args.sim_folder):
        parser.error('provide exactly one of --set or --sim-folder')

    if args.sim_folder:
        set_key = register_sim_folder(args.sim_folder, units=args.units)
    else:
        set_key = args.set
        # Auto-prefix fragments
        if set_key not in ('fl', 'fl_optimized', 'md', 'core', 'mka',
                           'norrm', 'noart', 'md3art'):
            if (CWD / 'fragments' / set_key).is_dir():
                set_key = f'frag_{set_key}'

    # ── Replot-only short-circuit ──
    if args.replot_landscapes:
        return _replot_landscapes_from_cache(set_key, args)

    # Resolve effective domain list: drop those not in the construct
    available_ranges = get_construct_domain_ranges(set_key, args.domains)
    effective_domains = list(available_ranges.keys())
    dropped = [d for d in args.domains if d not in effective_domains]

    print(f"=" * 60)
    print(f"Clustering {set_key}")
    print(f"=" * 60)
    print(f"  Domains requested: {args.domains}")
    print(f"  Domains in set:    {effective_domains}")
    if dropped:
        print(f"  Domains dropped (not in this construct): {dropped}")
    print(f"  Add Rg:      {args.add_rg}")
    print(f"  K:           {args.k}")

    if len(effective_domains) < 2:
        print(f"\nERROR: need at least 2 domains for pairwise distances, "
              f"got {effective_domains}")
        return 1

    # Use the effective domain list for everything downstream
    args.domains = effective_domains

    # ── Set up dated output directory per featurization+method ──
    if args.features == 'ca':
        feat_tag = f'ca{args.ca_stride}'
    elif args.features == 'linker_ca':
        feat_tag = f'linker{args.ca_stride}'
    elif args.features == 'interface_ca':
        feat_tag = f'iface{args.ca_stride}'
    elif args.features == 'segments':
        feat_tag = f'seg{args.ca_stride}'
    elif args.features == 'orient':
        feat_tag = 'orient'
    elif args.features == 'pose':
        feat_tag = 'pose'
    elif args.features == 'pose_iface':
        feat_tag = 'poseiface'
    elif args.features == 'torsion':
        feat_tag = 'torsion'
    elif args.features == 'inter_exposed':
        feat_tag = f'exposed{args.ca_stride}'
    elif args.features == 'site_orient':
        feat_tag = 'siteorient'
    else:
        feat_tag = 'com'
    if args.add_rg:
        feat_tag = f'{feat_tag}_rg'
    # A --domains subset must never collide (same cache file / output dir) with a
    # different subset run at the same features_mode on the same set -- e.g. all-11
    # domains vs. just the 4 active-site domains both under --features com. Append a
    # short fingerprint whenever the requested domains aren't the full default set.
    domains_suffix = _domains_suffix(effective_domains)
    feat_tag_full = f'{feat_tag}{domains_suffix}'
    run_tag = f'{set_key.replace("frag_","")}_{feat_tag_full}_{args.reduce}'
    import sys as _sys
    _sys.path.insert(0, str(CWD))
    from _fig_layout import get_fig_dir as _get_fig_dir
    global FIG_PATH
    FIG_PATH = _get_fig_dir('05_clustering', subname=run_tag)
    print(f"  Output dir: {FIG_PATH}")

    # ── Step 1: extract features ──
    print("\n[1] Extracting features per frame...")
    # Cache name encodes feature type + stride (so com and ca-stride-10 coexist),
    # plus the domains_suffix computed above (so a --domains subset never loads/
    # clobbers another subset's cache at the same features_mode).
    if args.features == 'com':
        feat_tag_cache = domains_suffix
    else:
        feat_tag_cache = f'_{feat_tag}{domains_suffix}'
    cache_npz = DATA_PATH / (
        f'cluster_features_{set_key.replace("frag_","")}'
        f'{feat_tag_cache}.npz')

    if cache_npz.exists():
        print(f"  Loading cached features from {cache_npz}")
        cached = np.load(cache_npz, allow_pickle=True)
        X_raw = cached['features']
        metadata = list(cached['metadata'])
        feature_names = list(cached['feature_names'])
        if 'rep_boundaries' in cached:
            rep_boundaries = [tuple(b) for b in cached['rep_boundaries']]
        else:
            # Reconstruct from metadata (seed, sample, frame_idx)
            rep_boundaries = []
            cur_key = None
            start = 0
            for i, m in enumerate(metadata):
                key = (int(m[0]), int(m[1]))
                if cur_key is None:
                    cur_key = key
                elif key != cur_key:
                    rep_boundaries.append((start, i))
                    start = i
                    cur_key = key
            rep_boundaries.append((start, len(metadata)))
    else:
        data = collect_features(set_key, args.domains, args.add_rg,
                                args.workers, args.features, args.ca_stride)
        if data is None:
            print("ERROR: no replicates found")
            return 1
        X_raw = data['features']
        metadata = data['metadata']
        feature_names = data['feature_names']
        rep_boundaries = data['rep_boundaries']
        np.savez(cache_npz,
                 features=X_raw,
                 metadata=np.array(metadata, dtype=object),
                 feature_names=np.array(feature_names),
                 rep_boundaries=np.array(rep_boundaries))
        print(f"  Saved cache: {cache_npz}")

    print(f"  Feature matrix: {X_raw.shape}")
    print(f"  Features: {feature_names}")

    # ── Step 2: standardize ──
    from sklearn.preprocessing import StandardScaler
    scaler = StandardScaler()
    X_std = scaler.fit_transform(X_raw)

    # ── Diagnostic: per-feature autocorrelation ──
    # Helps decide if features have slow dynamics worth clustering on.
    suffix_base = set_key.replace('frag_', '')
    if args.features == 'ca':
        suffix_base = f'{suffix_base}_ca{args.ca_stride}'
    print(f"\n[2pre] Computing per-feature autocorrelation (max lag "
          f"{args.acf_max_lag} frames = {args.acf_max_lag * 0.01:.1f} ns)...")
    plot_feature_autocorrelation(
        X_std, rep_boundaries, feature_names,
        max_lag=args.acf_max_lag, frame_dt_ns=0.01,
        outpath=FIG_PATH / 'acf',
        units='ps', n_show=args.n_acf_features)

    # ── Screen mode: PCA + TICA density plots with K-means overlay ──
    if args.screen:
        from sklearn.decomposition import PCA
        from sklearn.cluster import KMeans

        print(f"\n[SCREEN] Generating diagnostic density plots...")
        print(f"  Output suffix: {suffix_base}")

        # PCA → density + clusters
        n_pca = min(args.pca if args.pca else 10, X_std.shape[1])
        pca = PCA(n_components=n_pca, random_state=42)
        X_pca = pca.fit_transform(X_std)
        evr = pca.explained_variance_ratio_
        print(f"  PCA: {n_pca} components, "
              f"PC1={evr[0]*100:.1f}% / PC2={evr[1]*100:.1f}% / "
              f"total95={np.argmax(np.cumsum(evr) >= 0.95) + 1} comps")
        km_pca = KMeans(n_clusters=args.k, n_init=10, random_state=42)
        labels_pca = km_pca.fit_predict(X_pca)
        plot_density_with_clusters(
            X_pca, labels_pca, args.k,
            outpath=FIG_PATH / 'density_pca',
            axis_label='PC', method_name='PCA')

        # TICA → density + clusters (use --tica-lag, plus run ITS sweep)
        try:
            X_tica, eigvals, _ = compute_tica(
                X_std, rep_boundaries, lag=args.tica_lag,
                n_components=n_pca)
            timescales = implied_timescales(eigvals, args.tica_lag)
            print(f"  TICA (lag={args.tica_lag} frames "
                  f"= {args.tica_lag * 0.01:.2f} ns): "
                  f"top timescales {['%.2f ns' % t for t in timescales[:3]]}")
            km_tica = KMeans(n_clusters=args.k, n_init=10, random_state=42)
            labels_tica = km_tica.fit_predict(X_tica)
            plot_density_with_clusters(
                X_tica, labels_tica, args.k,
                outpath=FIG_PATH / 'density_tica',
                axis_label='IC', method_name='TICA')

            # ITS sweep
            if args.its_lags:
                its_lags = sorted(set(args.its_lags))
            else:
                # Lags = 0.2, 0.4, 0.8, 1.4, 2.0, 4.0 ns at 10 ps/frame
                its_lags = [20, 40, 80, 140, 200, 400]
            print(f"  ITS sweep at lags {its_lags} frames...")
            ts_matrix, _ = its_sweep(
                X_std, rep_boundaries, its_lags, n_components=n_pca)
            plot_its(its_lags, ts_matrix, frame_dt_ns=0.01,
                     outpath=FIG_PATH / 'its',
                     units='ps')
        except RuntimeError as e:
            print(f"  TICA failed: {e}")

        print(f"\n[SCREEN] Done. Diagnostic plots in {FIG_PATH}/:")
        print(f"  cluster_acf_{suffix_base}.png")
        print(f"  cluster_density_pca_{suffix_base}.png")
        print(f"  cluster_density_tica_{suffix_base}.png")
        print(f"  cluster_its_{suffix_base}.png")
        print(f"\nTo run full analysis after picking the best featurization:")
        print(f"  $PY cluster_states.py --set {set_key} --k {args.k} "
              f"--reduce {args.reduce}"
              + (f" --features ca --ca-stride {args.ca_stride}"
                  if args.features == 'ca' else "")
              + (f" --add-rg" if args.add_rg else ""))
        return 0

    # ── Step 3: dimensionality reduction (PCA or TICA) ──
    suffix = suffix_base
    reduce_method = args.reduce
    if reduce_method == 'tica':
        from sklearn.decomposition import PCA
        # First reduce to enough PCA components to capture variance (TICA
        # works on standardized data but benefits from a sane dimension)
        n_components = args.pca if args.pca else min(X_std.shape[1], 10)

        # ── ITS plot (multi-lag) ──
        # Default lags if --its-lags not provided: 4 lags spanning the
        # plausible range for 1500-frame trajectories (= 15 ns each)
        if args.its_lags:
            its_lags = sorted(set(args.its_lags))
        else:
            its_lags = [5, 10, 25, 50, 100, 200]
        print(f"\n[2a] MSM ITS validation: TICA across "
              f"{len(its_lags)} lag times {its_lags} frames "
              f"({[f'{l*0.01:.2f}ns' for l in its_lags]})...")
        ts_matrix, eig_matrix = its_sweep(
            X_std, rep_boundaries, its_lags,
            n_components=n_components)
        # Print converged-looking entries
        print(f"  {'Lag(ns)':<10} " + ' '.join(
            f'IC{i+1}(ns)' for i in range(min(5, n_components))))
        for i, lag in enumerate(its_lags):
            row = ts_matrix[i, :min(5, n_components)]
            row_str = ' '.join(f'{v:6.2f}' if not np.isnan(v) else '   nan'
                                for v in row)
            print(f"  {lag * 0.01:<10.2f} {row_str}")
        plot_its(its_lags, ts_matrix, frame_dt_ns=0.01,
                 outpath=FIG_PATH / 'its', units='ps')
        plot_its(its_lags, ts_matrix, frame_dt_ns=0.01,
                 outpath=FIG_PATH / 'its_ns', units='ns')
        # Save raw data
        np.savez(DATA_PATH / f'cluster_its_{suffix}.npz',
                 lags=its_lags,
                 ts_matrix=ts_matrix,
                 eig_matrix=eig_matrix,
                 frame_dt_ns=0.01)

        # Now compute TICA at the chosen --tica-lag for clustering
        print(f"\n[2b] TICA for clustering (lag={args.tica_lag} frames = "
              f"{args.tica_lag * 0.01:.2f} ns)")
        X_red, tica_eigvals, tica_evecs = compute_tica(
            X_std, rep_boundaries, lag=args.tica_lag,
            n_components=n_components)
        timescales = implied_timescales(tica_eigvals, args.tica_lag)
        print(f"  Eigenvalues:        "
              f"{['%.3f' % v for v in tica_eigvals[:5]]}")
        print(f"  Implied timescales: "
              f"{['%.2f ns' % t if not np.isnan(t) else 'n/a' for t in timescales[:5]]}")

        # Save TICA spectrum plot
        plot_tica_spectrum(tica_eigvals, timescales, args.tica_lag,
                            FIG_PATH / 'tica_spectrum')
    else:
        from sklearn.decomposition import PCA
        if args.pca:
            n_pca = min(args.pca, X_std.shape[1])
        else:
            pca_full = PCA().fit(X_std)
            n_pca = int(np.argmax(np.cumsum(pca_full.explained_variance_ratio_)
                                  >= 0.95) + 1)
            n_pca = max(2, n_pca)
        pca = PCA(n_components=n_pca, random_state=42)
        X_red = pca.fit_transform(X_std)
        print(f"\n[2] PCA: {n_pca} components, "
              f"explained variance = "
              f"{pca.explained_variance_ratio_.sum()*100:.1f}%")
        print(f"  Components: "
              f"{[f'{r*100:.1f}%' for r in pca.explained_variance_ratio_]}")

    X_pca = X_red  # downstream code uses this name
    # Tag plot suffix with reduction method
    suffix = f"{suffix}_{reduce_method}"

    # ── Step 4: silhouette sweep ──
    if not args.no_silhouette:
        print("\n[3] Silhouette k-sweep (k=2..10)...")
        k_range = range(2, 11)
        sil, inertias = silhouette_sweep(X_pca, k_range)
        best_k = list(k_range)[int(np.argmax(sil))]
        print(f"  Best k by silhouette: {best_k}")
        plot_silhouette(k_range, sil, inertias,
                         FIG_PATH / 'silhouette')

    # ── Step 5: Microstate clustering + (optional) PCCA+ coarse-graining ──
    # Default to 200 microstates if user enabled --pcca or --n-microstates
    use_microstate = (args.pcca or args.n_microstates is not None
                      or args.reduce == 'tica')  # TICA path → microstates
    if args.n_microstates is None:
        n_micro = 200 if use_microstate else args.k
    else:
        n_micro = args.n_microstates

    if use_microstate and n_micro > args.k:
        # MSM-style: microstate clustering → MSM → PCCA+ → macrostates
        if args.clustering == 'regspace':
            mode_str = ('regspace [2D placement]' if args.regspace_2d_only
                        else 'regspace [N-D placement]')
            print(f"\n[4a] Microstate clustering: {mode_str} "
                  f"targeting ~{n_micro} centers...")
            micro_labels, micro_centers, _ = cluster_regspace(
                X_pca, target_n=n_micro,
                use_2d_only=args.regspace_2d_only)
            n_actual = len(micro_centers)
            if n_actual != n_micro:
                print(f"  Got {n_actual} centers (binary search converged "
                      f"to ~{n_micro})")
                n_micro = n_actual
        else:
            print(f"\n[4a] Microstate clustering: K-means with "
                  f"k={n_micro}...")
            micro_labels, micro_centers, _ = cluster_kmeans(X_pca, n_micro)
        populations = np.bincount(micro_labels, minlength=n_micro)
        print(f"  Microstate populations: min={populations.min()}, "
              f"max={populations.max()}, mean={populations.mean():.1f}")

        # Microstate landscape plot
        plot_microstate_centers(X_pca[:, :2], micro_centers[:, :2],
                                  populations,
                                  outpath=FIG_PATH / 'microstates_landscape',
                                  axis_label='IC' if reduce_method == 'tica'
                                                 else 'PC',
                                  method_name=reduce_method.upper())

        # Build MSM
        msm_lag = args.msm_lag if args.msm_lag else args.tica_lag
        print(f"\n[4b] Building MSM at lag {msm_lag} frames "
              f"({msm_lag * 0.01:.2f} ns)...")
        T, pi, C = build_msm_transition_matrix(
            micro_labels, rep_boundaries, msm_lag, n_micro)
        print(f"  Transition matrix: {T.shape}, "
              f"min populated states: "
              f"{(C.sum(axis=1) > 0).sum()}/{n_micro}")
        np.savez(DATA_PATH / f'msm_{suffix}.npz',
                 T=T, pi=pi, C=C,
                 micro_labels=micro_labels,
                 micro_centers=micro_centers,
                 lag=msm_lag)

        # PCCA+ → macrostates
        print(f"\n[4c] PCCA+ coarse-graining to {args.k} macrostates...")
        macro_assign, pcca_eigvals = pcca_plus(T, args.k)
        print(f"  MSM eigenvalues (top {args.k}): "
              f"{['%.3f' % v for v in pcca_eigvals[:args.k]]}")
        # Imply timescales from MSM eigenvalues
        msm_ts = implied_timescales(pcca_eigvals[1:args.k], msm_lag)
        print(f"  Macrostate implied timescales: "
              f"{['%.2f ns' % t if not np.isnan(t) else 'n/a' for t in msm_ts]}")

        # Map every frame's microstate to its macrostate
        labels = macro_assign[micro_labels]

        # Plot microstate-to-macrostate decomposition
        plot_microstate_to_macro(
            X_pca[:, :2], micro_labels, labels,
            micro_centers[:, :2], macro_assign,
            outpath=FIG_PATH / 'pcca_decomposition',
            axis_label='IC' if reduce_method == 'tica' else 'PC',
            method_name=reduce_method.upper())

        # ── CK test (Chapman-Kolmogorov MSM validation) ──
        if args.ck_test:
            print(f"\n[4d] Chapman-Kolmogorov test "
                  f"(τ={msm_lag} frames, "
                  f"checking k=1..{args.ck_nlags})...")
            k_arr, T_obs, T_pred = ck_test(
                micro_labels, macro_assign, rep_boundaries,
                base_lag=msm_lag, n_lag_multiples=args.ck_nlags)

            T_obs_lower = T_obs_upper = None
            if args.ck_bootstrap > 0:
                print(f"  Bootstrap CK ({args.ck_bootstrap} resamples)...")
                T_obs_lower, T_obs_upper = ck_test_bootstrap(
                    micro_labels, macro_assign, rep_boundaries,
                    base_lag=msm_lag,
                    n_lag_multiples=args.ck_nlags,
                    n_bootstrap=args.ck_bootstrap)

            plot_ck_test(k_arr, T_obs, T_pred, msm_lag,
                          frame_dt_ns=0.01,
                          outpath=FIG_PATH / 'ck_test',
                          T_obs_lower=T_obs_lower,
                          T_obs_upper=T_obs_upper)
            print(f"  Saved: {FIG_PATH / 'ck_test.png'}")

            # Quick goodness-of-fit summary (mean absolute deviation
            # between observed and predicted transition probabilities)
            mad = np.nanmean(np.abs(T_obs - T_pred))
            print(f"  Mean |T_obs - T_pred| = {mad:.4f}  "
                  f"(<0.05 = excellent, <0.10 = good, >0.20 = MSM is invalid)")

        # Centroid frames per macrostate (closest to macrostate mean)
        centroid_indices = find_centroid_frames(X_pca, labels, args.k)
        centroids = np.array([X_pca[labels == c].mean(axis=0)
                               for c in range(args.k)])
    else:
        # Geometric mode: direct K-means with args.k (original behavior)
        print(f"\n[4] K-means (geometric, no MSM/PCCA+) with k={args.k}...")
        labels, centroids, km = cluster_kmeans(X_pca, args.k)
        print(f"  Cluster sizes: {np.bincount(labels)}")
        centroid_indices = find_centroid_frames(X_pca, labels, args.k)

    # Optional: replace density-centroid representatives with evenly-spread ones.
    # (Keeps any MSM/ITS/CK diagnostics above; only changes which frames are the
    #  representative states and re-tiles the landscape into non-overlapping cells.)
    if args.rep_mode == 'spread':
        print(f"\n[4b] Re-selecting {args.k} representatives by SPREAD "
              f"(farthest-point sampling + Voronoi tiling)...")
        centroid_indices, labels = find_spread_frames(X_pca[:, :2], args.k)
        centroids = X_pca[centroid_indices]
        print(f"  Spread state sizes: {np.bincount(labels, minlength=args.k)}")

    # Save assignments
    out_assign = DATA_PATH / f'cluster_assignments_{suffix}.npz'
    np.savez(out_assign,
             labels=labels,
             centroid_frames=np.array([
                 metadata[idx] if idx is not None else (-1, -1, -1)
                 for idx in centroid_indices]),
             feature_names=np.array(feature_names),
             X_raw=X_raw)
    print(f"  Saved: {out_assign}")

    # ── Step 6: plots ──
    print("\n[5] Generating plots...")
    axis_label = 'IC' if reduce_method == 'tica' else 'PC'
    method_name = 'TICA' if reduce_method == 'tica' else 'PCA'
    plot_pca_clusters(X_pca, labels, args.k, centroid_indices,
                       FIG_PATH / 'pca_clusters',
                       axis_label=axis_label)
    plot_state_distance_profile(X_raw, labels, feature_names, args.k,
                                 FIG_PATH / 'distance_profile')
    plot_population(labels, args.k,
                     FIG_PATH / 'population')
    plot_replicate_state_distribution(metadata, labels, args.k,
                                       FIG_PATH / 'replicate_state')

    # New landscape plot: density + state boundaries + centroids
    centroid_2d = np.zeros((args.k, 2))
    for c in range(args.k):
        m = labels == c
        if m.sum() > 0:
            centroid_2d[c] = X_pca[m, :2].mean(axis=0)
    plot_density_with_state_boundaries(
        X_pca[:, :2], labels, args.k, centroid_2d,
        outpath=FIG_PATH / 'landscape_with_states',
        axis_label=axis_label, method_name=method_name)

    # Save data needed for replot
    np.savez(DATA_PATH / f'landscape_{suffix}.npz',
             X_pca_2d=X_pca[:, :2],
             labels=labels,
             centroid_2d=centroid_2d,
             k=args.k,
             axis_label=axis_label,
             method_name=method_name)

    # ── Step 7: extract representative PDBs ──
    if not args.skip_pdb:
        print("\n[6] Extracting representative frames...")
        rep_dir = REP_PATH / _date.today().isoformat() / run_tag
        rep_dir.mkdir(parents=True, exist_ok=True)
        rep_info = []
        for c, idx in enumerate(centroid_indices):
            if idx is None:
                continue
            seed, sample, frame_idx = metadata[idx]
            pop = (labels == c).sum()
            pct = pop / len(labels) * 100
            out_pdb = rep_dir / f'state_{c+1}.pdb'
            ok = extract_pdb_for_frame(set_key, seed, sample, frame_idx,
                                        out_pdb)
            status = 'OK' if ok else 'FAIL'
            print(f"  State {c+1}: seed-{seed}_sample-{sample} "
                  f"frame {frame_idx} → {out_pdb.name} ({status}) "
                  f"[{pop} frames, {pct:.1f}%]")
            rep_info.append({
                'state': c + 1,
                'population': int(pop),
                'fraction': float(pct / 100),
                'representative_seed': int(seed),
                'representative_sample': int(sample),
                'representative_frame': int(frame_idx),
                'representative_pdb': str(out_pdb.relative_to(CWD)),
                'mean_features': {fn: float(X_raw[labels == c][:, j].mean())
                                  for j, fn in enumerate(feature_names)},
            })

        # Save state info JSON
        out_info = rep_dir / 'state_info.json'
        with open(out_info, 'w') as f:
            json.dump({
                'set': set_key,
                'k': args.k,
                'feature_names': feature_names,
                'states': rep_info,
            }, f, indent=2)
        print(f"  Saved: {out_info}")

    # ── Summary ──
    print(f"\n{'=' * 60}")
    print(f"State summary ({args.k} clusters, {len(labels)} frames total)")
    print(f"{'=' * 60}")
    print(f"  {'State':<7} {'N':>6} {'%':>6}  Mean features (nm)")
    for c in range(args.k):
        mask = labels == c
        n = mask.sum()
        pct = n / len(labels) * 100
        means = X_raw[mask].mean(axis=0)
        means_str = ' '.join(f'{fn[:8]}:{v:.2f}'
                              for fn, v in zip(feature_names, means))
        print(f"  {c+1:<7} {n:>6} {pct:>5.1f}%  {means_str}")


if __name__ == '__main__':
    main()
