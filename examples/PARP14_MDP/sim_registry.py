#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Shared simulation registry for PARP14 CALVADOS analysis scripts.

This is the single source of truth for:
  - canonical domain-unit / active-site / domain-boundary definitions
  - the named simulation sets (fl, md, core, mka, norrm, noart, md3art, fl_optimized)
  - FL -> construct residue-number mapping (so active sites / domains defined in
    full-length numbering get remapped into any sub-construct automatically)
  - auto-discovery of fragment simulations under fragments/
  - registration of ARBITRARY simulation folders via --sim-folder

Why this exists
---------------
Every analysis script used to re-implement the same DOMAIN_UNITS table, the same
`{set}/seed-{1-5}_sample-{0-4}/{sysname}.dcd` path logic, and the same FL->construct
remap. Adding a new construct meant editing every script. Now a NEW simulation can
be analyzed by ANY script with `--sim-folder PATH`, with no code edits.

A simulation folder is any directory laid out as:

    PATH/
      seed-{1-5}_sample-{0-4}/
        top.pdb
        <sysname>.dcd
      metadata.json        (optional but recommended)
      input/domains.yaml   (optional; restraint domains in construct numbering)

To map active sites / domains, the registry needs to know which FL domain UNITS the
construct contains. It resolves them in this order:
  1. explicit `units=[...]` (e.g. from --units on the CLI)
  2. metadata.json  -> "units": [...]
  3. if the folder basename matches a known named set -> that set's units
  4. give up with a clear error (pass --units or write metadata.json)

`write_metadata()` lets prepare scripts (or `stamp_metadata.py`) drop a metadata.json
into every folder so step 2 always works and nothing extra is needed at analysis time.
"""

import os
import glob
import json

try:
    import yaml
except Exception:  # yaml is optional (only needed to read domains.yaml)
    yaml = None

CWD = os.path.dirname(os.path.abspath(__file__))

SEEDS = list(range(1, 6))
SAMPLES = list(range(0, 5))

# ============================================================
# Canonical definitions (full-length numbering)
# ============================================================

# 11 grouped domain units used to build the combinatorial library.
DOMAIN_UNITS = {
    'rrm1': (1, 145), 'rrm2': (146, 224), 'rrm3': (225, 314),
    'kh1-kh6': (315, 737), 'kh7a': (738, 789),
    'md1l1': (790, 1004), 'md2': (1004, 1193), 'md3': (1207, 1388),
    'khb-kh8': (1389, 1533), 'wwe': (1534, 1602), 'art': (1603, 1801),
}

# Active sites in FULL-LENGTH numbering (parp14/input/active_sites.yaml).
ACTIVE_SITES_FL = {
    'MD1': {'catalytic': [831, 923, 962],
            'pocket': [822, 823, 824, 825, 826, 827, 828, 829, 830, 831, 832, 833,
                       834, 835, 836, 919, 920, 921, 922, 923, 924, 925, 926, 927,
                       961, 962, 966]},
    'MD2': {'catalytic': [1035, 1046, 1134, 1171],
            'pocket': [1021, 1022, 1023, 1024, 1034, 1035, 1036, 1037, 1038, 1039,
                       1040, 1041, 1042, 1043, 1044, 1045, 1046, 1047, 1130, 1131,
                       1132, 1133, 1134, 1135, 1136, 1137, 1138, 1139, 1140, 1141,
                       1170, 1171, 1175, 1178]},
    'MD3': {'catalytic': [1248, 1259, 1330, 1371],
            'pocket': [1235, 1236, 1237, 1247, 1248, 1249, 1250, 1251, 1252, 1253,
                       1254, 1255, 1256, 1257, 1258, 1259, 1260, 1261, 1302, 1303,
                       1304, 1324, 1325, 1326, 1327, 1328, 1329, 1330, 1331, 1332,
                       1333, 1334, 1335, 1336, 1337, 1369, 1370, 1371, 1375]},
    'ART': {'catalytic': [1684, 1705, 1706, 1722],
            'pocket': [1681, 1682, 1683, 1684, 1685, 1688, 1701, 1704, 1705, 1706,
                       1707, 1708, 1709, 1714, 1715, 1716, 1721, 1722, 1726, 1727,
                       1781]},
}

SITE_NAMES = ['MD1', 'MD2', 'MD3', 'ART']
SITE_COLORS = {'MD1': '#e6194b', 'MD2': '#3cb44b', 'MD3': '#4363d8', 'ART': '#f58231'}

# Maps each FL domain unit -> the active site it carries (for sub-construct detection).
UNIT_TO_SITE = {'md1l1': 'MD1', 'md2': 'MD2', 'md3': 'MD3', 'art': 'ART'}

# Energy-analysis domain boundaries (FL numbering): trimmed by 3 residues at
# zero/small-gap borders to prevent steric-clash artifacts between adjacent domains.
FL_DOMAINS = {
    'RRM1':    (6, 88),
    'RRM2':    (150, 223),
    'RRM3':    (227, 301),
    'KH1-6':   (315, 734),
    'KH7a':    (741, 786),
    'MD1':     (794, 978),
    'MD2':     (1003, 1190),
    'MD3':     (1216, 1384),
    'KHb-KH8': (1392, 1530),
    'WWE':     (1537, 1601),
    'ART':     (1605, 1801),
}

# ============================================================
# Named simulation sets
# ============================================================
# Each set: sysname (dcd basename), label, color, and the FL domain units present.
# `units` drives the FL->construct mapping for every derived quantity.

SETS = {
    'fl': {
        'sysname': 'parp14', 'label': 'Full Length (1-1801)', 'color': '#1f77b4',
        'units': list(DOMAIN_UNITS.keys()),
    },
    'norrm': {
        'sysname': 'parp14_norrm', 'label': 'KH1-ART (315-1801)', 'color': '#9467bd',
        'units': ['kh1-kh6', 'kh7a', 'md1l1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art'],
    },
    'noart': {
        'sysname': 'parp14_noart', 'label': 'KH1-WWE (315-1602)', 'color': '#8c564b',
        'units': ['kh1-kh6', 'kh7a', 'md1l1', 'md2', 'md3', 'khb-kh8', 'wwe'],
    },
    'core': {
        'sysname': 'parp14_core', 'label': 'KH7-ART (738-1801)', 'color': '#2ca02c',
        'units': ['kh7a', 'md1l1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art'],
    },
    'mka': {
        'sysname': 'parp14_mka', 'label': 'MD1-ART (790-1801)', 'color': '#d62728',
        'units': ['md1l1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art'],
    },
    'md': {
        'sysname': 'parp14_macrodomains', 'label': 'MD1-MD3 (790-1388)', 'color': '#ff7f0e',
        'units': ['md1l1', 'md2', 'md3'],
    },
    'md3art': {
        'sysname': 'parp14_md3art', 'label': 'MD3-ART (1207-1801)', 'color': '#17becf',
        'units': ['md3', 'khb-kh8', 'wwe', 'art'],
    },
    'fl_optimized': {
        'sysname': 'parp14', 'label': 'FL (optimized restraints)', 'color': '#aec7e8',
        'units': list(DOMAIN_UNITS.keys()),
    },
    'fl_go': {
        'sysname': 'parp14', 'label': 'FL (Go-model restraints, 2us extensions)', 'color': '#7f7f7f',
        'units': list(DOMAIN_UNITS.keys()),
    },
}

# Restraint-domain boundaries (CONSTRUCT numbering) for the named sets. Used by the
# FNC module, which reads restraint domains straight from how the sim was built.
# Most sets are derivable from `units`, but these are the exact as-simulated values.
CONSTRUCT_DOMAINS = {
    'fl':    {'parp14':              [[6, 88], [150, 223], [227, 301], [791, 978], [1003, 1190], [1216, 1387], [1523, 1601], [1605, 1801]]},
    'md':    {'parp14_macrodomains': [[2, 189], [214, 401], [414, 585]]},
    'core':  {'parp14_core':         [[54, 241], [266, 453], [466, 637], [773, 851], [855, 1051]]},
    'mka':   {'parp14_mka':          [[2, 189], [214, 401], [414, 585], [721, 799], [803, 999]]},
    'norrm': {'parp14_norrm':        [[477, 664], [689, 876], [889, 1060], [1196, 1274], [1278, 1474]]},
    'noart': {'parp14_noart':        [[477, 664], [689, 876], [889, 1060], [1196, 1274]]},
    'md3art':{'parp14_md3art':       [[10, 181], [317, 395], [399, 595]]},
    'fl_optimized': {'parp14': [[16, 78], [156, 214], [235, 304], [320, 379], [390, 449], [460, 515],
                                [531, 583], [599, 660], [676, 727], [738, 789], [800, 968],
                                [1015, 1183], [1217, 1378], [1389, 1461], [1462, 1533],
                                [1549, 1587], [1613, 1791]]},
    # fl_go's input/domains.yaml has these identical boundaries (same AF2 structure,
    # same restraint domain trims); it differs only in restraint TYPE (custom Go-model
    # contacts, input/custom_restraints_go.txt) rather than harmonic domain restraints.
    'fl_go': {'parp14': [[16, 78], [156, 214], [235, 304], [320, 379], [390, 449], [460, 515],
                         [531, 583], [599, 660], [676, 727], [738, 789], [800, 968],
                         [1015, 1183], [1217, 1378], [1389, 1461], [1462, 1533],
                         [1549, 1587], [1613, 1791]]},
}

# Holds absolute directories for folders registered via register_external_folder().
# Populated in a process before any ProcessPoolExecutor is created so forked
# workers inherit it (Linux fork start method).
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


# ============================================================
# FL -> construct residue mapping
# ============================================================

def compute_fl_blocks(unit_names):
    """Merge the FL residue ranges of the given units into contiguous blocks."""
    ranges = sorted([DOMAIN_UNITS[n] for n in unit_names])
    merged = [list(ranges[0])]
    for s, e in ranges[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]


def build_fl_to_construct_map(unit_names):
    """Return a function mapping an FL residue number to construct numbering
    (1-based, contiguous) for a construct containing `unit_names`."""
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


def get_units(set_key):
    """FL domain units present in a set/fragment/external folder."""
    return SETS[set_key]['units']


def get_construct_sites(set_key):
    """Active-site names present in a construct (subset of SITE_NAMES)."""
    units = get_units(set_key)
    return [UNIT_TO_SITE[u] for u in units if u in UNIT_TO_SITE]


def get_active_sites_for_set(set_key):
    """Active sites with residue numbers remapped into the set's construct numbering."""
    if set_key in ('fl', 'fl_optimized', 'fl_go'):
        # Both are the full uncompressed 1801-residue sequence (fl_optimized
        # differs only in per-domain restraint trims, not numbering) -- must
        # bypass build_fl_to_construct_map(), which compresses out gaps
        # between DOMAIN_UNITS (e.g. the 1194-1206 MD2/MD3 linker) under the
        # assumption those residues are genuinely absent, as they are for a
        # real sub-construct like 'noart'. That assumption is false here and
        # silently shifted every MD3/ART residue number by -13.
        return ACTIVE_SITES_FL
    fl_to_c = build_fl_to_construct_map(get_units(set_key))
    sites = {}
    for sname in get_construct_sites(set_key):
        data = ACTIVE_SITES_FL[sname]
        sites[sname] = {
            'catalytic': [r for r in (fl_to_c(x) for x in data['catalytic']) if r is not None],
            'pocket': [r for r in (fl_to_c(x) for x in data['pocket']) if r is not None],
        }
    return sites


def get_construct_domains_for_set(set_key):
    """FL_DOMAINS remapped into the set's construct numbering (energy boundaries)."""
    if set_key in ('fl', 'fl_optimized', 'fl_go'):
        return FL_DOMAINS
    fl_to_c = build_fl_to_construct_map(get_units(set_key))
    mapped = {}
    for dname, (fl_s, fl_e) in FL_DOMAINS.items():
        c_s = fl_to_c(fl_s)
        c_e = fl_to_c(fl_e)
        if c_s is not None and c_e is not None:
            mapped[dname] = (c_s, c_e)
    return mapped


# ============================================================
# Path resolution
# ============================================================

def get_sim_dir(set_key, seed, sample):
    """Absolute path to one replicate directory for any registered set."""
    if set_key in EXTERNAL_DIRS:
        base = EXTERNAL_DIRS[set_key]
    elif set_key.startswith('frag_'):
        base = os.path.join(CWD, 'fragments', set_key[5:])
    else:
        base = os.path.join(CWD, set_key)
    sim_dir = os.path.join(base, f'seed-{seed}_sample-{sample}')
    if set_key in EXTERNAL_DIRS and not os.path.isdir(sim_dir):
        if set_key not in _FLAT_DIRS_CACHE:
            _FLAT_DIRS_CACHE[set_key] = _flat_replicate_dirs(base, SETS[set_key]['sysname'])
        flat = _FLAT_DIRS_CACHE[set_key]
        idx = (seed - 1) * len(SAMPLES) + sample
        if idx < len(flat):
            sim_dir = flat[idx]
    return sim_dir


def get_sim_paths(set_key, seed, sample):
    """Return (pdb, dcd) for one replicate of any registered set."""
    sysname = SETS[set_key]['sysname']
    sim_dir = get_sim_dir(set_key, seed, sample)
    dcd = os.path.join(sim_dir, f'{sysname}.dcd')
    for name in ('top.pdb', 'restart.pdb', 'checkpoint.pdb'):
        pdb = os.path.join(sim_dir, name)
        if os.path.isfile(pdb):
            return pdb, dcd
    return os.path.join(sim_dir, 'top.pdb'), dcd


def replicate_jobs(active_sets, seeds=SEEDS, samples=SAMPLES):
    """All (set_key, seed, sample) tuples whose pdb+dcd exist on disk."""
    jobs = []
    for set_key in active_sets:
        for seed in seeds:
            for sample in samples:
                pdb, dcd = get_sim_paths(set_key, seed, sample)
                if os.path.isfile(pdb) and os.path.isfile(dcd):
                    jobs.append((set_key, seed, sample))
    return jobs


# ============================================================
# Metadata I/O
# ============================================================

def metadata_path(folder):
    return os.path.join(folder, 'metadata.json')


def read_metadata(folder):
    """Return the parsed metadata.json for a folder, or None."""
    p = metadata_path(folder)
    if not os.path.isfile(p):
        return None
    with open(p) as f:
        return json.load(f)


def write_metadata(folder, units, sysname=None, label=None, extra=None):
    """Write a metadata.json describing a simulation folder so analysis scripts
    can analyze it via --sim-folder with no other arguments.

    units   : list of FL domain-unit keys present (subset of DOMAIN_UNITS keys)
    sysname : dcd basename (auto-detected from the folder if omitted)
    """
    units = [u.lower() for u in units]
    bad = [u for u in units if u not in DOMAIN_UNITS]
    if bad:
        raise ValueError(f"unknown domain units {bad}; valid: {list(DOMAIN_UNITS)}")
    if sysname is None:
        sysname = detect_sysname(folder)
    meta = {
        'units': units,
        'sysname': sysname,
        'sites': [UNIT_TO_SITE[u] for u in units if u in UNIT_TO_SITE],
    }
    if label:
        meta['label'] = label
    if extra:
        meta.update(extra)
    with open(metadata_path(folder), 'w') as f:
        json.dump(meta, f, indent=2)
    return meta


def detect_sysname(folder):
    """Best-effort sysname (dcd basename) detection for a folder.

    Tries, in order: a seed-*/*.dcd basename, then the first key of
    input/domains.yaml, then 'parp14'.
    """
    for sd in sorted(glob.glob(os.path.join(folder, 'seed-*_sample-*'))):
        dcds = glob.glob(os.path.join(sd, '*.dcd'))
        if dcds:
            return os.path.splitext(os.path.basename(dcds[0]))[0]
    dyaml = os.path.join(folder, 'input', 'domains.yaml')
    if yaml is not None and os.path.isfile(dyaml):
        with open(dyaml) as f:
            d = yaml.safe_load(f) or {}
        if d:
            return next(iter(d))
    return 'parp14'


# ============================================================
# Discovery / registration
# ============================================================

def register_fragment(frag_name, metadata):
    """Register a contiguous fragment (under fragments/) as a set."""
    set_key = f'frag_{frag_name}'
    units = [u.lower() for u in metadata.get('units', [])]
    SETS[set_key] = {
        'sysname': metadata.get('sysname', 'parp14'),
        'label': frag_name,
        'color': '#888888',
        'units': units,
    }
    CONSTRUCT_DOMAINS[set_key] = {
        SETS[set_key]['sysname']: metadata.get('domain_ranges_construct', [])
    }
    return set_key


def discover_fragments():
    """Scan fragments/ and register every fragment that has a trajectory.
    Returns the list of registered set keys (with 'frag_' prefix)."""
    frag_dir = os.path.join(CWD, 'fragments')
    if not os.path.isdir(frag_dir):
        return []
    discovered = []
    for name in sorted(os.listdir(frag_dir)):
        fdir = os.path.join(frag_dir, name)
        meta = read_metadata(fdir)
        if meta is None:
            continue
        has_traj = any(
            glob.glob(os.path.join(fdir, entry, '*.dcd'))
            for entry in os.listdir(fdir) if entry.startswith('seed-')
        )
        if not has_traj:
            continue
        discovered.append(register_fragment(name, meta))
    return discovered


def _resolve_units(folder, units=None):
    """Determine FL domain units for a folder (see module docstring for order)."""
    if units:
        return [u.lower() for u in units]
    meta = read_metadata(folder)
    if meta and meta.get('units'):
        return [u.lower() for u in meta['units']]
    base = os.path.basename(os.path.normpath(folder))
    if base in SETS:
        return list(SETS[base]['units'])
    raise ValueError(
        f"cannot determine domain units for '{folder}'. Pass --units, or add a "
        f"metadata.json (see write_metadata / stamp_metadata.py). "
        f"Valid units: {list(DOMAIN_UNITS)}")


def register_external_folder(path, units=None, sysname=None, name=None,
                             label=None, color='#444444'):
    """Register an arbitrary simulation folder so any script can analyze it.

    Returns the set_key under which it is registered. The folder is expected to
    contain seed-{1-5}_sample-{0-4}/ replicate dirs with top.pdb + <sysname>.dcd.
    """
    folder = os.path.abspath(path)
    if not os.path.isdir(folder):
        raise FileNotFoundError(f"sim folder not found: {folder}")
    resolved_units = _resolve_units(folder, units)
    meta = read_metadata(folder) or {}
    if sysname is None:
        sysname = meta.get('sysname') or detect_sysname(folder)
    set_key = name or os.path.basename(os.path.normpath(folder))
    SETS[set_key] = {
        'sysname': sysname,
        'label': label or meta.get('label') or set_key,
        'color': color,
        'units': resolved_units,
    }
    # Restraint domains (construct numbering) for the FNC module, if available.
    dranges = meta.get('domain_ranges_construct')
    if dranges is None:
        dyaml = os.path.join(folder, 'input', 'domains.yaml')
        if yaml is not None and os.path.isfile(dyaml):
            with open(dyaml) as f:
                d = yaml.safe_load(f) or {}
            dranges = d.get(sysname) or (next(iter(d.values())) if d else None)
    if dranges is not None:
        CONSTRUCT_DOMAINS[set_key] = {sysname: dranges}
    EXTERNAL_DIRS[set_key] = folder
    return set_key


# ============================================================
# CLI helpers
# ============================================================

def add_sim_folder_args(parser):
    """Add the standard --sim-folder / --units arguments to an ArgumentParser."""
    parser.add_argument(
        '--sim-folder', nargs='+', default=None, metavar='PATH',
        help='One or more simulation-set folders to analyze (each containing '
             'seed-*_sample-*/ replicates). Domain units are read from the '
             "folder's metadata.json, or pass --units.")
    parser.add_argument(
        '--units', nargs='+', default=None, metavar='UNIT',
        help='FL domain units in the --sim-folder construct (e.g. md1l1 md2 md3). '
             'Required only if the folder has no metadata.json. '
             f'Valid: {", ".join(DOMAIN_UNITS)}')
    return parser


def register_cli_folders(args):
    """Register every folder passed via --sim-folder; return their set keys.

    Call this once, in the main process, before launching any worker pool.
    """
    keys = []
    folders = getattr(args, 'sim_folder', None)
    if not folders:
        return keys
    units = getattr(args, 'units', None)
    for folder in folders:
        keys.append(register_external_folder(folder, units=units))
    return keys
