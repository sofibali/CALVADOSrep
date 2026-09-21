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
    # 'md1' spans macrodomain 1 (UniProt Macro 1, 791-978) PLUS the linker that
    # runs from the end of that domain to the start of MD2. It was called
    # 'md1l1' while that linker was thought to be missing; the linker is in
    # fact present, so the unit is just MD1-with-its-linker. 'md1l1' is kept as
    # a read alias below so the existing metadata.json files and
    # representative_frames/ directory names still resolve.
    #
    # It ends at 1003, not 1004: MD2 starts at 1004, and a residue cannot
    # belong to two units (the old 790-1004 / 1004-1193 pair double-counted
    # 1004, so the 11 units summed to 1789 residues but covered only 1788).
    'md1': (790, 1003), 'md2': (1004, 1193), 'md3': (1207, 1388),
    'khb-kh8': (1389, 1533), 'wwe': (1534, 1602), 'art': (1603, 1801),
}

# Deprecated unit spellings -> canonical key. Historical metadata.json files,
# fragment names and representative_frames/ directories still carry 'md1l1';
# resolving rather than rewriting them keeps that corpus readable.
UNIT_ALIASES = {'md1l1': 'md1'}


def canonical_unit(name):
    """Canonical DOMAIN_UNITS key for a possibly-deprecated unit spelling."""
    n = str(name).lower()
    return UNIT_ALIASES.get(n, n)


def canonical_units(names):
    return [canonical_unit(n) for n in names]

# Active sites in FULL-LENGTH numbering.
#
# Loaded from parp14/input/active_sites.yaml rather than duplicated here. That
# file is the single source of truth; it was previously copied into six
# different modules, which is exactly how they drifted out of sync (and how the
# mislabelled catalytic residues survived so long -- see
# examples/PARP14_MDP/docs/NUMBERING_AUDIT.md).
_ACTIVE_SITES_YAML = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    'parp14', 'input', 'active_sites.yaml')


def _load_active_sites(path=_ACTIVE_SITES_YAML):
    import yaml as _yaml
    with open(path) as fh:
        raw = _yaml.safe_load(fh)
    out = {}
    for name, info in raw.items():
        if not isinstance(info, dict) or 'catalytic_residues' not in info:
            continue
        out[name] = {'catalytic': list(info['catalytic_residues']),
                     'pocket': sorted(set(info.get('pocket_residues', []))),
                     'domain_range': info.get('domain_range')}
    return out


ACTIVE_SITES_FL = _load_active_sites()

SITE_NAMES = ['MD1', 'MD2', 'MD3', 'WWE', 'ART']
SITE_COLORS = {'MD1': '#e6194b', 'MD2': '#3cb44b', 'MD3': '#4363d8', 'WWE': '#911eb4', 'ART': '#f58231'}

# Maps each FL domain unit -> the active site it carries (for sub-construct detection).
UNIT_TO_SITE = {'md1': 'MD1', 'md2': 'MD2', 'md3': 'MD3', 'wwe': 'WWE', 'art': 'ART'}

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
        'units': ['kh1-kh6', 'kh7a', 'md1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art'],
    },
    'noart': {
        'sysname': 'parp14_noart', 'label': 'KH1-WWE (315-1602)', 'color': '#8c564b',
        'units': ['kh1-kh6', 'kh7a', 'md1', 'md2', 'md3', 'khb-kh8', 'wwe'],
    },
    'core': {
        'sysname': 'parp14_core', 'label': 'KH7-ART (738-1801)', 'color': '#2ca02c',
        'units': ['kh7a', 'md1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art'],
    },
    'mka': {
        'sysname': 'parp14_mka', 'label': 'MD1-ART (790-1801)', 'color': '#d62728',
        'units': ['md1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art'],
    },
    'md': {
        'sysname': 'parp14_macrodomains', 'label': 'MD1-MD3 (790-1388)', 'color': '#ff7f0e',
        'units': ['md1', 'md2', 'md3'],
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
        'sysname': 'parp14', 'label': 'FL (Go-model restraints, 2us extensions)', 'color': '#CC79A7',
        'units': list(DOMAIN_UNITS.keys()),
    },
    'md_full': {
        'sysname': 'parp14_md1md3_full', 'label': 'MD1-MD3 contiguous (790-1388, 599 res)',
        'color': '#c8831c',
        'units': ['md1', 'md2', 'md3'],
    },
    'mka_full': {
        'sysname': 'parp14_mka_full', 'label': 'MD1-ART contiguous (790-1801, 1012 res)',
        'color': '#2845bd',
        'units': ['md1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art'],
    },
    'core_full_go': {
        'sysname': 'parp14_core_full_go',
        'label': 'KH7a-ART contiguous + Go-model KH7a-KHb restraints (738-1801, 1064 res)',
        'color': '#E69F00',
        'units': ['kh7a', 'md1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art'],
    },
    'kh1_art_full': {
        'sysname': 'parp14_kh1_art_full',
        'label': 'KH1-6-ART contiguous + Go-model KH7a-KHb restraints (315-1801, 1487 res)',
        'color': '#56B4E9',
        'units': ['kh1-kh6', 'kh7a', 'md1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art'],
    },
    'kh1_wwe_full': {
        'sysname': 'parp14_kh1_wwe_full',
        'label': 'KH1-6-WWE contiguous + Go-model KH7a-KHb restraints (315-1602, 1288 res, no ART)',
        'color': '#009E73',
        'units': ['kh1-kh6', 'kh7a', 'md1', 'md2', 'md3', 'khb-kh8', 'wwe'],
    },
    'md2_art_full': {
        'sysname': 'parp14_md2_art_full',
        'label': 'MD2-ART contiguous (1004-1801, 798 res)',
        'color': '#eb79eb',
        'units': ['md2', 'md3', 'khb-kh8', 'wwe', 'art'],
    },
    'md2_wwe_full': {
        'sysname': 'parp14_md2_wwe_full',
        'label': 'MD2-WWE contiguous (1004-1602, 599 res, no ART)',
        'color': '#a82366',
        'units': ['md2', 'md3', 'khb-kh8', 'wwe'],
    },
    'md3_art_full': {
        'sysname': 'parp14_md3_art_full',
        'label': 'MD3-ART contiguous (1207-1801, 595 res)',
        'color': '#1cc895',
        'units': ['md3', 'khb-kh8', 'wwe', 'art'],
    },
    'md3_wwe_full': {
        'sysname': 'parp14_md3_wwe_full',
        'label': 'MD3-WWE contiguous (1207-1602, 396 res, no ART)',
        'color': '#2387a8',
        'units': ['md3', 'khb-kh8', 'wwe'],
    },
    'core_wwe_full_go': {
        'sysname': 'parp14_core_wwe_full_go',
        'label': 'KH7a-WWE contiguous + Go-model KH7a-KHb restraints (738-1602, 865 res, no ART)',
        'color': '#0072B2',
        'units': ['kh7a', 'md1', 'md2', 'md3', 'khb-kh8', 'wwe'],
    },
    'mka_wwe_full': {
        'sysname': 'parp14_mka_wwe_full',
        'label': 'MD1L1-WWE contiguous (790-1602, 813 res, no ART)',
        'color': '#8a62e8',
        'units': ['md1', 'md2', 'md3', 'khb-kh8', 'wwe'],
    },
    'fl_wwe_full_go': {
        'sysname': 'parp14_fl_wwe_full_go',
        'label': 'FL-WWE contiguous + Go-model KH7a-KHb restraints (1-1602, 1602 res, no ART)',
        'color': '#D55E00',
        'units': ['rrm1', 'rrm2', 'rrm3', 'kh1-kh6', 'kh7a', 'md1', 'md2', 'md3', 'khb-kh8', 'wwe'],
    },
    'md1_md2': {
        'sysname': 'parp14_md1_md2',
        'label': 'MD1L1-MD2 contiguous, 2-domain isolation test (790-1193, 404 res)',
        'color': '#c8831c',
        'units': ['md1', 'md2'],
    },
    'md2_md3': {
        'sysname': 'parp14_md2_md3',
        'label': 'MD2-MD3 contiguous, 2-domain isolation test (1004-1388, 385 res)',
        'color': '#2845bd',
        'units': ['md2', 'md3'],
    },
}

# Colors above for the 13 "contiguous full-length family" sets + fl_go (the
# original 25-rep FL Go-restraint run, reused as this family's FL reference)
# come from sim_analysis/parp14_mdp_colors.py's CONSTRUCT_GO/CONSTRUCT_NOGO --
# validated colorblind-safe (OKLab CVD Delta E, --pairs all) within each of two
# groups that are never plotted in the same figure: Go-restraint sets
# (core_full_go, kh1_art_full, kh1_wwe_full, core_wwe_full_go, fl_wwe_full_go,
# fl_go) and non-Go sets (md_full, mka_full, md2_art_full, md2_wwe_full,
# md3_art_full, md3_wwe_full, mka_wwe_full). The old tab20-derived hexes here
# failed that validation outright (off the lightness band, chroma floor, and
# CVD separation checks all at once). md1_md2/md2_md3 reuse two of the
# non-Go-group hexes since they're always their own separate small figure.

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
    # From prepare_md_full.py / prepare_mka_full.py's own printed FL->local mapping
    # (single contiguous offset, local = FL - 789 -- see CONTIGUOUS_FL_RANGE above).
    'md_full': {'parp14_md1md3_full': [[2, 189], [214, 401], [427, 598]]},
    'mka_full': {'parp14_mka_full': [[2, 189], [214, 401], [427, 598], [734, 812], [816, 1012]]},
    # From prepare_core_full_go.py's own printed FL->local mapping (local = FL - 737).
    'core_full_go': {'parp14_core_full_go': [[54, 241], [266, 453], [479, 650], [786, 864], [868, 1064]]},
    # From prepare_extra_full_constructs.py's own printed FL->local mappings.
    'kh1_art_full': {'parp14_kh1_art_full': [[477, 664], [689, 876], [902, 1073], [1209, 1287], [1291, 1487]]},
    'kh1_wwe_full': {'parp14_kh1_wwe_full': [[477, 664], [689, 876], [902, 1073], [1209, 1287]]},
    'md2_art_full': {'parp14_md2_art_full': [[1, 187], [213, 384], [520, 598], [602, 798]]},
    'md2_wwe_full': {'parp14_md2_wwe_full': [[1, 187], [213, 384], [520, 598]]},
    'md3_art_full': {'parp14_md3_art_full': [[10, 181], [317, 395], [399, 595]]},
    'md3_wwe_full': {'parp14_md3_wwe_full': [[10, 181], [317, 395]]},
    'core_wwe_full_go': {'parp14_core_wwe_full_go': [[54, 241], [266, 453], [479, 650], [786, 864]]},
    'mka_wwe_full': {'parp14_mka_wwe_full': [[2, 189], [214, 401], [427, 598], [734, 812]]},
    'fl_wwe_full_go': {'parp14_fl_wwe_full_go': [[6, 88], [150, 223], [227, 301], [791, 978],
                                                  [1003, 1190], [1216, 1387], [1523, 1601]]},
    'md1_md2': {'parp14_md1_md2': [[2, 189], [214, 401]]},
    'md2_md3': {'parp14_md2_md3': [[1, 187], [213, 384]]},
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

# Sets that are a genuinely CONTIGUOUS slice of the FL sequence (no residues
# actually excised), keyed by set_key -> (fl_start, fl_end). This can't be
# inferred from `units` alone: a construct's unit list only says which
# DOMAIN_UNITS boundaries it covers, not whether the real simulated structure
# also kept the inter-unit linkers DOMAIN_UNITS doesn't assign to any unit
# (e.g. the 1194-1206 MD2/MD3 linker). 'fl'/'fl_optimized'/'fl_go' happen to
# cover ALL units so build_fl_to_construct_map()'s own all-units check catches
# them; md_full/mka_full are genuine sub-ranges (3 and 6 units) that were
# deliberately built to keep every linker within their span, so they need an
# explicit override here or compute_fl_blocks' gap-merging would wrongly
# compress out that same 1194-1206 linker for them too.
CONTIGUOUS_FL_RANGE = {
    'md_full': (790, 1388),
    'mka_full': (790, 1801),
    'core_full_go': (738, 1801),
    'kh1_art_full': (315, 1801),
    'kh1_wwe_full': (315, 1602),
    'md2_art_full': (1004, 1801),
    'md2_wwe_full': (1004, 1602),
    'md3_art_full': (1207, 1801),
    'md3_wwe_full': (1207, 1602),
    'core_wwe_full_go': (738, 1602),
    'mka_wwe_full': (790, 1602),
    'fl_wwe_full_go': (1, 1602),
    'md1_md2': (790, 1193),
    'md2_md3': (1004, 1388),
}


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


def build_fl_to_construct_map(unit_names, set_key=None):
    """Return a function mapping an FL residue number to construct numbering
    (1-based, contiguous) for a construct containing `unit_names`.

    If unit_names covers ALL 11 known domain units, the construct is
    full-length -- nothing was actually excised, so this must be the
    identity map. Without this check, compute_fl_blocks' gap-merging
    heuristic (designed for genuine sub-constructs that really do delete
    inter-unit linkers, e.g. 'md'/'core'/'norrm') incorrectly treats the
    real, un-excised inter-unit linkers within a full-length sequence (e.g.
    the 13-residue MD2-MD3 linker, residues 1194-1206, which isn't part of
    either the md2 or md3 unit's own boundaries) as deleted, shifting every
    downstream unit's mapped range by the gap size. This bit any
    --sim-folder-registered full-length sim (e.g. fl_go) even though the
    named 'fl'/'fl_optimized' sets were already correctly special-cased
    elsewhere to skip this mapping entirely.

    `set_key` lets a construct that ISN'T full-length (so the all-units check
    above doesn't fire) still declare itself contiguous via CONTIGUOUS_FL_RANGE
    (e.g. md_full/mka_full: 3 or 6 units, but deliberately built to keep every
    linker within their span -- the same gap-compression bug, just for a
    genuine sub-range instead of the whole sequence).
    """
    if set(unit_names) == set(DOMAIN_UNITS.keys()):
        return lambda fl_resid: fl_resid if 1 <= fl_resid <= 1801 else None
    if set_key in CONTIGUOUS_FL_RANGE:
        fl_start, fl_end = CONTIGUOUS_FL_RANGE[set_key]
        offset = 1 - fl_start
        return lambda fl_resid: (fl_resid + offset
                                  if fl_start <= fl_resid <= fl_end else None)
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


def build_construct_to_fl_map(unit_names, set_key=None):
    """Inverse of build_fl_to_construct_map: construct-local residue number ->
    FL residue number. Needed by anything that overlays construct-local data
    (e.g. a per-replicate contact map) onto a shared FL-numbered reference
    grid/domain-boundary set -- without this, local resi 1 of a sub-construct
    like 'md' (which is really FL 790, MD1L1) gets plotted at FL position 1
    (RRM1) instead, squeezing the whole construct into one corner of the grid.
    """
    if set(unit_names) == set(DOMAIN_UNITS.keys()):
        return lambda local_resid: local_resid if 1 <= local_resid <= 1801 else None
    if set_key in CONTIGUOUS_FL_RANGE:
        fl_start, fl_end = CONTIGUOUS_FL_RANGE[set_key]
        offset = 1 - fl_start
        return lambda local_resid: (local_resid - offset
                                     if fl_start + offset <= local_resid <= fl_end + offset
                                     else None)
    fl_blocks = compute_fl_blocks(unit_names)
    segments = []
    construct_pos = 1
    for fl_start, fl_end in fl_blocks:
        offset = construct_pos - fl_start
        segments.append((fl_start + offset, fl_end + offset, offset))
        construct_pos += (fl_end - fl_start + 1)

    def map_resid(local_resid):
        for local_s, local_e, off in segments:
            if local_s <= local_resid <= local_e:
                return local_resid - off
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
    fl_to_c = build_fl_to_construct_map(get_units(set_key), set_key=set_key)
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
    fl_to_c = build_fl_to_construct_map(get_units(set_key), set_key=set_key)
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
    units = canonical_units(units)
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
    units = canonical_units(metadata.get('units', []))
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
        return canonical_units(units)
    meta = read_metadata(folder)
    if meta and meta.get('units'):
        return canonical_units(meta['units'])
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
        help='FL domain units in the --sim-folder construct (e.g. md1 md2 md3). '
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
