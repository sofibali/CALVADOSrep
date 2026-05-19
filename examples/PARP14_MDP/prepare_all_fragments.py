#!/usr/bin/env python3
"""
Prepare CALVADOS simulations for all 66 contiguous PARP14 fragments.

Uses the same per-domain trim values as fl_optimized (RRM1 trim 10, RRM2 trim
10, ..., WWE trim 15, ART trim 10) but mapped to each construct's local
residue numbering. Inter-domain KH7a-KHb custom restraints are added when
BOTH KH7a and KHb are in the construct.

For each contiguous fragment present in alphafold_outputs/, prepares 25
replicates (5 seeds x 5 samples) using the seed-1_sample-0 AF3 model as the
input structure.

Output:
    fragments/{fragment_name}/seed-{1..5}_sample-{0..4}/
    fragments/{fragment_name}/input/  (shared domains, custom restraints)
    fragments/launch_commands.txt     (one launch command per fragment)

Usage:
    python prepare_all_fragments.py            # prepare all that have AF3
    python prepare_all_fragments.py --list     # list fragments + AF3 status
    python prepare_all_fragments.py --fragments md1l1_md2 md3_wwe_art
                                               # only specific fragments
"""

import os
import sys
import json
import yaml
import shutil
import argparse
import numpy as np
from pathlib import Path

# ============================================================
# Configuration
# ============================================================

CWD = Path(__file__).resolve().parent
PARP14_DIR = Path('/home/sbali/CALVADOS/parp14')
AF3_BASE = PARP14_DIR / 'alphafold_outputs'
FRAGMENTS_DIR = CWD / 'fragments'

FL_RESIDUES = CWD / 'input' / 'residues_CALVADOS3.csv'
DEFAULT_CONFIG = CWD.parent.parent / 'calvados' / 'data' / 'default_config.yaml'
FL_PDB = CWD / 'input' / 'parp14.pdb'

# Simulation parameters (match fl_optimized)
N_STEPS = 2_000_000          # 20 ns
WFREQ = 1_000
TEMP = 293
IONIC = 0.19
PH = 7.0
PLATFORM = 'CPU'
THREADS = 4
K_HARMONIC = 700.0
K_CUSTOM = 350.0
CUTOFF_RESTR = 0.9
N_CUSTOM_RES_CUTOFF_NM = 0.9   # for KH7a-KHb pair generation

SEEDS = range(1, 6)
SAMPLES = range(0, 5)

# ============================================================
# Domain definitions (FL numbering, must match prepare_fl_optimized.py)
# ============================================================

# Order matters for contiguous fragment generation
DOMAIN_UNITS = [
    ('rrm1',    (1, 145)),
    ('rrm2',    (146, 224)),
    ('rrm3',    (225, 314)),
    ('kh1-kh6', (315, 737)),
    ('kh7a',    (738, 789)),
    ('md1l1',   (790, 1004)),
    ('md2',     (1004, 1193)),
    ('md3',     (1207, 1388)),
    ('khb-kh8', (1389, 1533)),
    ('wwe',     (1534, 1602)),
    ('art',     (1603, 1801)),
]
UNIT_FL_RANGE = dict(DOMAIN_UNITS)
UNIT_ORDER = [u for u, _ in DOMAIN_UNITS]

# Restraint domain extents (full extents we trim from)
DOMAIN_EXTENTS = {
    'RRM1': (1, 145), 'RRM2': (146, 224), 'RRM3': (225, 314),
    'KH1': (315, 384), 'KH2': (385, 454), 'KH3': (455, 520),
    'KH4': (521, 593), 'KH5': (594, 665), 'KH6': (666, 737),
    'KH7a': (738, 789),
    'MD1L1': (790, 978),
    'MD2': (1005, 1193),
    'MD3': (1207, 1388),
    'KHb': (1389, 1461), 'KH8': (1462, 1533),
    'WWE': (1534, 1602), 'ART': (1603, 1801),
}

# RRM1 has a custom structured core (excludes disordered N/C)
STRUCTURED_EXTENTS = {'RRM1': (6, 88)}

# Per-domain trim values (same as fl_optimized)
DOMAIN_TRIMS = {
    'RRM1':  10, 'RRM2':  10, 'RRM3':  10,
    'KH1':    5, 'KH2':    5, 'KH3':    5,
    'KH4':   10, 'KH5':    5, 'KH6':   10,
    'KH7a':   0,
    'MD1L1': 10, 'MD2':   10, 'MD3':   10,
    'KHb':    0, 'KH8':    0,
    'WWE':   15, 'ART':   10,
}

# Map each restraint domain to its parent DOMAIN_UNIT (lowercase)
# If parent unit is not in the construct, the restraint domain is excluded.
RESTRAINT_TO_UNIT = {
    'RRM1': 'rrm1', 'RRM2': 'rrm2', 'RRM3': 'rrm3',
    'KH1': 'kh1-kh6', 'KH2': 'kh1-kh6', 'KH3': 'kh1-kh6',
    'KH4': 'kh1-kh6', 'KH5': 'kh1-kh6', 'KH6': 'kh1-kh6',
    'KH7a': 'kh7a',
    'MD1L1': 'md1l1', 'MD2': 'md2', 'MD3': 'md3',
    'KHb': 'khb-kh8', 'KH8': 'khb-kh8',
    'WWE': 'wwe', 'ART': 'art',
}


# ============================================================
# Helpers
# ============================================================

def all_contiguous_fragments():
    """Generate all 66 contiguous unit subsets."""
    fragments = []
    n = len(UNIT_ORDER)
    for i in range(n):
        for j in range(i, n):
            fragments.append(UNIT_ORDER[i:j+1])
    return fragments


def fragment_name(units):
    """Canonical name for AF3 directory lookup."""
    return '_'.join(units)


def af3_dir_for_fragment(units):
    """Locate AF3 output directory for this fragment. Returns Path or None.

    Prefers the canonical name (no timestamp). Falls back to most recent
    timestamped variant.
    """
    name = fragment_name(units)
    candidate = AF3_BASE / name
    if candidate.is_dir() and (candidate / 'seed-1_sample-0' /
                               'model.cif').is_file():
        return candidate

    # Look for timestamped variants
    matches = sorted(AF3_BASE.glob(f'{name}_2026*'))
    for m in matches:
        if (m / 'seed-1_sample-0' / 'model.cif').is_file():
            return m
    return None


def compute_fl_blocks(units):
    """Continuous FL residue blocks merging unit ranges (handle overlaps)."""
    ranges = sorted([UNIT_FL_RANGE[u] for u in units])
    merged = [list(ranges[0])]
    for s, e in ranges[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]


def build_fl_to_construct_map(units):
    """Return (segments, n_residues, fl_to_c).

    segments: list of (fl_start, fl_end, offset_to_construct).
    fl_to_c(fl_resid) -> construct 1-indexed residue or None.
    """
    fl_blocks = compute_fl_blocks(units)
    segments = []
    construct_pos = 1
    for fl_s, fl_e in fl_blocks:
        offset = construct_pos - fl_s
        segments.append((fl_s, fl_e, offset))
        construct_pos += (fl_e - fl_s + 1)
    n_residues = construct_pos - 1

    def fl_to_c(fl_resid):
        for fl_s, fl_e, off in segments:
            if fl_s <= fl_resid <= fl_e:
                return fl_resid + off
        return None

    return segments, n_residues, fl_to_c


def compute_construct_domains(units):
    """Map per-domain trim boundaries from FL to construct numbering.

    Returns:
        domain_ranges: list of [c_start, c_end]
        labels: list of restraint domain names (parallel to domain_ranges)
        n_residues: construct length
    """
    segments, n_residues, fl_to_c = build_fl_to_construct_map(units)

    construct_fl_set = set()
    for fl_s, fl_e, _ in segments:
        construct_fl_set.update(range(fl_s, fl_e + 1))

    unit_set = set(units)

    domain_ranges = []
    labels = []

    for rname, parent_unit in RESTRAINT_TO_UNIT.items():
        # Skip if parent unit is not in this fragment
        if parent_unit not in unit_set:
            continue

        # Compute trimmed FL boundaries
        ext_s, ext_e = STRUCTURED_EXTENTS.get(rname, DOMAIN_EXTENTS[rname])
        trim = DOMAIN_TRIMS[rname]
        fl_s = ext_s + trim
        fl_e = ext_e - trim
        if fl_e - fl_s + 1 < 10:
            continue

        # Clip to construct's FL coverage
        while fl_s <= fl_e and fl_s not in construct_fl_set:
            fl_s += 1
        while fl_e >= fl_s and fl_e not in construct_fl_set:
            fl_e -= 1
        if fl_s > fl_e:
            continue

        c_s = fl_to_c(fl_s)
        c_e = fl_to_c(fl_e)
        if c_s is None or c_e is None:
            continue
        if c_e - c_s + 1 < 10:
            continue

        domain_ranges.append([c_s, c_e])
        labels.append(rname)

    return domain_ranges, labels, n_residues


def generate_kh7a_khb_pairs_for_fragment(units, pdb_path,
                                         cutoff_nm=N_CUSTOM_RES_CUTOFF_NM):
    """Generate KH7a-KHb pairs in CONSTRUCT numbering.

    Reads CA coords from the construct's AF3 model (which has residues
    numbered 1..n_residues), and identifies pairs that correspond to
    FL KH7a (738-789) and FL KHb (1389-1461).
    """
    if 'kh7a' not in units or 'khb-kh8' not in units:
        return []

    segments, n_residues, fl_to_c = build_fl_to_construct_map(units)

    # Map KH7a (FL 738-789) and KHb (FL 1389-1461) to construct numbering
    kh7a_construct = [fl_to_c(r) for r in range(738, 790)]
    kh7a_construct = [r for r in kh7a_construct if r is not None]
    khb_construct = [fl_to_c(r) for r in range(1389, 1462)]
    khb_construct = [r for r in khb_construct if r is not None]

    if not kh7a_construct or not khb_construct:
        return []

    # Load CIF and get CA coords (residue numbering in CIF is 1..N construct)
    from Bio.PDB import MMCIFParser, PDBParser
    if str(pdb_path).endswith('.cif'):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    structure = parser.get_structure('model', str(pdb_path))

    ca = {}
    for r in structure[0].get_residues():
        if r.id[0] != ' ':
            continue
        cas = [a for a in r if a.name == 'CA']
        if cas:
            ca[r.id[1]] = cas[0].get_vector().get_array()

    cutoff_ang = cutoff_nm * 10
    pairs = []
    for r1 in kh7a_construct:
        if r1 not in ca:
            continue
        for r2 in khb_construct:
            if r2 not in ca:
                continue
            d = float(np.linalg.norm(ca[r1] - ca[r2]))
            if d <= cutoff_ang:
                pairs.append((r1, r2, d / 10.0))
    return pairs


def write_run_py(sim_dir):
    with open(sim_dir / 'run.py', 'w') as f:
        f.write("""from calvados import sim
from argparse import ArgumentParser

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--path', nargs='?', default='.', const='.', type=str)
    parser.add_argument('--config', nargs='?', default='config.yaml', const='config.yaml', type=str)
    parser.add_argument('--components', nargs='?', default='components.yaml', const='components.yaml', type=str)

    args = parser.parse_args()
    sim.run(path=args.path, fconfig=args.config, fcomponents=args.components)
""")


def convert_cif_to_pdb(cif_path, pdb_path):
    """Convert mmCIF to PDB (backbone only)."""
    from Bio.PDB import MMCIFParser, PDBIO, Select
    class BackboneSelect(Select):
        def accept_atom(self, atom):
            return atom.get_name() in ('N', 'CA', 'C', 'O')
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('m', str(cif_path))
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(pdb_path), BackboneSelect())


# ============================================================
# Prepare one fragment (all 25 replicates)
# ============================================================

def prepare_fragment(units, af3_dir, verbose=True):
    """Prepare 25 replicates for one fragment. Returns success/skip reason."""
    name = fragment_name(units)
    frag_dir = FRAGMENTS_DIR / name
    shared_input = frag_dir / 'input'
    shared_input.mkdir(parents=True, exist_ok=True)

    # Convert AF3 CIF -> PDB (use seed-1_sample-0 as the structure)
    cif_src = af3_dir / 'seed-1_sample-0' / 'model.cif'
    pdb_dst = shared_input / 'parp14.pdb'  # name must match sysname
    if not pdb_dst.exists():
        try:
            convert_cif_to_pdb(cif_src, pdb_dst)
        except Exception as e:
            return False, f'CIF conversion failed: {e}'

    # Symlink residues file
    res_dst = shared_input / 'residues_CALVADOS3.csv'
    if not res_dst.exists():
        res_dst.symlink_to(FL_RESIDUES.resolve())

    # Compute construct domains
    domain_ranges, labels, n_residues = compute_construct_domains(units)

    # Write domains.yaml
    with open(shared_input / 'domains.yaml', 'w') as f:
        yaml.dump({'parp14': domain_ranges}, f, default_flow_style=True)

    # Generate KH7a-KHb custom restraints if applicable
    custom_pairs = generate_kh7a_khb_pairs_for_fragment(units, pdb_dst)
    use_custom = len(custom_pairs) > 0
    if use_custom:
        with open(shared_input / 'custom_restraints.txt', 'w') as f:
            for r1, r2, d_nm in sorted(custom_pairs):
                f.write(f'parp14 1 {r1} | parp14 1 {r2} | '
                        f'{d_nm:.3f} {K_CUSTOM:.1f}\n')

    # Choose box size: larger for larger constructs
    box_nm = max(40, int(n_residues * 0.15) + 50)
    box_nm = min(box_nm, 300)  # cap

    # Save metadata
    meta = {
        'fragment': name,
        'units': units,
        'n_residues': n_residues,
        'n_restraint_domains': len(domain_ranges),
        'n_restrained_residues': sum(de - ds + 1 for ds, de in domain_ranges),
        'restraint_labels': labels,
        'domain_ranges_construct': [list(r) for r in domain_ranges],
        'n_custom_pairs': len(custom_pairs),
        'use_custom_kh7a_khb': use_custom,
        'box_nm': box_nm,
        'af3_source': str(cif_src),
    }
    with open(frag_dir / 'metadata.json', 'w') as f:
        json.dump(meta, f, indent=2)

    # Load CALVADOS default config
    with open(DEFAULT_CONFIG) as f:
        base_config = yaml.safe_load(f)

    # Create 25 replicate dirs
    for seed in SEEDS:
        for sample in SAMPLES:
            rep_name = f'seed-{seed}_sample-{sample}'
            sim_dir = frag_dir / rep_name
            input_dir = sim_dir / 'input'
            input_dir.mkdir(parents=True, exist_ok=True)

            # Symlink shared inputs
            for fname in ['parp14.pdb', 'residues_CALVADOS3.csv',
                          'domains.yaml']:
                dst = input_dir / fname
                src = shared_input / fname
                if not dst.exists() and src.exists():
                    dst.symlink_to(src.resolve())
            if use_custom:
                dst = input_dir / 'custom_restraints.txt'
                src = shared_input / 'custom_restraints.txt'
                if not dst.exists() and src.exists():
                    dst.symlink_to(src.resolve())

            # Per-replicate config
            random_seed = seed * 1000 + sample
            cfg = dict(base_config)
            cfg.update({
                'sysname': 'parp14',
                'box': [box_nm, box_nm, box_nm],
                'temp': TEMP, 'ionic': IONIC, 'pH': PH,
                'topol': 'center',
                'steps': N_STEPS, 'wfreq': WFREQ, 'runtime': 0,
                'platform': PLATFORM, 'threads': THREADS,
                'restart': 'checkpoint', 'frestart': 'restart.chk',
                'verbose': True,
                'random_number_seed': random_seed,
                'custom_restraints': use_custom,
                'custom_restraint_type': 'harmonic',
                'fcustom_restraints':
                    'input/custom_restraints.txt' if use_custom else
                    'custom_restraints.txt',
            })
            with open(sim_dir / 'config.yaml', 'w') as f:
                yaml.dump(cfg, f, default_flow_style=False)

            # components.yaml
            components = {
                'defaults': {
                    'molecule_type': 'protein', 'nmol': 1,
                    'charge_termini': 'both', 'alpha': 0,
                    'ffasta': 'fastabib.fasta', 'kb': 8033.0,
                    'ext_restraint': True,
                    'restraint': True,
                    'cutoff_restr': CUTOFF_RESTR,
                    'pdb_folder': str(input_dir),
                    'restraint_type': 'harmonic',
                    'k_harmonic': K_HARMONIC,
                    'fdomains': str(input_dir / 'domains.yaml'),
                    'k_go': 15.0, 'use_com': True, 'periodic': False,
                    'colabfold': 2,   # AF3 colabfold mode
                    'bfac_shift': 0.8, 'bfac_width': 50.0,
                    'pae_shift': 0.3, 'pae_width': 15.0,
                    'rna_kb1': 8033.0, 'rna_kb2': 8033.0,
                    'rna_ka': 7.24, 'rna_pa': 3.14,
                    'rna_nb_sigma': 0.4, 'rna_nb_scale': 15,
                    'rna_nb_cutoff': 0.6, 'n_ends': 1,
                    'ptm_name': 'example_ptm', 'ptm_locations': [],
                    'fresidues': str(FL_RESIDUES),
                },
                'system': {'parp14': {}},
            }
            with open(sim_dir / 'components.yaml', 'w') as f:
                yaml.dump(components, f, default_flow_style=False)

            write_run_py(sim_dir)

    return True, meta


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--list', action='store_true',
                        help='List all 66 fragments with AF3 status')
    parser.add_argument('--fragments', nargs='+',
                        help='Only prepare specific fragments (by canonical '
                             'name, e.g. md1l1_md2)')
    parser.add_argument('--skip-existing', action='store_true',
                        help='Skip fragments already prepared')
    args = parser.parse_args()

    fragments = all_contiguous_fragments()
    print(f"All contiguous fragments: {len(fragments)}")

    # Check AF3 availability
    available = []
    missing = []
    for units in fragments:
        name = fragment_name(units)
        af3 = af3_dir_for_fragment(units)
        if af3:
            available.append((units, af3))
        else:
            # Special case: full-length uses AF2 PDB
            if units == UNIT_ORDER:
                available.append((units, 'AF2'))
            else:
                missing.append(name)

    print(f"  Available (AF3 or AF2):  {len(available)}")
    print(f"  Missing structure:       {len(missing)}")

    if args.list:
        print(f"\n{'#':>2} {'Fragment':<60} {'AF3 status':<20} {'nUnits':>6}")
        print("-" * 95)
        for i, units in enumerate(fragments, 1):
            name = fragment_name(units)
            af3 = af3_dir_for_fragment(units)
            if af3:
                status = af3.name
            elif units == UNIT_ORDER:
                status = 'AF2 (full-length)'
            else:
                status = '--- missing ---'
            print(f"{i:>2} {name:<60} {status:<20} {len(units):>6}")
        if missing:
            print(f"\n{len(missing)} fragments missing AF3 structures:")
            for n in missing:
                print(f"  {n}")
        return

    # Filter to specific fragments
    if args.fragments:
        wanted = set(args.fragments)
        available = [(u, a) for u, a in available
                     if fragment_name(u) in wanted]
        print(f"  Filtered to {len(available)} fragments")

    FRAGMENTS_DIR.mkdir(exist_ok=True)
    print(f"\nPreparing into {FRAGMENTS_DIR}/")

    success = []
    failed = []
    skipped = []

    for units, af3 in available:
        name = fragment_name(units)
        frag_dir = FRAGMENTS_DIR / name

        if args.skip_existing and (frag_dir / 'metadata.json').exists():
            skipped.append(name)
            continue

        # Use AF2 PDB for full-length
        if af3 == 'AF2':
            # Special handling: copy AF2 PDB instead of converting AF3 CIF
            frag_dir.mkdir(exist_ok=True)
            shared = frag_dir / 'input'
            shared.mkdir(exist_ok=True)
            pdb_dst = shared / 'parp14.pdb'
            if not pdb_dst.exists():
                shutil.copy2(FL_PDB, pdb_dst)
            # Create a fake af3_dir pointer so prepare_fragment can find the PDB
            # We'll directly call prepare with af3 as None and handle FL case
            # Simpler: skip if fl_optimized already exists
            print(f"  [SKIP] {name}: use fl_optimized/ for full-length")
            skipped.append(name)
            continue

        try:
            ok, msg = prepare_fragment(units, af3)
            if ok:
                success.append((name, msg))
                print(f"  [OK] {name}: {msg['n_residues']} res, "
                      f"{msg['n_restraint_domains']} domains, "
                      f"{msg['n_custom_pairs']} custom pairs")
            else:
                failed.append((name, msg))
                print(f"  [FAIL] {name}: {msg}")
        except Exception as e:
            failed.append((name, str(e)))
            print(f"  [ERROR] {name}: {e}")

    # Write launch commands file
    cmd_file = FRAGMENTS_DIR / 'launch_commands.txt'
    with open(cmd_file, 'w') as f:
        f.write("# Launch commands for PARP14 fragment simulations\n")
        f.write("# Each line launches all 25 replicates of one fragment in\n")
        f.write("# parallel (8 jobs). Uncomment lines you want to run.\n")
        f.write(f"# Python env: /home/sbali/miniconda3/envs/CALVADOS/bin/python\n")
        f.write(f"#\n")
        f.write(f"# Total: {len(success)} fragments prepared\n")
        f.write(f"# Each: 25 replicates x {N_STEPS * 0.01 / 1000:.0f} ns = "
                f"{25 * N_STEPS * 0.01 / 1000:.0f} ns total\n\n")

        py = "/home/sbali/miniconda3/envs/CALVADOS/bin/python"
        for name, meta in success:
            f.write(f"# {name} ({meta['n_residues']} res, "
                    f"{meta['n_restraint_domains']} domains, "
                    f"{meta['n_custom_pairs']} custom pairs)\n")
            f.write(f"#cd {FRAGMENTS_DIR}/{name} && "
                    f"ls -d seed-*_sample-* | "
                    f"xargs -I{{}} -P 8 bash -c "
                    f"'cd {{}} && {py} run.py > run.log 2>&1'\n\n")

    print(f"\n{'='*60}")
    print(f"Prepared: {len(success)}  Skipped: {len(skipped)}  "
          f"Failed: {len(failed)}")
    print(f"{'='*60}")
    print(f"Launch commands: {cmd_file}")
    print(f"  (uncomment the lines you want to run)")


if __name__ == '__main__':
    main()
