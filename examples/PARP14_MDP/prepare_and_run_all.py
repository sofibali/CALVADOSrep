#!/usr/bin/env python3
"""
Prepare and launch all PARP14 CALVADOS simulations.

Six simulation sets, 25 replicates each:
  1. fl     — Full-length (1801 res), EBI AF2 structure, 20 ns
  2. md     — Macrodomains only: MD1L1+MD2+MD3 (586 res), 20 ns
  3. core   — Core: KH7a+MD1L1+MD2+MD3+KHb-KH8+WWE+ART (1051 res), 20 ns
  4. mka    — MD1L1+MD2+MD3+KHb-KH8+WWE+ART (999 res), 20 ns
  5. norrm  — No-RRM: KH1-6 through ART (1474 res), 20 ns
  6. noart  — No-ART: KH1-6 through WWE (1275 res), 20 ns

Directory layout (nested):
  {set}/input/                  # shared inputs (domains.yaml, residues)
  {set}/seed-{N}_sample-{M}/   # per-replicate simulation directory

Domain restraint boundaries use the FL PARP14 simulation domains.yaml
(structured cores only), mapped to construct numbering. Linker regions
between domains are left flexible (unrestrained).

Usage:
    python prepare_and_run_all.py --prepare              # prepare all inputs
    python prepare_and_run_all.py --prepare --set norrm  # prepare one set only
    python prepare_and_run_all.py --write-launcher       # write bash launcher
    python prepare_and_run_all.py --run                  # run all sequentially
"""

import os
import sys
import shutil
import yaml
import subprocess
from argparse import ArgumentParser
from Bio.PDB import MMCIFParser, PDBIO, Select
from calvados.cfg import Config, Components

# ============================================================
# Configuration
# ============================================================

CWD = os.path.dirname(os.path.abspath(__file__))
AF3_BASE = os.path.join('/home/sbali/CALVADOS/parp14/alphafold_outputs')

# Simulation parameters (shared)
N_STEPS = 500_000       # 5 ns per replicate
WFREQ = 1_000           # -> 500 frames per sim
TEMP = 293              # K
IONIC = 0.19            # M
PH = 7.0
THREADS = 4
PLATFORM = 'CPU'
K_HARMONIC = 700.0      # kJ/mol/nm^2

SEEDS = range(1, 6)
SAMPLES = range(0, 5)

# Shared residues file
FL_RESIDUES = os.path.join(CWD, 'input', 'residues_CALVADOS3.csv')

# ============================================================
# FL restraint domain boundaries (from input/domains.yaml)
# These define only the structured cores; linkers are left flexible.
# ============================================================

# Each restraint domain is mapped to its primary domain unit.
# A restraint domain is only included in a construct if its primary unit
# is in the construct's unit list (prevents edge-overlap artifacts).
FL_RESTRAINT_DOMAINS = [
    # (name, fl_start, fl_end, primary_unit)
    ('RRM1', 6, 88,     'rrm1'),
    ('RRM2', 150, 223,  'rrm2'),
    ('RRM3', 227, 301,  'rrm3'),
    ('MD1',  791, 978,  'md1l1'),
    ('MD2',  1003, 1190, 'md2'),
    ('MD3',  1216, 1387, 'md3'),
    ('WWE',  1523, 1601, 'wwe'),
    ('ART',  1605, 1801, 'art'),
]

# ============================================================
# Domain units used in AF3 combinatorial library
# These define which FL residues are included in each construct.
# Overlapping boundaries (e.g. MD1L1 ends at 1004, MD2 starts at 1004)
# are merged into continuous blocks.
# ============================================================

DOMAIN_UNITS = {
    'rrm1':     (1, 145),
    'rrm2':     (146, 224),
    'rrm3':     (225, 314),
    'kh1-kh6':  (315, 737),
    'kh7a':     (738, 789),
    'md1l1':    (790, 1004),
    'md2':      (1004, 1193),
    'md3':      (1207, 1388),
    'khb-kh8':  (1389, 1533),
    'wwe':      (1534, 1602),
    'art':      (1603, 1801),
}

# ============================================================
# Construct definitions
# ============================================================

CONSTRUCTS = {
    'md': {
        'label': 'Macrodomains only (MD1L1+MD2+MD3)',
        'units': ['md1l1', 'md2', 'md3'],
        'af3_dir': 'md1l1_md2_md3',
        'sysname': 'parp14_macrodomains',
        'box': 80,
    },
    'core': {
        'label': 'Core (KH7a+MD1L1+MD2+MD3+KHb-KH8+WWE+ART)',
        'units': ['kh7a', 'md1l1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art'],
        'af3_dir': 'kh7a_md1l1_md2_md3_khb-kh8_wwe_art',
        'sysname': 'parp14_core',
        'box': 120,
    },
    'mka': {
        'label': 'MD1L1+MD2+MD3+KHb-KH8+WWE+ART',
        'units': ['md1l1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art'],
        'af3_dir': 'md1l1_md2_md3_khb-kh8_wwe_art',
        'sysname': 'parp14_mka',
        'box': 100,
    },
    'norrm': {
        'label': 'No-RRM (KH1-6+KH7a+MD1L1+MD2+MD3+KHb-KH8+WWE+ART)',
        'units': ['kh1-kh6', 'kh7a', 'md1l1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art'],
        'af3_dir': 'kh1-kh6_kh7a_md1l1_md2_md3_khb-kh8_wwe_art',
        'sysname': 'parp14_norrm',
        'box': 250,
    },
    'noart': {
        'label': 'No-ART (KH1-6+KH7a+MD1L1+MD2+MD3+KHb-KH8+WWE)',
        'units': ['kh1-kh6', 'kh7a', 'md1l1', 'md2', 'md3', 'khb-kh8', 'wwe'],
        'af3_dir': 'kh1-kh6_kh7a_md1l1_md2_md3_khb-kh8_wwe',
        'sysname': 'parp14_noart',
        'box': 200,
    },
    'md3art': {
        'label': 'MD3-ART (MD3+KHb-KH8+WWE+ART)',
        'units': ['md3', 'khb-kh8', 'wwe', 'art'],
        'af3_dir': 'md3_khb-kh8_wwe_art_20260320_175029',
        'sysname': 'parp14_md3art',
        'box': 80,
    },
}


# ============================================================
# FL-to-construct domain boundary mapping
# ============================================================

def compute_fl_blocks(unit_names):
    """
    Given a list of domain unit names, compute the continuous FL residue
    blocks (merging overlaps), and return the blocks as sorted
    [(fl_start, fl_end), ...].
    """
    # Collect all FL residue ranges
    ranges = []
    for name in unit_names:
        s, e = DOMAIN_UNITS[name]
        ranges.append((s, e))
    ranges.sort()

    # Merge overlapping/adjacent ranges
    merged = [list(ranges[0])]
    for s, e in ranges[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])

    return [(s, e) for s, e in merged]


def build_fl_to_construct_map(fl_blocks):
    """
    Build a mapping function from FL residue number to construct residue number.
    fl_blocks: sorted list of (fl_start, fl_end) continuous blocks.
    Returns: function(fl_resid) -> construct_resid or None
    """
    segments = []  # (fl_start, fl_end, offset_to_add)
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


def compute_construct_domains(unit_names):
    """
    Map FL restraint domain boundaries to construct numbering.
    Only includes domains whose FL residue range overlaps with the construct.
    Clips boundaries to the construct's residue range.
    Returns: list of [construct_start, construct_end] for domains.yaml
    """
    fl_blocks = compute_fl_blocks(unit_names)
    fl_to_c = build_fl_to_construct_map(fl_blocks)

    # Set of all FL residues in the construct
    construct_fl_residues = set()
    for s, e in fl_blocks:
        construct_fl_residues.update(range(s, e + 1))

    unit_set = set(unit_names)

    domains = []
    domain_labels = []
    for name, fl_start, fl_end, primary_unit in FL_RESTRAINT_DOMAINS:
        # Skip if the primary domain unit is not in this construct
        if primary_unit not in unit_set:
            continue
        # Clip to construct range
        clipped_start = fl_start
        clipped_end = fl_end
        while clipped_start <= fl_end and clipped_start not in construct_fl_residues:
            clipped_start += 1
        while clipped_end >= clipped_start and clipped_end not in construct_fl_residues:
            clipped_end -= 1

        if clipped_start > clipped_end:
            continue  # domain not in this construct

        c_start = fl_to_c(clipped_start)
        c_end = fl_to_c(clipped_end)
        if c_start is not None and c_end is not None:
            domains.append([c_start, c_end])
            domain_labels.append(name)

    return domains, domain_labels


# ============================================================
# Helpers
# ============================================================

class BackboneSelect(Select):
    """Select only backbone atoms for PDB output."""
    def accept_atom(self, atom):
        return atom.get_name() in ['N', 'CA', 'C', 'O']


def convert_cif_to_pdb(cif_path, pdb_path):
    """Convert mmCIF to PDB format (backbone only)."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('model', cif_path)
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_path, BackboneSelect())


def create_run_script(sim_dir):
    """Write the standard CALVADOS run.py script."""
    run_py = os.path.join(sim_dir, 'run.py')
    with open(run_py, 'w') as f:
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


# ============================================================
# Prepare full-length replicates
# ============================================================

def prepare_fl(n_steps=N_STEPS):
    """Prepare 25 full-length replicates (same structure, different random seeds)."""
    print("\n" + "=" * 70)
    print(f"Preparing: Full-Length PARP14 (1801 res, 25 x {n_steps * 0.01 / 1000:.0f} ns)")
    print("=" * 70)

    set_dir = os.path.join(CWD, 'fl')
    os.makedirs(set_dir, exist_ok=True)

    # metadata.json so analysis scripts can use `--sim-folder fl` (full-length
    # contains all 11 domain units). See sim_registry.py / stamp_metadata.py.
    try:
        import json as _json
        with open(os.path.join(set_dir, 'metadata.json'), 'w') as f:
            _json.dump({'units': list(DOMAIN_UNITS.keys()), 'sysname': 'parp14',
                        'sites': ['MD1', 'MD2', 'MD3', 'ART'],
                        'label': 'Full Length (1-1801)'}, f, indent=2)
    except Exception as _e:
        print(f"  WARN: could not write metadata.json ({_e})")

    # Shared input — use fl/input/
    fl_input = os.path.join(set_dir, 'input')
    os.makedirs(fl_input, exist_ok=True)

    # Ensure shared files are in fl/input/
    src_input = os.path.join(CWD, 'input')
    for fname in ['parp14.pdb', 'domains.yaml', 'residues_CALVADOS3.csv']:
        dest = os.path.join(fl_input, fname)
        src = os.path.join(src_input, fname)
        if not os.path.isfile(dest) and os.path.isfile(src):
            shutil.copy2(src, dest)

    fl_pdb = os.path.join(fl_input, 'parp14.pdb')
    fl_domains = os.path.join(fl_input, 'domains.yaml')
    fl_residues = os.path.join(fl_input, 'residues_CALVADOS3.csv')

    n_prepared = 0
    for seed in SEEDS:
        for sample in SAMPLES:
            rep_name = f'seed-{seed}_sample-{sample}'
            sim_dir = os.path.join(set_dir, rep_name)
            input_dir = os.path.join(sim_dir, 'input')

            print(f"  Preparing: fl/{rep_name}")
            os.makedirs(input_dir, exist_ok=True)

            # Copy FL PDB
            pdb_dest = os.path.join(input_dir, 'parp14.pdb')
            if not os.path.isfile(pdb_dest):
                shutil.copy2(fl_pdb, pdb_dest)

            # Config with unique random seed
            random_seed = seed * 1000 + sample
            config = Config(
                sysname='parp14',
                box=[300, 300, 300],
                temp=TEMP, ionic=IONIC, pH=PH,
                topol='center',
                wfreq=WFREQ, steps=n_steps, runtime=0,
                platform=PLATFORM, threads=THREADS,
                restart='checkpoint', frestart='restart.chk',
                verbose=True,
                random_number_seed=random_seed,
            )
            config.write(sim_dir, name='config.yaml')

            components = Components(
                molecule_type='protein', nmol=1,
                restraint=True, charge_termini='both',
                fresidues=fl_residues,
                fdomains=fl_domains,
                pdb_folder=input_dir,
                restraint_type='harmonic', use_com=True,
                colabfold=1, k_harmonic=K_HARMONIC,
            )
            components.add(name='parp14')
            components.write(sim_dir, name='components.yaml')

            create_run_script(sim_dir)
            n_prepared += 1

    print(f"\n  Prepared {n_prepared} full-length replicates")
    return n_prepared


# ============================================================
# Prepare AF3 construct replicates (generic)
# ============================================================

def prepare_construct(key, n_steps=N_STEPS):
    """Prepare 25 replicates for a given construct."""
    info = CONSTRUCTS[key]
    sysname = info['sysname']
    af3_dir = os.path.join(AF3_BASE, info['af3_dir'])
    box = info['box']

    # Compute domain boundaries from FL mapping
    fl_blocks = compute_fl_blocks(info['units'])
    n_residues = sum(e - s + 1 for s, e in fl_blocks)
    domains, domain_labels = compute_construct_domains(info['units'])

    ns = n_steps * 0.01 / 1000
    print(f"\n{'='*70}")
    print(f"Preparing: {info['label']} ({n_residues} res, 25 x {ns:.0f} ns)")
    print(f"{'='*70}")

    # Create nested set dir and shared input dir
    set_dir = os.path.join(CWD, key)
    shared_input = os.path.join(set_dir, 'input')
    os.makedirs(shared_input, exist_ok=True)

    # Write domains.yaml
    domains_yaml = os.path.join(shared_input, 'domains.yaml')
    with open(domains_yaml, 'w') as f:
        yaml.dump({sysname: domains}, f, default_flow_style=False, sort_keys=False)

    # Write metadata.json so the analysis scripts can analyze this set via
    # `--sim-folder <set>` with no other arguments (units -> active-site/domain
    # remapping). See sim_registry.py / stamp_metadata.py.
    try:
        import json as _json
        _meta = {
            'units': list(info['units']),
            'sysname': sysname,
            'sites': [s for u, s in [('md1l1', 'MD1'), ('md2', 'MD2'),
                                     ('md3', 'MD3'), ('art', 'ART')]
                      if u in info['units']],
            'label': info.get('label', key),
            'domain_ranges_construct': domains,
        }
        with open(os.path.join(set_dir, 'metadata.json'), 'w') as f:
            _json.dump(_meta, f, indent=2)
    except Exception as _e:
        print(f"  WARN: could not write metadata.json ({_e})")

    print(f"  Domains (FL boundaries -> construct numbering):")
    for label, bounds in zip(domain_labels, domains):
        print(f"    {label}: {bounds}")

    # Compute unrestrained (flexible) regions
    restrained = set()
    for s, e in domains:
        restrained.update(range(s, e + 1))
    flexible = set(range(1, n_residues + 1)) - restrained
    if flexible:
        flex_sorted = sorted(flexible)
        ranges = []
        start = flex_sorted[0]
        prev = start
        for r in flex_sorted[1:]:
            if r == prev + 1:
                prev = r
            else:
                ranges.append(f"{start}-{prev}" if start != prev else str(start))
                start = r
                prev = r
        ranges.append(f"{start}-{prev}" if start != prev else str(start))
        print(f"  Flexible (unrestrained) regions: {', '.join(ranges)}")

    # Copy residues file
    residues_dest = os.path.join(shared_input, 'residues_CALVADOS3.csv')
    if not os.path.isfile(residues_dest):
        shutil.copy2(FL_RESIDUES, residues_dest)

    n_prepared = 0
    n_skipped = 0

    for seed in SEEDS:
        for sample in SAMPLES:
            seed_sample = f'seed-{seed}_sample-{sample}'
            rep_name = seed_sample
            sim_dir = os.path.join(set_dir, rep_name)
            input_dir = os.path.join(sim_dir, 'input')

            cif_path = os.path.join(af3_dir, seed_sample, 'model.cif')
            conf_json = os.path.join(af3_dir, seed_sample, 'confidences.json')

            if not os.path.isfile(cif_path):
                print(f"  SKIP: {seed_sample} - model.cif not found")
                n_skipped += 1
                continue

            print(f"  Preparing: {key}/{rep_name}")
            os.makedirs(input_dir, exist_ok=True)

            # Convert CIF -> PDB
            pdb_path = os.path.join(input_dir, f'{sysname}.pdb')
            convert_cif_to_pdb(cif_path, pdb_path)

            # Copy PAE JSON
            pae_path = os.path.join(input_dir, f'{sysname}.json')
            shutil.copy2(conf_json, pae_path)

            # Config
            config = Config(
                sysname=sysname,
                box=[box, box, box],
                temp=TEMP, ionic=IONIC, pH=PH,
                topol='center',
                wfreq=WFREQ, steps=n_steps, runtime=0,
                platform=PLATFORM, threads=THREADS,
                restart='checkpoint', frestart='restart.chk',
                verbose=True,
            )
            config.write(sim_dir, name='config.yaml')

            # Components
            components = Components(
                molecule_type='protein', nmol=1,
                restraint=True, charge_termini='both',
                fresidues=os.path.join(shared_input, 'residues_CALVADOS3.csv'),
                fdomains=domains_yaml,
                pdb_folder=input_dir,
                restraint_type='harmonic', use_com=True,
                colabfold=2, k_harmonic=K_HARMONIC,
            )
            components.add(name=sysname)
            components.write(sim_dir, name='components.yaml')

            create_run_script(sim_dir)
            n_prepared += 1

    print(f"\n  Prepared {n_prepared} replicates (skipped {n_skipped})")
    return n_prepared


# ============================================================
# Launcher script
# ============================================================

def write_launcher():
    """Write a bash script that runs all simulations."""
    all_sets = ['fl', 'md', 'core', 'mka', 'norrm', 'noart', 'md3art']
    launcher_path = os.path.join(CWD, 'run_all.sh')
    with open(launcher_path, 'w') as f:
        f.write(f"""#!/bin/bash
# Run all PARP14 CALVADOS simulations
# 4 sets x 25 replicates x 5 ns = 500 ns total (125 ns per set)
#
# Usage:
#   bash run_all.sh              # run all sequentially
#   bash run_all.sh parallel     # run all in parallel (background)
#   bash run_all.sh fl           # run only full-length
#   bash run_all.sh md           # run only macrodomains
#   bash run_all.sh core         # run only core construct
#   bash run_all.sh mka          # run only MD+KH+ART construct

set -e
cd "$(dirname "$0")"

MODE="${{1:-sequential}}"
PIDS=()

run_sim() {{
    local dir="$1"
    if [ ! -f "$dir/config.yaml" ]; then
        echo "SKIP: $dir (no config.yaml)"
        return
    fi
    if [ -f "$dir/restart.chk" ]; then
        echo "RUN (restart): $dir"
    else
        echo "RUN (fresh):   $dir"
    fi

    if [ "$MODE" = "parallel" ]; then
        (cd "$dir" && python run.py) &
        PIDS+=($!)
    else
        (cd "$dir" && python run.py)
    fi
}}

run_set() {{
    local prefix="$1"
    local label="$2"
    echo "=================================================="
    echo "$label (25 replicates x 5 ns = 125 ns)"
    echo "=================================================="
    for seed in 1 2 3 4 5; do
        for sample in 0 1 2 3 4; do
            run_sim "${{prefix}}/seed-${{seed}}_sample-${{sample}}"
        done
    done
}}

# Determine which sets to run
case "$MODE" in
    fl)    run_set "fl" "Full-Length PARP14 (1801 res)" ;;
    md)    run_set "md" "Macrodomains only (586 res)" ;;
    core)  run_set "core" "Core construct (1051 res)" ;;
    mka)   run_set "mka" "MD+KHb+ART construct (930 res)" ;;
    norrm) run_set "norrm" "No-RRM (1474 res)" ;;
    noart)  run_set "noart" "No-ART (1275 res)" ;;
    md3art) run_set "md3art" "MD3-ART (595 res)" ;;
    *)
        run_set "fl" "Full-Length PARP14 (1801 res)"
        run_set "md" "Macrodomains only (586 res)"
        run_set "core" "Core construct (1051 res)"
        run_set "mka" "MD+KHb+ART construct (930 res)"
        run_set "norrm" "No-RRM (1474 res)"
        run_set "noart" "No-ART (1275 res)"
        run_set "md3art" "MD3-ART (595 res)"
        ;;
esac

# Wait for parallel jobs
if [ "$MODE" = "parallel" ] && [ ${{#PIDS[@]}} -gt 0 ]; then
    echo ""
    echo "Waiting for ${{#PIDS[@]}} parallel jobs..."
    for pid in "${{PIDS[@]}}"; do
        wait "$pid"
    done
fi

echo ""
echo "=================================================="
echo "SIMULATIONS COMPLETE"
echo "=================================================="
echo "Protocol: 25 replicates x 5 ns = 125 ns per set"
echo "Equilibration: discard first 0.5 ns of each replicate"
echo "Effective: 25 x 4.5 ns = 112.5 ns compiled per set"
""")

    os.chmod(launcher_path, 0o755)
    print(f"\n  Wrote: {launcher_path}")
    print(f"  Usage:")
    print(f"    bash run_all.sh              # all sets, sequential")
    print(f"    bash run_all.sh parallel     # all sets, parallel")
    print(f"    bash run_all.sh fl           # full-length only")
    print(f"    bash run_all.sh md           # macrodomains only")
    print(f"    bash run_all.sh core         # core construct only")
    print(f"    bash run_all.sh mka          # MD+KHb+ART only")


# ============================================================
# Main
# ============================================================

def main():
    parser = ArgumentParser(description='Prepare and run all PARP14 simulations')
    parser.add_argument('--prepare', action='store_true', help='Prepare simulation inputs')
    parser.add_argument('--run', action='store_true', help='Run all simulations sequentially')
    parser.add_argument('--write-launcher', action='store_true', help='Write bash launcher script')
    parser.add_argument('--set', type=str, default=None,
                        choices=['fl', 'md', 'core', 'mka', 'norrm', 'noart', 'md3art'],
                        help='Only prepare/run a specific set')
    args = parser.parse_args()

    if not any([args.prepare, args.run, args.write_launcher]):
        parser.print_help()
        print("\n  Specify at least one of: --prepare, --run, --write-launcher")
        sys.exit(1)

    all_sets = ['fl', 'md', 'core', 'mka', 'norrm', 'noart', 'md3art']
    active_sets = [args.set] if args.set else all_sets

    print("=" * 70)
    print("PARP14 CALVADOS Simulations")
    print("=" * 70)
    print(f"  Per replicate:  {N_STEPS:,} steps = {N_STEPS * 0.01 / 1000:.1f} ns")
    print(f"  Replicates:     25 per set (5 seeds x 5 samples)")
    print(f"  Total per set:  125 ns")
    print(f"  Equilibration:  discard first 0.5 ns (50 frames)")
    print(f"  Sets:           {', '.join(active_sets)}")

    if args.prepare:
        total = 0
        for key in active_sets:
            if key == 'fl':
                total += prepare_fl()
            else:
                total += prepare_construct(key)
        print(f"\n  TOTAL PREPARED: {total} simulations")

    if args.write_launcher:
        write_launcher()

    if args.run:
        print("\n--- Running simulations sequentially ---")
        for key in active_sets:
            for seed in SEEDS:
                for sample in SAMPLES:
                    rep_name = f'seed-{seed}_sample-{sample}'
                    sim_dir = os.path.join(CWD, key, rep_name)
                    if os.path.isfile(os.path.join(sim_dir, 'config.yaml')):
                        print(f"\n  Running: {key}/{rep_name}")
                        subprocess.run([sys.executable, 'run.py'],
                                       cwd=sim_dir, check=True)

    print("\n" + "=" * 70)
    print("DONE")
    print("=" * 70)


if __name__ == '__main__':
    main()
