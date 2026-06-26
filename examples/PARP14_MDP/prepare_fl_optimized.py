#!/usr/bin/env python3
"""
Prepare 25-replicate full-length PARP14 CALVADOS simulations with the
custom per-domain trim values chosen from the restraint_tests sweep.

Output: fl_optimized/seed-{1..5}_sample-{0..4}/ (25 total)
Shared inputs: fl_optimized/input/

Per-domain trim values (residues removed from each side of starting extent):
    RRM1   trim 10  (6-88   -> 16-78)
    RRM2   trim 10
    RRM3   trim 10
    KH1    trim 5
    KH2    trim 5
    KH3    trim 5
    KH4    trim 10
    KH5    trim 5
    KH6    trim 10
    KH7a   trim 0   (full extent — needed for KH7a-KHb custom restraints)
    MD1L1  trim 10
    MD2    trim 10
    MD3    trim 10
    KHb    trim 0   (full extent — needed for KH7a-KHb custom restraints)
    KH8    trim 0
    WWE    trim 15
    ART    trim 10

KH7a-KHb custom inter-domain restraints (k=350) generated from the AF2
structure at trim 0 boundaries (full KH7a 738-789 vs full KHb 1389-1461),
141 pairs within 0.9 nm.

Usage:
    python prepare_fl_optimized.py            # prepare 25 replicate dirs
    bash run_fl_optimized.sh                  # launch all 25 (parallel)
"""

import os
import sys
import json
import yaml
import shutil
import numpy as np
from pathlib import Path

# ============================================================
# Configuration
# ============================================================

CWD = Path(__file__).resolve().parent
SET_DIR = CWD / 'fl_optimized'
SHARED_INPUT = SET_DIR / 'input'

FL_PDB = CWD / 'input' / 'parp14.pdb'
FL_RESIDUES = CWD / 'input' / 'residues_CALVADOS3.csv'
DEFAULT_CONFIG = CWD.parent.parent / 'calvados' / 'data' / 'default_config.yaml'

# Simulation parameters (match other fl sims)
N_STEPS = 2_000_000      # 20 ns at dt=0.01 ps (matches fl/ set)
WFREQ = 1_000            # 2000 frames
TEMP = 293
IONIC = 0.19
PH = 7.0
PLATFORM = 'CPU'
THREADS = 4
K_HARMONIC = 700.0
K_CUSTOM = 350.0
CUTOFF_RESTR = 0.9

SEEDS = range(1, 6)
SAMPLES = range(0, 5)

# ============================================================
# Domain extents (where to start trimming from)
# ============================================================

DOMAIN_EXTENTS = {
    'RRM1': (1, 145), 'RRM2': (146, 224), 'RRM3': (225, 314),
    'KH1': (315, 384), 'KH2': (385, 454), 'KH3': (455, 520),
    'KH4': (521, 593), 'KH5': (594, 665), 'KH6': (666, 737),
    'KH7a': (738, 789),
    'MD1L1': (790, 978),   # 3VFQ boundary
    'MD2': (1005, 1193),   # 3VFQ boundary
    'MD3': (1207, 1388),
    'KHb': (1389, 1461), 'KH8': (1462, 1533),
    'WWE': (1534, 1602), 'ART': (1603, 1801),
}

# RRM1 has a custom structured core (excludes disordered N-term)
STRUCTURED_EXTENTS = {
    'RRM1': (6, 88),
}

# Per-domain trim values (user-specified)
DOMAIN_TRIMS = {
    'RRM1':  10,
    'RRM2':  10,
    'RRM3':  10,
    'KH1':    5,
    'KH2':    5,
    'KH3':    5,
    'KH4':   10,
    'KH5':    5,
    'KH6':   10,
    'KH7a':   0,   # full extent for KH7a-KHb interface
    'MD1L1': 10,
    'MD2':   10,
    'MD3':   10,
    'KHb':    0,   # full extent for KH7a-KHb interface
    'KH8':    0,
    'WWE':   15,
    'ART':   10,
}


def compute_domain_boundaries():
    """Apply each domain's trim value to its starting extent."""
    result = {}
    for dname, trim in DOMAIN_TRIMS.items():
        ext_s, ext_e = STRUCTURED_EXTENTS.get(dname, DOMAIN_EXTENTS[dname])
        new_s = ext_s + trim
        new_e = ext_e - trim
        if new_e <= new_s:
            print(f"  WARNING: {dname} trim={trim} produces empty range, skipping")
            continue
        result[dname] = (new_s, new_e)
    return result


def generate_kh7a_khb_pairs(pdb_path, kh7a_range=(738, 789),
                            khb_range=(1389, 1461), cutoff_nm=0.9):
    """Generate KH7a-KHb inter-domain restraint pairs from PDB at trim 0
    boundaries (full extents)."""
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('af2', str(pdb_path))
    ca_by_resid = {}
    for r in structure[0].get_residues():
        if r.id[0] != ' ':
            continue
        cas = [a for a in r if a.name == 'CA']
        if cas:
            ca_by_resid[r.id[1]] = cas[0].get_vector().get_array()

    cutoff_ang = cutoff_nm * 10
    pairs = []
    for r1 in range(kh7a_range[0], kh7a_range[1] + 1):
        if r1 not in ca_by_resid:
            continue
        for r2 in range(khb_range[0], khb_range[1] + 1):
            if r2 not in ca_by_resid:
                continue
            d = float(np.linalg.norm(ca_by_resid[r1] - ca_by_resid[r2]))
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


def write_launcher_script():
    script = SET_DIR / 'run_fl_optimized.sh'
    with open(script, 'w') as f:
        f.write("""#!/bin/bash
# Launch all 25 fl_optimized replicates in parallel (or use GNU parallel)
set -e

CWD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_EXE="/home/sbali/miniconda3/envs/CALVADOS/bin/python"
PARALLEL=${1:-4}  # number of parallel jobs (default 4)

cd "$CWD"

# Collect replicate dirs
DIRS=$(ls -d seed-*_sample-* | sort)

# Run with xargs -P for parallelism
echo "$DIRS" | xargs -I {} -P $PARALLEL bash -c '
    cd "{}" || exit 1
    if [ -f "*.dcd" ] 2>/dev/null; then
        echo "{}: already has trajectory, skipping"
        exit 0
    fi
    echo "  Starting {}..."
    '"$PYTHON_EXE"' run.py > run.log 2>&1 || echo "  {}: FAILED"
    echo "  {}: done"
'

echo "All 25 replicates finished"
""")
    script.chmod(0o755)
    return script


def main():
    print("=" * 70)
    print(f"Preparing fl_optimized: 25 replicates with custom per-domain trims")
    print("=" * 70)

    # Make directories
    SET_DIR.mkdir(exist_ok=True)
    SHARED_INPUT.mkdir(exist_ok=True)

    # Copy/symlink shared files
    for fname, src in [('parp14.pdb', FL_PDB),
                       ('residues_CALVADOS3.csv', FL_RESIDUES)]:
        dest = SHARED_INPUT / fname
        if not dest.exists() and src.exists():
            shutil.copy2(src, dest)

    # Compute boundaries
    boundaries = compute_domain_boundaries()
    print("\nDomain boundaries after trim:")
    print(f"  {'Domain':<8} {'Trim':>4} {'Range':>14} {'Size':>5}")
    print("  " + "-" * 38)
    domain_list = []
    for dname in DOMAIN_TRIMS:
        if dname in boundaries:
            ds, de = boundaries[dname]
            trim = DOMAIN_TRIMS[dname]
            print(f"  {dname:<8} {trim:>4} {f'{ds}-{de}':>14} {de-ds+1:>5}")
            domain_list.append([ds, de])

    # Write shared domains.yaml
    with open(SHARED_INPUT / 'domains.yaml', 'w') as f:
        yaml.dump({'parp14': domain_list}, f, default_flow_style=True)
    print(f"\n  Wrote {SHARED_INPUT / 'domains.yaml'} "
          f"({len(domain_list)} domains)")

    # Generate KH7a-KHb custom restraints at trim 0
    pairs = generate_kh7a_khb_pairs(FL_PDB)
    cres_path = SHARED_INPUT / 'custom_restraints.txt'
    with open(cres_path, 'w') as f:
        for r1, r2, d_nm in sorted(pairs):
            f.write(f'parp14 1 {r1} | parp14 1 {r2} | '
                    f'{d_nm:.3f} {K_CUSTOM:.1f}\n')
    print(f"  Wrote {cres_path} ({len(pairs)} KH7a-KHb pairs at k={K_CUSTOM})")

    # Save summary metadata
    summary = {
        'set': 'fl_optimized',
        'n_replicates': len(list(SEEDS)) * len(list(SAMPLES)),
        'n_steps_per_replicate': N_STEPS,
        'k_harmonic': K_HARMONIC,
        'k_custom': K_CUSTOM,
        'domain_trims': DOMAIN_TRIMS,
        'domain_boundaries': {k: list(v) for k, v in boundaries.items()},
        'n_restraint_domains': len(domain_list),
        'n_restrained_residues': sum(de-ds+1 for ds, de in boundaries.values()),
        'n_kh7a_khb_custom_pairs': len(pairs),
        'kh7a_extent_for_pairs': [738, 789],
        'khb_extent_for_pairs': [1389, 1461],
        'structured_extents': STRUCTURED_EXTENTS,
    }
    with open(SET_DIR / 'simulation_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"  Wrote {SET_DIR / 'simulation_summary.json'}")

    # metadata.json so analysis scripts can use `--sim-folder fl_optimized`
    # (full-length, all 11 domain units). See sim_registry.py / stamp_metadata.py.
    with open(SET_DIR / 'metadata.json', 'w') as f:
        json.dump({'units': ['rrm1', 'rrm2', 'rrm3', 'kh1-kh6', 'kh7a', 'md1l1',
                             'md2', 'md3', 'khb-kh8', 'wwe', 'art'],
                   'sysname': 'parp14',
                   'sites': ['MD1', 'MD2', 'MD3', 'ART'],
                   'label': 'FL (optimized restraints)'}, f, indent=2)
    print(f"  Wrote {SET_DIR / 'metadata.json'}")

    # Load base config from CALVADOS defaults
    with open(DEFAULT_CONFIG) as f:
        base_config = yaml.safe_load(f)

    # Create 25 replicate directories
    print(f"\nCreating {len(list(SEEDS)) * len(list(SAMPLES))} replicates...")
    n_prepared = 0
    for seed in SEEDS:
        for sample in SAMPLES:
            rep_name = f'seed-{seed}_sample-{sample}'
            sim_dir = SET_DIR / rep_name
            input_dir = sim_dir / 'input'
            input_dir.mkdir(parents=True, exist_ok=True)

            # Symlink shared inputs
            for fname in ['parp14.pdb', 'residues_CALVADOS3.csv',
                          'domains.yaml', 'custom_restraints.txt']:
                dest = input_dir / fname
                src = SHARED_INPUT / fname
                if not dest.exists() and src.exists():
                    dest.symlink_to(src.resolve())

            # config.yaml — unique random seed per replicate
            random_seed = seed * 1000 + sample
            config = dict(base_config)
            config.update({
                'sysname': 'parp14',
                'box': [300, 300, 300],
                'temp': TEMP, 'ionic': IONIC, 'pH': PH,
                'topol': 'center',
                'steps': N_STEPS, 'wfreq': WFREQ, 'runtime': 0,
                'platform': PLATFORM, 'threads': THREADS,
                'restart': 'checkpoint', 'frestart': 'restart.chk',
                'verbose': True,
                'random_number_seed': random_seed,
                'custom_restraints': True,
                'custom_restraint_type': 'harmonic',
                'fcustom_restraints': 'input/custom_restraints.txt',
            })
            with open(sim_dir / 'config.yaml', 'w') as f:
                yaml.dump(config, f, default_flow_style=False)

            # components.yaml
            components_dict = {
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
                    'colabfold': 1,
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
                yaml.dump(components_dict, f, default_flow_style=False)

            write_run_py(sim_dir)
            n_prepared += 1

    print(f"  Prepared {n_prepared} replicates in {SET_DIR}/")

    # Write launcher
    launcher = write_launcher_script()
    print(f"\n  Wrote launcher: {launcher}")
    print(f"\nTo launch: bash {launcher} [N_PARALLEL]")
    print(f"  Default: 4 parallel jobs")
    print(f"\nProtocol: 25 x 20 ns = 500 ns total")
    print(f"  k_harmonic = {K_HARMONIC} kJ/mol/nm² (intra-domain)")
    print(f"  k_custom   = {K_CUSTOM} kJ/mol/nm² (KH7a-KHb {len(pairs)} pairs)")


if __name__ == '__main__':
    main()
