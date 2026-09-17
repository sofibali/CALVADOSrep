#!/usr/bin/env python
"""
Prepare 9 additional contiguous CALVADOS constructs (no residues excised
within their span -- same convention as prepare_md_full.py/
prepare_mka_full.py/prepare_core_full_go.py), sliced from the shared FL AF2
model (input/parp14.pdb), 5 replicates x 1 us each:

  kh1_art_full    FL 315-1801   (KH1-6 through ART)
  kh1_wwe_full    FL 315-1602   (KH1-6 through WWE, no ART)
  md2_art_full    FL 1004-1801  (MD2 through ART)
  md2_wwe_full    FL 1004-1602  (MD2 through WWE, no ART)
  md3_art_full    FL 1207-1801  (MD3 through ART)
  md3_wwe_full    FL 1207-1602  (MD3 through WWE, no ART)
  core_wwe_full_go FL 738-1602  (KH7a through WWE, no ART -- truncated core_full_go)
  mka_wwe_full    FL 790-1602   (MD1L1 through WWE, no ART -- truncated mka_full)
  fl_wwe_full_go  FL 1-1602     (RRM1 through WWE, no ART -- truncated FL)

Go-model KH7a-KHb custom restraint (same 141 pairs/distances as fl_go, k=15,
remapped to local numbering) is added to any construct spanning BOTH KH7a
(738-789) and KHb-KH8 (1389-1533) together -- kh1_art_full, kh1_wwe_full,
core_wwe_full_go, fl_wwe_full_go. The others only contain one KH piece (or
neither), so plain harmonic domain restraints suffice (same reasoning as
prepare_core_full_go.py's docstring).

Restraint cores are FL_RESTRAINT_DOMAINS (RRM1/RRM2/RRM3/MD1/MD2/MD3/WWE/ART)
CLIPPED to each construct's own span -- necessary because a couple of the
restraint-core boundaries fall slightly outside their nominal DOMAIN_UNITS
range (e.g. the MD2 restraint core starts at FL 1003, one residue before the
'md2' unit's own FL 1004 start), which would otherwise map to an invalid
(zero or negative) local residue index for constructs that start exactly at
a unit boundary like md2_art_full.

    python prepare_extra_full_constructs.py             # slice + prepare all 9
    bash <construct>/run_all.sh parallel                 # launch one
"""
import json
import os
import shutil
import yaml

from prepare_and_run_all import (
    Config, Components, create_run_script,
    FL_RESTRAINT_DOMAINS, CWD,
    TEMP, IONIC, PH, PLATFORM, WFREQ,
)

SEEDS = range(1, 6)
N_STEPS = 100_000_000  # 1 us at dt=0.01 ps
THREADS = 16

GO_SOURCE = os.path.join(CWD, 'fl_go', 'state-1_tica_seed-1_sample-1_fr777',
                          'input', 'custom_restraints_go.txt')

CONSTRUCTS = [
    # key, fl_start, fl_end, box_nm, units_present, needs_go
    {'key': 'kh1_art_full', 'fl_start': 315, 'fl_end': 1801, 'box': 250,
     'units': ('kh1-kh6', 'kh7a', 'md1l1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art'),
     'needs_go': True},
    {'key': 'kh1_wwe_full', 'fl_start': 315, 'fl_end': 1602, 'box': 200,
     'units': ('kh1-kh6', 'kh7a', 'md1l1', 'md2', 'md3', 'khb-kh8', 'wwe'),
     'needs_go': True},
    {'key': 'md2_art_full', 'fl_start': 1004, 'fl_end': 1801, 'box': 90,
     'units': ('md2', 'md3', 'khb-kh8', 'wwe', 'art'), 'needs_go': False},
    {'key': 'md2_wwe_full', 'fl_start': 1004, 'fl_end': 1602, 'box': 80,
     'units': ('md2', 'md3', 'khb-kh8', 'wwe'), 'needs_go': False},
    {'key': 'md3_art_full', 'fl_start': 1207, 'fl_end': 1801, 'box': 80,
     'units': ('md3', 'khb-kh8', 'wwe', 'art'), 'needs_go': False},
    {'key': 'md3_wwe_full', 'fl_start': 1207, 'fl_end': 1602, 'box': 70,
     'units': ('md3', 'khb-kh8', 'wwe'), 'needs_go': False},
    {'key': 'core_wwe_full_go', 'fl_start': 738, 'fl_end': 1602, 'box': 110,
     'units': ('kh7a', 'md1l1', 'md2', 'md3', 'khb-kh8', 'wwe'), 'needs_go': True},
    {'key': 'mka_wwe_full', 'fl_start': 790, 'fl_end': 1602, 'box': 90,
     'units': ('md1l1', 'md2', 'md3', 'khb-kh8', 'wwe'), 'needs_go': False},
    {'key': 'fl_wwe_full_go', 'fl_start': 1, 'fl_end': 1602, 'box': 280,
     'units': ('rrm1', 'rrm2', 'rrm3', 'kh1-kh6', 'kh7a', 'md1l1', 'md2', 'md3',
               'khb-kh8', 'wwe'), 'needs_go': True},
]


def slice_source_pdb(src_pdb, dst_pdb, fl_start, fl_end):
    n_atoms = 0
    resids = set()
    with open(src_pdb) as f, open(dst_pdb, 'w') as out:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                resi = int(line[22:26])
                if fl_start <= resi <= fl_end:
                    new_resi = resi - (fl_start - 1)
                    out.write(f"{line[:22]}{new_resi:4d}{line[26:]}")
                    n_atoms += 1
                    resids.add(new_resi)
            elif line.startswith(('TER', 'END')):
                out.write(line)
    n_res = fl_end - fl_start + 1
    assert resids == set(range(1, n_res + 1)), \
        f"sliced PDB residue range mismatch: got {min(resids)}-{max(resids)}, expected 1-{n_res}"
    return n_atoms


def restraint_cores(fl_start, fl_end, units_present):
    """FL_RESTRAINT_DOMAINS boxes whose primary_unit is in units_present,
    CLIPPED to [fl_start, fl_end] (see module docstring for why clipping is
    needed), mapped to local (1-based) numbering."""
    domains, labels = [], []
    for name, fl_s, fl_e, unit in FL_RESTRAINT_DOMAINS:
        if unit not in units_present:
            continue
        clipped_s = max(fl_s, fl_start)
        clipped_e = min(fl_e, fl_end)
        if clipped_s > clipped_e:
            continue
        local_s = clipped_s - (fl_start - 1)
        local_e = clipped_e - (fl_start - 1)
        domains.append([local_s, local_e])
        labels.append(name)
    return domains, labels


def remap_custom_restraints(src_path, dst_path, fl_start, fl_end, sysname):
    n = 0
    with open(src_path) as f, open(dst_path, 'w') as out:
        for line in f:
            left, right, params = line.strip().split('|')
            _, chain, r1 = left.split()
            _, _, r2 = right.split()
            r1, r2 = int(r1), int(r2)
            if not (fl_start <= r1 <= fl_end and fl_start <= r2 <= fl_end):
                continue  # pair falls outside this (shorter) construct's span
            r1_local = r1 - (fl_start - 1)
            r2_local = r2 - (fl_start - 1)
            out.write(f'{sysname} {chain} {r1_local} | {sysname} {chain} {r2_local} | {params.strip()}\n')
            n += 1
    return n


def prepare_one(spec):
    key = spec['key']
    fl_start, fl_end = spec['fl_start'], spec['fl_end']
    n_res = fl_end - fl_start + 1
    sysname = f'parp14_{key}'
    box = spec['box']
    units_present = spec['units']

    set_dir = os.path.join(CWD, key)
    shared_input = os.path.join(set_dir, 'input')
    os.makedirs(shared_input, exist_ok=True)

    domains, labels = restraint_cores(fl_start, fl_end, units_present)
    print(f"\n=== {key}: FL {fl_start}-{fl_end} ({n_res} res) ===")
    print(f"  Restraint cores: {list(zip(labels, domains))}")

    domains_yaml = os.path.join(shared_input, 'domains.yaml')
    with open(domains_yaml, 'w') as f:
        yaml.dump({sysname: domains}, f, default_flow_style=False, sort_keys=False)

    residues_dest = os.path.join(shared_input, 'residues_CALVADOS3.csv')
    if not os.path.isfile(residues_dest):
        shutil.copy2(os.path.join(CWD, 'input', 'residues_CALVADOS3.csv'), residues_dest)

    shared_pdb = os.path.join(shared_input, f'{sysname}.pdb')
    n_atoms = slice_source_pdb(os.path.join(CWD, 'input', 'parp14.pdb'), shared_pdb,
                                fl_start, fl_end)
    print(f"  Sliced shared structure: {n_atoms} atoms -> {shared_pdb}")

    cres_path = None
    n_pairs = 0
    if spec['needs_go']:
        cres_path = os.path.join(shared_input, 'custom_restraints_go.txt')
        n_pairs = remap_custom_restraints(GO_SOURCE, cres_path, fl_start, fl_end, sysname)
        print(f"  Remapped {n_pairs} KH7a-KHb Go-model pairs -> {cres_path}")

    n_prepared = 0
    for seed in SEEDS:
        rep = f'seed-{seed}_sample-0'
        sim_dir = os.path.join(set_dir, rep)
        input_dir = os.path.join(sim_dir, 'input')
        os.makedirs(input_dir, exist_ok=True)

        pdb_path = os.path.join(input_dir, f'{sysname}.pdb')
        if not os.path.isfile(pdb_path):
            shutil.copy2(shared_pdb, pdb_path)
        if cres_path is not None:
            cres_dest = os.path.join(input_dir, 'custom_restraints_go.txt')
            if not os.path.isfile(cres_dest):
                shutil.copy2(cres_path, cres_dest)

        config_kwargs = dict(
            sysname=sysname, box=[box, box, box],
            temp=TEMP, ionic=IONIC, pH=PH, topol='center',
            wfreq=WFREQ, steps=N_STEPS, runtime=0,
            platform=PLATFORM, threads=THREADS,
            restart='checkpoint', frestart='restart.chk', verbose=True,
            random_number_seed=seed * 1000,
        )
        if spec['needs_go']:
            config_kwargs.update(custom_restraints=True, custom_restraint_type='harmonic',
                                  fcustom_restraints='input/custom_restraints_go.txt')
        config = Config(**config_kwargs)
        config.write(sim_dir, name='config.yaml')

        components = Components(
            molecule_type='protein', nmol=1,
            restraint=True, charge_termini='both',
            fresidues=residues_dest, fdomains=domains_yaml,
            pdb_folder=input_dir,
            restraint_type='harmonic', use_com=True,
            colabfold=1,
        )
        components.add(name=sysname)
        components.write(sim_dir, name='components.yaml')
        create_run_script(sim_dir)
        n_prepared += 1

    print(f"  Prepared {n_prepared} replicates ({N_STEPS * 0.01 / 1000:.0f} ns / {N_STEPS} steps each).")

    with open(os.path.join(set_dir, 'metadata.json'), 'w') as f:
        json.dump({'units': list(units_present), 'sysname': sysname,
                   'label': f'{key} (FL {fl_start}-{fl_end}, {n_res} res)'}, f, indent=2)

    run_sh = os.path.join(set_dir, 'run_all.sh')
    with open(run_sh, 'w') as f:
        f.write(f'#!/bin/bash\n# Launch all prepared {key} replicates. Arg: "parallel" or "serial" (default serial).\n')
        f.write('set -e\nHERE="$(cd "$(dirname "$0")" && pwd)"\nMODE="${1:-serial}"\n')
        f.write('PY=/home/sbali/miniconda3/envs/CALVADOS/bin/python\n')
        f.write('for d in "$HERE"/seed-*_sample-*/; do\n')
        f.write('  [ -f "$d/run.py" ] || continue\n')
        f.write('  if [ "$MODE" = parallel ]; then ( cd "$d" && $PY run.py ) & else ( cd "$d" && $PY run.py ); fi\n')
        f.write(f'done\n[ "$MODE" = parallel ] && wait\necho "{key} done."\n')
    os.chmod(run_sh, 0o755)
    print(f"  Launcher: {run_sh}")


def main():
    for spec in CONSTRUCTS:
        prepare_one(spec)


if __name__ == "__main__":
    main()
