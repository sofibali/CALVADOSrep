#!/usr/bin/env python
"""
Prepare the COMPLETE MD1-MD3 CALVADOS run (FL 790-1388, 599 residues) -- the
contiguous macrodomain construct that INCLUDES the FL 1194-1206 MD2-MD3 linker
that the old `md` set dropped.

Why this exists: the old `md` set (and every fragment, incl. md1l1_md2_md3) was
built from DOMAIN_UNITS where md2=(1004,1193), md3=(1207,1388) -- so FL 1194-1206
belongs to no unit and is deleted (586 res, not 599). This script uses a single
contiguous FL block 790-1388 instead, so nothing is missing.

Mirrors the `md` set parameters (prepare_and_run_all.py): 25 replicates
(5 seeds x 5 samples), N_STEPS, harmonic restraints (k=700) on the folded cores,
colabfold=2 (AF3), box 80 nm. Extend to 20 ns via checkpoint restart like the md set.

Numbering: contiguous, local = FL - 789 for all residues. Restraint cores
(FL_RESTRAINT_DOMAINS mapped): MD1 791-978 -> 2-189, MD2 1003-1190 -> 214-401,
MD3 1216-1387 -> 427-598.

PREREQUISITE: run AF3 first (run_af3_md1md3_full.sh) so the structures exist at
  parp14/alphafold_outputs/md1l1_md2_md3_full/seed-{1-5}_sample-{0-4}/model.cif
Replicates without a model.cif are skipped (re-run this after AF3 finishes).

    python prepare_md_full.py          # prepare all available replicates
    bash md_full/run_all.sh parallel   # launch (won't auto-submit; see below)
"""
import os
import shutil
import yaml

# Reuse the md-set machinery + constants verbatim
from prepare_and_run_all import (
    Config, Components, convert_cif_to_pdb, create_run_script,
    FL_RESTRAINT_DOMAINS, SEEDS, SAMPLES, CWD,
    TEMP, IONIC, PH, THREADS, PLATFORM, WFREQ, N_STEPS, K_HARMONIC,
)

# --- the complete contiguous construct ---
FL_START, FL_END = 790, 1388          # contiguous, includes 1194-1206 linker
N_RES = FL_END - FL_START + 1         # 599
SYSNAME = 'parp14_md1md3_full'
BOX = 80                              # nm (mirror md set)
KEY = 'md_full'
AF3_DIR = '/home/sbali/CALVADOS/parp14/alphafold_outputs/md1l1_md2_md3_full'

def fl_to_local(fl):
    return fl - (FL_START - 1)        # single offset; no gaps

def macrodomain_cores():
    """MD1/MD2/MD3 folded restraint cores mapped to contiguous local numbering."""
    domains, labels = [], []
    for name, fl_s, fl_e, unit in FL_RESTRAINT_DOMAINS:
        if unit not in ('md1l1', 'md2', 'md3'):
            continue
        domains.append([fl_to_local(fl_s), fl_to_local(fl_e)])
        labels.append(name)
    return domains, labels

def main():
    set_dir = os.path.join(CWD, KEY)
    shared_input = os.path.join(set_dir, 'input')
    os.makedirs(shared_input, exist_ok=True)

    domains, labels = macrodomain_cores()
    print(f"Construct: FL {FL_START}-{FL_END} contiguous ({N_RES} res), incl. 1194-1206 linker")
    print("Restraint cores (FL -> local):")
    for (name, fl_s, fl_e, _), d in zip(
            [x for x in FL_RESTRAINT_DOMAINS if x[3] in ('md1l1','md2','md3')], domains):
        print(f"  {name}: FL {fl_s}-{fl_e} -> local {d[0]}-{d[1]}")

    # domains.yaml
    domains_yaml = os.path.join(shared_input, 'domains.yaml')
    with open(domains_yaml, 'w') as f:
        yaml.dump({SYSNAME: domains}, f, default_flow_style=False, sort_keys=False)

    # residues file (reuse md set's)
    residues_dest = os.path.join(shared_input, 'residues_CALVADOS3.csv')
    if not os.path.isfile(residues_dest):
        shutil.copy2(os.path.join(CWD, 'md', 'input', 'residues_CALVADOS3.csv'), residues_dest)

    n_prepared = n_skipped = 0
    for seed in SEEDS:
        for sample in SAMPLES:
            rep = f'seed-{seed}_sample-{sample}'
            sim_dir = os.path.join(set_dir, rep)
            input_dir = os.path.join(sim_dir, 'input')
            cif = os.path.join(AF3_DIR, rep, 'model.cif')
            conf = os.path.join(AF3_DIR, rep, 'confidences.json')
            if not os.path.isfile(cif):
                n_skipped += 1
                continue
            os.makedirs(input_dir, exist_ok=True)
            pdb_path = os.path.join(input_dir, f'{SYSNAME}.pdb')
            convert_cif_to_pdb(cif, pdb_path)
            if os.path.isfile(conf):
                shutil.copy2(conf, os.path.join(input_dir, f'{SYSNAME}.json'))

            config = Config(
                sysname=SYSNAME, box=[BOX, BOX, BOX],
                temp=TEMP, ionic=IONIC, pH=PH, topol='center',
                wfreq=WFREQ, steps=N_STEPS, runtime=0,
                platform=PLATFORM, threads=THREADS,
                restart='checkpoint', frestart='restart.chk', verbose=True,
                random_number_seed=seed * 1000 + sample,
            )
            config.write(sim_dir, name='config.yaml')

            components = Components(
                molecule_type='protein', nmol=1,
                restraint=True, charge_termini='both',
                fresidues=residues_dest, fdomains=domains_yaml,
                pdb_folder=input_dir,
                restraint_type='harmonic', use_com=True,
                colabfold=2, k_harmonic=K_HARMONIC,
            )
            components.add(name=SYSNAME)
            components.write(sim_dir, name='components.yaml')
            create_run_script(sim_dir)
            print(f"  prepared {KEY}/{rep}")
            n_prepared += 1

    print(f"\nPrepared {n_prepared} replicates, skipped {n_skipped} (no AF3 model.cif yet).")
    if n_skipped:
        print("Run AF3 (run_af3_md1md3_full.sh), then re-run this script for the rest.")

    # simple launcher
    run_sh = os.path.join(set_dir, 'run_all.sh')
    with open(run_sh, 'w') as f:
        f.write('#!/bin/bash\n# Launch all prepared md_full replicates. Arg: "parallel" or "serial" (default serial).\n')
        f.write('set -e\nHERE="$(cd "$(dirname "$0")" && pwd)"\nMODE="${1:-serial}"\n')
        f.write('for d in "$HERE"/seed-*_sample-*/; do\n')
        f.write('  [ -f "$d/run.py" ] || continue\n')
        f.write('  if [ "$MODE" = parallel ]; then ( cd "$d" && python run.py & ); else ( cd "$d" && python run.py ); fi\n')
        f.write('done\nwait\n')
    os.chmod(run_sh, 0o755)
    print(f"Launcher: {run_sh}  (bash run_all.sh parallel)")


if __name__ == "__main__":
    main()
