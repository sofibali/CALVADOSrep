#!/usr/bin/env python
"""
Prepare the COMPLETE MD1-ART CALVADOS run (FL 790-1801, 1012 residues) -- the
contiguous MD1-through-ART construct that INCLUDES the FL 1194-1206 MD2-MD3
linker that the old `mka` set dropped (same gap prepare_md_full.py fixed for
MD1-MD3; see that script's docstring for why the gap exists at all).

Structure source: sliced directly from input/parp14.pdb (the same full-length
AF2 model used by the `fl`/`fl_optimized` sets), residues 790-1801 renumbered
1-1012. NOT a fresh AF3 prediction -- this machine has no GPU, and AF3
inference needs one. Reusing one shared slice across all replicates mirrors
the `fl` set's own convention for full-length-scale constructs (colabfold=1,
identical starting PDB every replicate, only the Langevin RNG seed differs;
see prepare_and_run_all.py's prepare_fl()). The FL AF2 model already has real
predicted coordinates for 1194-1206 (it was never actually missing/disordered,
just excluded from the AF3 combinatorial sub-construct library).

5 replicates only (seed 1-5, sample-0), each run FRESH to steps=100,000,000
(1 us at dt=0.01 ps) in one continuous run -- not a checkpoint extension of
anything pre-existing, since no mka_full replicates existed before this.
threads=16 (up from the project's usual 4) since this machine has 256 cores
and only ~10 jobs (this + the 5 md_full extensions) run concurrently.

Restraint cores: MD1/MD2/MD3/WWE/ART structured cores (FL_RESTRAINT_DOMAINS
filtered to the units this construct actually has -- md1l1/md2/md3/wwe/art;
no khb-kh8 core, matching how the existing `mka` set also left khb-kh8
unrestrained), mapped through a SINGLE contiguous offset (local = FL - 789,
no gaps to skip).

    python prepare_mka_full.py           # slice PDB + prepare all 5 replicates
    bash mka_full/run_all.sh parallel    # launch
"""
import os
import yaml

from prepare_and_run_all import (
    Config, Components, create_run_script,
    FL_RESTRAINT_DOMAINS, CWD,
    TEMP, IONIC, PH, PLATFORM, WFREQ,
)

FL_START, FL_END = 790, 1801           # contiguous, includes 1194-1206 linker
N_RES = FL_END - FL_START + 1          # 1012
SYSNAME = 'parp14_mka_full'
BOX = 100                              # nm (mirror mka set's box)
KEY = 'mka_full'
SEEDS = range(1, 6)                    # 5 replicates only (not the usual 25)
N_STEPS = 100_000_000                  # 1 us at dt=0.01 ps
THREADS = 16
UNITS_PRESENT = ('md1l1', 'md2', 'md3', 'wwe', 'art')  # no khb-kh8 core, matches `mka`


def fl_to_local(fl):
    return fl - (FL_START - 1)


def slice_source_pdb(src_pdb, dst_pdb):
    n_atoms = 0
    resids = set()
    with open(src_pdb) as f, open(dst_pdb, 'w') as out:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                resi = int(line[22:26])
                if FL_START <= resi <= FL_END:
                    new_resi = resi - (FL_START - 1)
                    out.write(f"{line[:22]}{new_resi:4d}{line[26:]}")
                    n_atoms += 1
                    resids.add(new_resi)
            elif line.startswith(('TER', 'END')):
                out.write(line)
    assert resids == set(range(1, N_RES + 1)), \
        f"sliced PDB residue range mismatch: got {min(resids)}-{max(resids)}, expected 1-{N_RES}"
    return n_atoms


def restraint_cores():
    domains, labels = [], []
    for name, fl_s, fl_e, unit in FL_RESTRAINT_DOMAINS:
        if unit not in UNITS_PRESENT:
            continue
        domains.append([fl_to_local(fl_s), fl_to_local(fl_e)])
        labels.append(name)
    return domains, labels


def main():
    set_dir = os.path.join(CWD, KEY)
    shared_input = os.path.join(set_dir, 'input')
    os.makedirs(shared_input, exist_ok=True)

    domains, labels = restraint_cores()
    print(f"Construct: FL {FL_START}-{FL_END} contiguous ({N_RES} res), incl. 1194-1206 linker")
    print("Restraint cores (FL -> local):")
    for (name, fl_s, fl_e, _), d in zip(
            [x for x in FL_RESTRAINT_DOMAINS if x[3] in UNITS_PRESENT], domains):
        print(f"  {name}: FL {fl_s}-{fl_e} -> local {d[0]}-{d[1]}")

    domains_yaml = os.path.join(shared_input, 'domains.yaml')
    with open(domains_yaml, 'w') as f:
        yaml.dump({SYSNAME: domains}, f, default_flow_style=False, sort_keys=False)

    residues_dest = os.path.join(shared_input, 'residues_CALVADOS3.csv')
    if not os.path.isfile(residues_dest):
        import shutil
        shutil.copy2(os.path.join(CWD, 'input', 'residues_CALVADOS3.csv'), residues_dest)

    shared_pdb = os.path.join(shared_input, f'{SYSNAME}.pdb')
    n_atoms = slice_source_pdb(os.path.join(CWD, 'input', 'parp14.pdb'), shared_pdb)
    print(f"Sliced shared structure: {n_atoms} atoms -> {shared_pdb}")

    n_prepared = 0
    for seed in SEEDS:
        rep = f'seed-{seed}_sample-0'
        sim_dir = os.path.join(set_dir, rep)
        input_dir = os.path.join(sim_dir, 'input')
        os.makedirs(input_dir, exist_ok=True)

        pdb_path = os.path.join(input_dir, f'{SYSNAME}.pdb')
        if not os.path.isfile(pdb_path):
            import shutil
            shutil.copy2(shared_pdb, pdb_path)

        config = Config(
            sysname=SYSNAME, box=[BOX, BOX, BOX],
            temp=TEMP, ionic=IONIC, pH=PH, topol='center',
            wfreq=WFREQ, steps=N_STEPS, runtime=0,
            platform=PLATFORM, threads=THREADS,
            restart='checkpoint', frestart='restart.chk', verbose=True,
            random_number_seed=seed * 1000,
        )
        config.write(sim_dir, name='config.yaml')

        components = Components(
            molecule_type='protein', nmol=1,
            restraint=True, charge_termini='both',
            fresidues=residues_dest, fdomains=domains_yaml,
            pdb_folder=input_dir,
            restraint_type='harmonic', use_com=True,
            colabfold=1,  # shared FL-model slice, no per-replicate AF3 confidence file
        )
        components.add(name=SYSNAME)
        components.write(sim_dir, name='components.yaml')
        create_run_script(sim_dir)
        print(f"  prepared {KEY}/{rep}")
        n_prepared += 1

    print(f"\nPrepared {n_prepared} replicates ({N_STEPS * 0.01 / 1000:.0f} ns / {N_STEPS} steps each).")

    run_sh = os.path.join(set_dir, 'run_all.sh')
    with open(run_sh, 'w') as f:
        f.write('#!/bin/bash\n# Launch all prepared mka_full replicates. Arg: "parallel" or "serial" (default serial).\n')
        f.write('set -e\nHERE="$(cd "$(dirname "$0")" && pwd)"\nMODE="${1:-serial}"\n')
        f.write('PY=/home/sbali/miniconda3/envs/CALVADOS/bin/python\n')
        f.write('for d in "$HERE"/seed-*_sample-*/; do\n')
        f.write('  [ -f "$d/run.py" ] || continue\n')
        f.write('  if [ "$MODE" = parallel ]; then ( cd "$d" && $PY run.py ) & else ( cd "$d" && $PY run.py ); fi\n')
        f.write('done\n[ "$MODE" = parallel ] && wait\necho "mka_full done."\n')
    os.chmod(run_sh, 0o755)
    print(f"Launcher: {run_sh}  (bash run_all.sh parallel)")


if __name__ == "__main__":
    main()
