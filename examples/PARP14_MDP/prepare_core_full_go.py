#!/usr/bin/env python
"""
Prepare the contiguous KH7a-ART CALVADOS run (FL 738-1801, 1064 residues) --
the same domain span as the existing `core` set (KH7a+MD1L1+MD2+MD3+KHb-KH8+
WWE+ART) but CONTIGUOUS (no inter-unit linkers excised), mirroring how
prepare_md_full.py/prepare_mka_full.py made contiguous versions of `md`/`mka`.

PLUS Go-model custom restraints between KH7a and KHb-KH8: these two domains
are adjacent in the folded structure but split apart in sequence (KH7a at
738-789, KHb-KH8 at 1389-1533), so a plain per-domain harmonic restraint
(which only restrains residues WITHIN a domain) can't hold them together --
that's the same reason fl_optimized/fl_go both carry the identical 141-pair
KH7a-KHb custom restraint file. Since core_full_go's 738-1801 span contains
BOTH domains (unlike md_full/mka_full, which start at 790 and never include
KH7a), it needs the same custom restraint. `md`/`core`/`mka` (as originally
built) don't include it only because nobody added it there, not because
they don't need it structurally -- core_full_go corrects that.

The pairs are reused as-is from fl_go's custom_restraints_go.txt (same 141
CA-CA pairs, same target distances, same Go-model k=15 kJ/mol/nm^2 -- soft
compared to fl_optimized's k=350 harmonic strength), just remapped from FL
numbering to this construct's own local numbering (local = FL - 737). All
141 pairs fall within 738-1454, safely inside the 738-1801 span, so the
remap is a simple integer offset with no pairs dropped.

Structure source: sliced directly from input/parp14.pdb (the same full-length
AF2 model used by fl/fl_optimized/md_full/mka_full), residues 738-1801
renumbered 1-1064. Reusing one shared slice across all replicates mirrors the
fl/md_full/mka_full convention (colabfold=1, identical starting PDB every
replicate, only the Langevin RNG seed differs).

5 replicates (seed 1-5, sample-0), each run FRESH to steps=100,000,000 (1 us
at dt=0.01 ps) -- matches mka_full's replicate scheme.

    python prepare_core_full_go.py            # slice PDB + prepare all 5 replicates
    bash core_full_go/run_all.sh parallel      # launch
"""
import json
import os
import yaml

from prepare_and_run_all import (
    Config, Components, create_run_script,
    FL_RESTRAINT_DOMAINS, CWD,
    TEMP, IONIC, PH, PLATFORM, WFREQ,
)

FL_START, FL_END = 738, 1801            # contiguous KH7a-ART span
N_RES = FL_END - FL_START + 1           # 1064
SYSNAME = 'parp14_core_full_go'
BOX = 120                               # nm (mirror `core` set's box)
KEY = 'core_full_go'
SEEDS = range(1, 6)                     # 5 replicates, matches mka_full
N_STEPS = 100_000_000                   # 1 us at dt=0.01 ps
THREADS = 16
UNITS_PRESENT = ('kh7a', 'md1l1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art')
# FL_RESTRAINT_DOMAINS has no kh7a/khb-kh8 entries project-wide (no structured
# restraint core defined for either KH domain anywhere), so only
# md1l1/md2/md3/wwe/art actually produce restraint boxes below -- matches
# how `core`/`mka`/`mka_full` also leave both KH domains unrestrained.

GO_SOURCE = os.path.join(CWD, 'fl_go', 'state-1_tica_seed-1_sample-1_fr777',
                          'input', 'custom_restraints_go.txt')


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


def remap_custom_restraints(src_path, dst_path):
    n = 0
    with open(src_path) as f, open(dst_path, 'w') as out:
        for line in f:
            left, right, params = line.strip().split('|')
            _, chain, r1 = left.split()
            _, _, r2 = right.split()
            r1_local = fl_to_local(int(r1))
            r2_local = fl_to_local(int(r2))
            assert 1 <= r1_local <= N_RES and 1 <= r2_local <= N_RES, \
                f"custom restraint pair ({r1},{r2}) falls outside {FL_START}-{FL_END}"
            out.write(f'{SYSNAME} {chain} {r1_local} | {SYSNAME} {chain} {r2_local} | {params.strip()}\n')
            n += 1
    return n


def main():
    set_dir = os.path.join(CWD, KEY)
    shared_input = os.path.join(set_dir, 'input')
    os.makedirs(shared_input, exist_ok=True)

    domains, labels = restraint_cores()
    print(f"Construct: FL {FL_START}-{FL_END} contiguous ({N_RES} res), incl. KH7a + KHb-KH8")
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

    cres_path = os.path.join(shared_input, 'custom_restraints_go.txt')
    n_pairs = remap_custom_restraints(GO_SOURCE, cres_path)
    print(f"Remapped {n_pairs} KH7a-KHb Go-model pairs (FL -> local) -> {cres_path}")

    n_prepared = 0
    for seed in SEEDS:
        rep = f'seed-{seed}_sample-0'
        sim_dir = os.path.join(set_dir, rep)
        input_dir = os.path.join(sim_dir, 'input')
        os.makedirs(input_dir, exist_ok=True)

        import shutil
        pdb_path = os.path.join(input_dir, f'{SYSNAME}.pdb')
        if not os.path.isfile(pdb_path):
            shutil.copy2(shared_pdb, pdb_path)
        cres_dest = os.path.join(input_dir, 'custom_restraints_go.txt')
        if not os.path.isfile(cres_dest):
            shutil.copy2(cres_path, cres_dest)

        config = Config(
            sysname=SYSNAME, box=[BOX, BOX, BOX],
            temp=TEMP, ionic=IONIC, pH=PH, topol='center',
            wfreq=WFREQ, steps=N_STEPS, runtime=0,
            platform=PLATFORM, threads=THREADS,
            restart='checkpoint', frestart='restart.chk', verbose=True,
            random_number_seed=seed * 1000,
            custom_restraints=True, custom_restraint_type='harmonic',
            fcustom_restraints='input/custom_restraints_go.txt',
        )
        config.write(sim_dir, name='config.yaml')

        components = Components(
            molecule_type='protein', nmol=1,
            restraint=True, charge_termini='both',
            fresidues=residues_dest, fdomains=domains_yaml,
            pdb_folder=input_dir,
            restraint_type='harmonic', use_com=True,
            colabfold=1,
        )
        components.add(name=SYSNAME)
        components.write(sim_dir, name='components.yaml')
        create_run_script(sim_dir)
        print(f"  prepared {KEY}/{rep}")
        n_prepared += 1

    print(f"\nPrepared {n_prepared} replicates ({N_STEPS * 0.01 / 1000:.0f} ns / {N_STEPS} steps each).")

    with open(os.path.join(set_dir, 'metadata.json'), 'w') as f:
        json.dump({'units': list(UNITS_PRESENT),
                   'sysname': SYSNAME,
                   'sites': ['MD1', 'MD2', 'MD3', 'ART'],
                   'label': 'KH7a-ART contiguous + Go-model KH7a-KHb restraints (738-1801, 1064 res)'},
                  f, indent=2)

    run_sh = os.path.join(set_dir, 'run_all.sh')
    with open(run_sh, 'w') as f:
        f.write('#!/bin/bash\n# Launch all prepared core_full_go replicates. Arg: "parallel" or "serial" (default serial).\n')
        f.write('set -e\nHERE="$(cd "$(dirname "$0")" && pwd)"\nMODE="${1:-serial}"\n')
        f.write('PY=/home/sbali/miniconda3/envs/CALVADOS/bin/python\n')
        f.write('for d in "$HERE"/seed-*_sample-*/; do\n')
        f.write('  [ -f "$d/run.py" ] || continue\n')
        f.write('  if [ "$MODE" = parallel ]; then ( cd "$d" && $PY run.py ) & else ( cd "$d" && $PY run.py ); fi\n')
        f.write('done\n[ "$MODE" = parallel ] && wait\necho "core_full_go done."\n')
    os.chmod(run_sh, 0o755)
    print(f"Launcher: {run_sh}  (bash run_all.sh parallel)")


if __name__ == "__main__":
    main()
