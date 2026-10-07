#!/usr/bin/env python3
"""Null controls: are the inter-chain contacts binding, or just collisions?

The campaign measures bound fractions of 0.06-0.18 with contact lifetimes of
~0.25 ns. That is consistent with weak binding AND with pure excluded-volume
collision statistics, and the two are not distinguishable from the production
runs alone. These controls separate them.

The lever is the per-residue lambda in residues_CALVADOS3.csv. In CALVADOS's
Ashbaugh-Hatch form (interactions.py:30) lambda scales the ATTRACTIVE tail:

    eps * select(step(r - 2^(1/6)s),
                 4*l*((s/r)^12 - (s/r)^6 - shift),          <- r > sigma: tail
                 4*((s/r)^12 - (s/r)^6 - l*shift) + (1-l))  <- r < sigma: core

Set every lambda to 0 and the tail vanishes: pure WCA repulsion, excluded volume
preserved. NOTE the config key `fixed_lambda` does NOT do this -- line 32 reads
`l = select(id1+id2, (id1*id2)*0.5*(l1+l2), fixed_lambda)`, so fixed_lambda only
applies when BOTH bead ids are 0, never for protein beads.

Two variants, which bracket the answer:

  null_lambda  lambda = 0, charges kept.  Removes hydrophobic attraction but
               leaves Debye-Huckel electrostatics, so oppositely charged patches
               can still attract. Isolates the hydrophobic contribution.
  null_steric  lambda = 0 AND q = 0.  Pure excluded volume (plus bonds and the
               domain restraints). Whatever contact remains here is the
               collision background and nothing else.

Read the result as: production - null_lambda = hydrophobic binding;
null_lambda - null_steric = electrostatic steering; null_steric = geometry.
If production ~= null_steric, there is no binding to discuss.

Usage:  python prepare_null_controls.py [--sets ...] [--nreps 5]
"""
import shutil
from pathlib import Path
from argparse import ArgumentParser
import pandas as pd

HERE = Path(__file__).resolve().parent
SRC = HERE / 'binding'
DST = HERE / 'binding_null'

# Span the measured range: strongest homotypic, a strong heterotypic, a weak one.
DEFAULT_SETS = ['p14_homo', 'p14_dtx3l', 'p9_dtx3l']
VARIANTS = {'null_lambda': dict(zero_lambda=True, zero_charge=False),
            'null_steric': dict(zero_lambda=True, zero_charge=True)}


def make_residues(src_csv, out_csv, zero_lambda, zero_charge):
    d = pd.read_csv(src_csv)
    if zero_lambda:
        d['lambdas'] = 0.0
    if zero_charge:
        d['q'] = 0.0
    d.to_csv(out_csv, index=False)
    return len(d)


def main():
    ap = ArgumentParser()
    ap.add_argument('--sets', nargs='*', default=DEFAULT_SETS)
    ap.add_argument('--nreps', type=int, default=5)
    a = ap.parse_args()

    made = 0
    for var, opt in VARIANTS.items():
        for s in a.sets:
            sd = SRC / s
            if not sd.is_dir():
                print(f'  skip {s}: not found'); continue
            for rep in sorted(sd.glob('rep-*'),
                              key=lambda p: int(p.name.split('-')[1]))[:a.nreps]:
                out = DST / f'{var}_{s}' / rep.name
                if (out / 'config.yaml').exists():
                    continue
                out.mkdir(parents=True, exist_ok=True)
                # copy dereferenced so the null tree is self-contained
                for item in rep.iterdir():
                    if item.name in ('config.yaml', 'components.yaml', 'run.py'):
                        shutil.copy2(item, out / item.name)
                    elif item.name == 'input':
                        shutil.copytree(item, out / 'input', symlinks=False,
                                        dirs_exist_ok=True)
                n = make_residues(rep / 'input' / 'residues_CALVADOS3.csv',
                                  out / 'input' / 'residues_CALVADOS3.csv', **opt)
                made += 1
            print(f'  {var}_{s}: {a.nreps} replicates ({n} residue types zeroed)')
    print(f'\nprepared {made} replicate directories under {DST}')


if __name__ == '__main__':
    main()
