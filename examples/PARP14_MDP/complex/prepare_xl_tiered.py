#!/usr/bin/env python3
"""Crosslink-restrained CALVADOS runs, rebuilt on the CORRECT data.

WHAT WAS WRONG WITH THE FIRST SET (binding_go, already run, 20 replicates)
--------------------------------------------------------------------------
Two errors, both now fixed:

1. It used 3 crosslinks. The supplement holds 51 inter-protein crosslinks; the
   first extraction read only the full-length sheet and hard-coded a column
   index, so the two truncated-DTX3L sheets (26 and 22 crosslinks) silently
   returned nothing. With the per-replicate 'x' marks now readable, 30 of the
   51 are seen in >=2 of the 4 BS3 replicates.

2. Every restraint got the same r = 2.0 nm. Following IMProv (Ziemianowicz et
   al., Mol Cell Proteomics 2021;20:100139), a crosslink formed in a flexible
   region can be kinetically trapped and reports a transient excursion rather
   than the equilibrium structure, so it warrants a looser threshold than one
   joining two stable folds. IMProv's offset scheme: 25 A stable-stable,
   40 A when a partner is unstable.

TIERING RULE
------------
An end is "rigid" if it falls inside a restrained block of
binding/p9_dtx3l/input/domains.yaml, "flexible" otherwise.
    both rigid  -> r = 2.5 nm   (IMProv stable-stable)
    otherwise   -> r = 4.0 nm   (IMProv unstable)

Restraints are HARMONIC so the 85 pre-existing intra-PARP9 split-KH1 restraints
stay byte-identical to the unrestrained control -- CALVADOS has one global
custom_restraint_type, and switching to Go would reinterpret those as 144 kT
wells ~300x stiffer than the harmonic they replace. The control for the
difference map is the existing unrestrained binding/p9_dtx3l.

Usage:  python prepare_xl_tiered.py [--nreps 5] [--ks 20]
"""
import shutil
from pathlib import Path
from argparse import ArgumentParser
import pandas as pd, yaml

HERE = Path(__file__).resolve().parent
SRC = HERE / 'binding' / 'p9_dtx3l'
DST = HERE / 'binding_xl'
XL = HERE / 'data' / 'published_xlinks_v2.csv'
R_RIGID, R_FLEX = 2.5, 4.0


def rigid(dom, chain, resid):
    return any(lo <= resid <= hi for lo, hi in dom[chain])


def main():
    ap = ArgumentParser()
    ap.add_argument('--nreps', type=int, default=5)
    ap.add_argument('--ks', type=float, nargs='*', default=[20.0])
    a = ap.parse_args()

    dom = yaml.safe_load(open(SRC / 'input' / 'domains.yaml'))
    d = pd.read_csv(XL)
    d = d[(d.block == 'P9-Dtx3L') & (d.repro_2of4.fillna(False).astype(bool))]
    d = d.drop_duplicates(subset=['res_a', 'res_b']).reset_index(drop=True)
    tiers = []
    for _, r in d.iterrows():
        both = rigid(dom, 'parp9', r.res_a) and rigid(dom, 'dtx3l', r.res_b)
        tiers.append((int(r.res_a), int(r.res_b), R_RIGID if both else R_FLEX,
                      'rigid-rigid' if both else 'flexible'))
    nr = sum(1 for t in tiers if t[3] == 'rigid-rigid')
    print(f'{len(tiers)} unique reproducible inter-protein crosslinks: '
          f'{nr} rigid-rigid (r={R_RIGID} nm), {len(tiers)-nr} flexible (r={R_FLEX} nm)')

    made = 0
    for k in a.ks:
        tag = f'p9_dtx3l_xl30_k{int(k)}'
        for rep in sorted(SRC.glob('rep-*'),
                          key=lambda p: int(p.name.split('-')[1]))[:a.nreps]:
            out = DST / tag / rep.name
            if (out / 'config.yaml').exists():
                continue
            out.mkdir(parents=True, exist_ok=True)
            for item in rep.iterdir():
                if item.name in ('components.yaml', 'run.py'):
                    shutil.copy2(item, out / item.name)
                elif item.name == 'input':
                    shutil.copytree(item, out / 'input', symlinks=False, dirs_exist_ok=True)
            f = out / 'input' / 'custom_restraints.txt'
            lines = [l for l in f.read_text().split('\n') if l.strip()] if f.exists() else []
            n_intra = len(lines)
            for p9, dx, r0, _ in tiers:
                lines.append(f'parp9 1 {p9} | dtx3l 1 {dx} | {r0} {k}')
            f.write_text('\n'.join(lines) + '\n')
            cfg = yaml.safe_load(open(rep / 'config.yaml'))
            cfg.update(custom_restraints=True, custom_restraint_type='harmonic',
                       fcustom_restraints='input/custom_restraints.txt',
                       restart='checkpoint', frestart='restart.chk')
            yaml.dump(cfg, open(out / 'config.yaml', 'w'),
                      default_flow_style=False, sort_keys=False)
            made += 1
        print(f'  {tag}: {a.nreps} reps, {len(tiers)} tiered crosslink restraints '
              f'at k={k} kJ/mol/nm^2 (+{n_intra} intra, unchanged)')
    print(f'\nprepared {made} replicate directories under {DST}')
    print('control for the difference map: unrestrained binding/p9_dtx3l (7 reps, done)')


if __name__ == '__main__':
    main()
