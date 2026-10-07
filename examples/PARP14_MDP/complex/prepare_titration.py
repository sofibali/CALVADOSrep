#!/usr/bin/env python3
"""Concentration series: the in-silico dilution experiment.

WHY THIS AND NOT A FORCE-FIELD NULL
-----------------------------------
Switching off the Ashbaugh-Hatch attraction tests whether CALVADOS's potential
produces the contacts. That is a statement about the model, and there is no
bench experiment corresponding to it. A dilution series is the opposite: it is
exactly what an ITC/SPR/AUC titration does, so the result is comparable to
something measurable.

THE TEST
--------
One chain of each type in a cubic box of side L is a concentration
C = 1/(N_A L^3). For a 1:1 equilibrium the apparent dissociation constant is

    Kd = C (1 - f)^2 / f              f = bound fraction

If the association is a genuine thermodynamic equilibrium, Kd is INVARIANT
across box sizes: dilution lowers f in exactly the way mass action predicts.
If Kd drifts systematically with L, what is being measured is confinement --
two chains in a small box collide because they have nowhere else to be -- and
the "bound fraction" is an artefact of the box, not a property of the pair.

This also fixes the weakness found in the production runs: Kd there varies
2.6-6.4x with the contact cutoff, so a single box and a single cutoff cannot
distinguish a complex from a collision. Invariance across concentration is a
much stronger claim than a plateau in cutoff.

BOX FLOOR
---------
L must exceed the chain's own span plus cutoff_yu (4.0 nm), or a chain sees its
own periodic image. 99th-percentile spans: PARP14 32.0, DTX3L 26.8, PARP9 20.8
nm -> floor 31 nm for parp9/dtx3l pairs, 36 nm for anything containing PARP14.

Separated arm only: the docked sets restart from fixed coordinates
(restart: pdb), so rescaling their box would not rebuild the start state.

Usage:  python prepare_titration.py [--nreps 5]
"""
import shutil, math
from pathlib import Path
from argparse import ArgumentParser
import yaml

HERE = Path(__file__).resolve().parent
SRC = HERE / 'binding'
DST = HERE / 'binding_titration'
NA = 6.02214076e23

# set -> box sizes (nm). 40 already exists in the production runs.
PLAN = {
    'p9_dtx3l':  [35, 50, 60, 70],   # floor 31 nm
    'p14_dtx3l': [50, 70],           # floor 36 nm
}


def conc_uM(L_nm):
    return 1.0 / (NA * (L_nm * 1e-9) ** 3 * 1000) * 1e6


def main():
    ap = ArgumentParser()
    ap.add_argument('--nreps', type=int, default=5)
    a = ap.parse_args()
    made = 0
    print(f'{"set":<12}{"box nm":>8}{"conc uM":>10}  reps')
    for s, boxes in PLAN.items():
        sd = SRC / s
        if not sd.is_dir():
            print(f'  skip {s}'); continue
        for L in boxes:
            tag = f'{s}_box{L}'
            n = 0
            for rep in sorted(sd.glob('rep-*'),
                              key=lambda p: int(p.name.split('-')[1]))[:a.nreps]:
                out = DST / tag / rep.name
                if (out / 'config.yaml').exists():
                    n += 1; continue
                out.mkdir(parents=True, exist_ok=True)
                for item in rep.iterdir():
                    if item.name in ('components.yaml', 'run.py'):
                        shutil.copy2(item, out / item.name)
                    elif item.name == 'input':
                        shutil.copytree(item, out / 'input', symlinks=False,
                                        dirs_exist_ok=True)
                cfg = yaml.safe_load(open(rep / 'config.yaml'))
                cfg['box'] = [float(L), float(L), float(L)]
                # fresh start: grid placement is rebuilt from the new box
                cfg['restart'] = 'checkpoint'
                cfg['frestart'] = 'restart.chk'
                yaml.dump(cfg, open(out / 'config.yaml', 'w'),
                          default_flow_style=False, sort_keys=False)
                n += 1; made += 1
            print(f'{s:<12}{L:>8}{conc_uM(L):>10.1f}  {n}')
    print(f'\nprepared {made} replicate directories under {DST}')
    print('production runs sit at 40 nm = %.1f uM' % conc_uM(40))


if __name__ == '__main__':
    main()
