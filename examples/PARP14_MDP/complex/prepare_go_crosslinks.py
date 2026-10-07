#!/usr/bin/env python3
"""Go-restrained PARP9-DTX3L runs driven by the PUBLISHED crosslinks.

WHY
---
CALVADOS puts PARP9-DTX3L at an apparent Kd of 263 uM and spreads the interface
over thousands of equally weak contacts. The experiment (Ashok et al. 2022,
Biochem J 479:289-304, doi 10.1042/BCJ20210722) finds a nanomolar 1:1
heterodimer via DTX3L D3 (230-510) and exactly THREE reproducible BS3
inter-protein crosslinks in the full-length + full-length condition:

    PARP9 K557 - DTX3L K401
    PARP9 K632 - DTX3L K401
    PARP9 K557 - DTX3L K363

Those three pairs are the experimental restraint set. Using them is integrative
modelling: the geometry comes from data, not from a predicted pose (the AF3
PARP9-DTX3L models sit at ipTM 0.30-0.51, too low to restrain against).

A HARD CONSTRAINT IN CALVADOS, FOUND WHILE BUILDING THIS
--------------------------------------------------------
There is ONE global `custom_restraint_type` per run, so harmonic and Go cannot
be mixed. That matters because this set already carries 85 custom restraints,
all intra-PARP9, bridging residues 60-96 to 516-539 -- they hold PARP9's SPLIT
KH1 domain together as one fold. They cannot simply be dropped.

So there are exactly two buildable options, and neither is free:

  HARMONIC inter-chain (primary here)
      The 85 intra restraints stay byte-identical to the control, so the ONLY
      difference from the unrestrained run is the three added crosslinks. Clean
      difference map. Cost: a harmonic restraint is unbounded, so the pair is
      tethered and can never dissociate.

  GO inter-chain (one condition, for contrast)
      Reversible -- finite well depth, so the complex can come apart. Cost: the
      85 intra restraints are reinterpreted as Go too. At k=350 that is a 144 kT
      well whose curvature at the minimum is 116,667 kJ/mol/nm^2 versus the
      harmonic's 350 -- i.e. PARP9's split KH1 becomes ~300x MORE rigid than in
      the control, which is a confound in every downstream comparison.

The harmonic sweep is primary because its comparison to the existing control is
exact. The single Go condition tests whether reversibility changes the answer.

WHY GO IS STILL THE MORE PHYSICAL FORM
--------------------------------------
interactions.init_restraints:
  harmonic -> HarmonicBondForce, k in kJ/mol/nm^2. Quadratic, never releases;
              it WELDS the pair. Bound fraction becomes 1 by construction and
              the observable is destroyed.
  go       -> k*(5*(s/r)^12 - 6*(s/r)^10), k in kJ/mol = WELL DEPTH, minimum at
              r = s, decaying to zero. Finite and reversible, so the complex can
              still dissociate and exchange -- the binding free energy becomes a
              tunable parameter instead of an imposed constraint.

r is set to 2.0 nm, not to a tight contact distance. The data says "within BS3
reach", so the well minimum belongs at BS3 reach; putting it at 0.8 nm would
assert a specific register the three crosslinks do not determine.

k IS SWEPT, because the right value is unknown and is the thing being fitted.
At 293 K, kT = 2.44 kJ/mol:
    k =  5 kJ/mol  ~2 kT   marginal, largely transient
    k = 15 kJ/mol  ~6 kT   the project's k_go default
    k = 50 kJ/mol ~20 kT   strong, near-permanent

WHAT IS AND IS NOT MEASURABLE AFTERWARDS
----------------------------------------
Once restrained, binding is an INPUT. Bound fraction and Kd are no longer
results. What the runs can answer: given that the D3-PARP9 interface is held,
(a) does the rest of the interface concentrate onto the published crosslinks, and
(b) do other contacts appear as a consequence. (b) needs the unrestrained
control -- this project has already been bitten once, when Go-KH restraints
manufactured the dominant MD3-linker contact (+25 contacts/frame under Go,
nothing without it). Always difference-map against the control.

NOTE add_custom_restraints calls add_exclusions by default, so AH and Yukawa are
REMOVED between restrained pairs: the restraint replaces the native interaction
rather than adding to it. Relevant when interpreting the well depth.

Usage:  python prepare_go_crosslinks.py [--nreps 5] [--ks 5 15 50]
"""
import shutil
from pathlib import Path
from argparse import ArgumentParser
import yaml

HERE = Path(__file__).resolve().parent
SRC = HERE / 'binding' / 'p9_dtx3l'
DST = HERE / 'binding_go'

# (parp9 residue, dtx3l residue) -- Ashok 2022, FL PARP9 + FL DTX3L, BS3
XLINKS = [(557, 401), (632, 401), (557, 363)]
R_NM = 2.0

# (tag suffix, restraint type, k). Harmonic k is kJ/mol/nm^2 (stiffness);
# Go k is kJ/mol (well depth). At 293 K kT = 2.44 kJ/mol, so the harmonic
# tethers allow excursions of sqrt(kT/k) = 0.70 / 0.35 / 0.16 nm about r.
CONDITIONS = [('xl_h5', 'harmonic', 5.0),
              ('xl_h20', 'harmonic', 20.0),
              ('xl_h100', 'harmonic', 100.0),
              ('xl_go15', 'go', 15.0)]


def main():
    ap = ArgumentParser()
    ap.add_argument('--nreps', type=int, default=5)
    ap.add_argument('--only', nargs='*', default=None,
                    help='restrict to these condition tags')
    a = ap.parse_args()
    if not SRC.is_dir():
        raise SystemExit(f'{SRC} not found')
    made = 0
    for suffix, rtype, k in CONDITIONS:
        if a.only and suffix not in a.only:
            continue
        tag = f'p9_dtx3l_{suffix}'
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
                    shutil.copytree(item, out / 'input', symlinks=False,
                                    dirs_exist_ok=True)
            # the experimental restraint set, on top of the existing
            # intra-domain custom restraints this set already carries
            existing = (out / 'input' / 'custom_restraints.txt')
            lines = existing.read_text().rstrip('\n').split('\n') if existing.exists() else []
            n_intra = len([l for l in lines if l.strip()])
            for (p9, dx) in XLINKS:
                lines.append(f'parp9 1 {p9} | dtx3l 1 {dx} | {R_NM} {k}')
            existing.write_text('\n'.join(l for l in lines if l.strip()) + '\n')
            cfg = yaml.safe_load(open(rep / 'config.yaml'))
            cfg['custom_restraints'] = True
            cfg['custom_restraint_type'] = rtype
            cfg['fcustom_restraints'] = 'input/custom_restraints.txt'
            cfg['restart'] = 'checkpoint'
            cfg['frestart'] = 'restart.chk'
            yaml.dump(cfg, open(out / 'config.yaml', 'w'),
                      default_flow_style=False, sort_keys=False)
            made += 1
        unit = 'kJ/mol/nm^2' if rtype == 'harmonic' else 'kJ/mol (well depth)'
        print(f'  {tag}: {a.nreps} reps, {len(XLINKS)} inter-chain {rtype} restraints '
              f'at r={R_NM} nm, k={k} {unit}  (+{n_intra} intra, also {rtype})')
    print(f'\nprepared {made} replicate directories under {DST}')
    print('Control for the difference map: the existing unrestrained binding/p9_dtx3l '
          '(7 replicates, already run).')
    print('In the xl_go15 condition the 85 intra-PARP9 split-KH1 restraints are Go too, '
          'making that fold ~300x stiffer than the control -- expected, documented, and '
          'the reason the harmonic conditions are primary.')


if __name__ == '__main__':
    main()
