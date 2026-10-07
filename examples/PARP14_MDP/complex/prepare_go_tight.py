#!/usr/bin/env python3
"""Go-restrained runs on the 3 FL+FL crosslinks, at tighter target distances.

WHY THESE THREE AND NOT THE 21
------------------------------
The 21-crosslink run (binding_xl) failed: 1.9% agreement overall, 0.008% on the
reproducible subset, and only 1 of its 21 restrained pairs ever reached within
3.0 nm -- most sat 7-12 nm from target while being pulled. The restraints were
applied correctly (106 parsed, no errors, chains bound at 0.52 nm); they were
mutually unsatisfiable.

The reason is that those 21 pool crosslinks from THREE constructs. In the D3
construct DTX3L is residues 230-510 only, so PARP9 can reach D3 surfaces that
full-length DTX3L occludes with its N-terminal region. Pooling them asks for a
geometry no single full-length arrangement can provide. The 3 FL+FL crosslinks,
by contrast, were measured on the same construct being simulated, and restraining
them alone reached 30.4% while pulling 43 unrestrained crosslinks from 0.12% to
17.6%.

WHY TIGHTER
-----------
The earlier runs put the minimum at r = 2.0 nm, and the harmonic has no flat
bottom, so it pulls TO 2.0 nm -- only just inside the 3.0 nm scoring cutoff. A
genuinely crosslinked lysine pair sits well below the BS3 ceiling. r = 1.2 and
1.5 nm target the interior of the allowed range instead of its edge.

WHY GO
------
The Go form k*(5(s/r)^12 - 6(s/r)^10) has a finite well depth and a decaying
tail, so the pair can still dissociate -- binding free energy stays a tunable
parameter rather than an imposed constraint, which harmonic cannot offer. The
one previous Go condition used k = 15 kJ/mol (6 kT at 293 K) and reached 17.1%,
below harmonic k=20 at 30.4%, so the well was probably too shallow. This sweeps
deeper: 50 kJ/mol (20 kT) and 100 kJ/mol (41 kT).

CARRIED-OVER CONFOUND, unchanged from the earlier Go run so the two compare:
CALVADOS has one global custom_restraint_type, so the 85 pre-existing
intra-PARP9 split-KH1 restraints become Go too. At k = 350 that is a 144 kT well
whose curvature at the minimum is ~116,700 kJ/mol/nm^2, against the harmonic's
350 -- PARP9's split KH1 is roughly 300x more rigid than in the harmonic runs.

Control for the difference map: unrestrained binding/p9_dtx3l (7 reps, done).

Usage:  python prepare_go_tight.py
"""
import shutil
from pathlib import Path
from argparse import ArgumentParser
import yaml

HERE = Path(__file__).resolve().parent
SRC = HERE / 'binding' / 'p9_dtx3l'
DST = HERE / 'binding_gotight'
XLINKS = [(557, 401), (632, 401), (557, 363)]     # Ashok 2022, FL PARP9 + FL DTX3L
CONDITIONS = [(1.5, 50.0), (1.5, 100.0), (1.2, 50.0), (1.2, 100.0)]
kT = 8.314e-3 * 293


def main():
    ap = ArgumentParser(); ap.add_argument('--nreps', type=int, default=5)
    a = ap.parse_args()
    made = 0
    for r0, k in CONDITIONS:
        tag = f'p9_dtx3l_got_r{str(r0).replace(".","")}_k{int(k)}'
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
            for p9, dx in XLINKS:
                lines.append(f'parp9 1 {p9} | dtx3l 1 {dx} | {r0} {k}')
            f.write_text('\n'.join(lines) + '\n')
            cfg = yaml.safe_load(open(rep / 'config.yaml'))
            cfg.update(custom_restraints=True, custom_restraint_type='go',
                       fcustom_restraints='input/custom_restraints.txt',
                       restart='checkpoint', frestart='restart.chk')
            yaml.dump(cfg, open(out / 'config.yaml', 'w'),
                      default_flow_style=False, sort_keys=False)
            made += 1
        print(f'  {tag}: {a.nreps} reps, 3 inter-chain Go restraints '
              f'r={r0} nm, depth {k} kJ/mol ({k/kT:.0f} kT)  (+{n_intra} intra, also Go)')
    print(f'\nprepared {made} replicate directories under {DST}')


if __name__ == '__main__':
    main()
