#!/usr/bin/env python3
"""Crosslink satisfaction: % of frames each published crosslink is within threshold.

THE METRIC
----------
The integrative-modelling field's standard score for a structure or ensemble
against XL-MS data: the fraction of restraints satisfied within a distance
threshold. Reporting it lets CALVADOS and HyRes be compared to each other and
to the experiment on the same footing, instead of by eye.

DATA
----
Ashok Y, Vela-Rodriguez C, Yang C, Alanen HI, Liu F, Paschal BM, Lehtio L.
"Reconstitution of the DTX3L-PARP9 complex reveals determinants for
high-affinity heterodimerization and multimeric assembly."
Biochem J. 2022 Feb 11;479(3):289-304. doi:10.1042/BCJ20210722. PMID 35037691.
BS3, four replicates, full-length PARP9 + full-length DTX3L, three
inter-protein crosslinks.

TIERED THRESHOLDS
-----------------
Not all three crosslinks carry the same information. Following IMProv
(Ziemianowicz et al., Mol Cell Proteomics 2021;20:100139,
doi:10.1016/j.mcpro.2021.100139), a crosslink formed in a flexible region can
be kinetically trapped and records a transient excursion rather than the
equilibrium structure, so it deserves a looser threshold than one joining two
stable folds. IMProv's offset scheme uses 25 A between stable regions and 40 A
between unstable ones.

Classifying the three against this project's restrained-domain definitions:

  K557-K401   PARP9 KH2 (restrained)  - DTX3L KH4 (restrained)  RIGID-RIGID
  K632-K401   PARP9 linker            - DTX3L KH4 (restrained)  MIXED
  K557-K363   PARP9 KH2 (restrained)  - DTX3L linker            MIXED

Only K557-K401 constrains the relative orientation of two folded domains. It is
the restraint that should hold in any correct equilibrium structure, and it is
therefore reported separately as well as in the total.

All three thresholds are reported regardless, so nothing hinges on the tiering:
  25 A  IMProv "stable-stable"
  30 A  the conventional BS3 CA-CA cutoff, and this project's LYS_LYS_CUTOFF
  40 A  IMProv "unstable-unstable"

REACTIVITY FILTER
-----------------
BS3 acylates solvent-exposed lysine primary amines; a buried lysine is
unreactive at any distance, so a model that buries one cannot satisfy its
crosslink even in principle. Exposure is reported per lysine as the fraction of
frames it sits below its chain's median CA coordination -- a CG proxy, not SASA.
A lysine exposed in ANY appreciable fraction of frames counts as reactive,
because crosslinking is a slow chemistry that integrates over the ensemble.

Usage:
    python crosslink_satisfaction.py                      # CALVADOS p9_dtx3l
    python crosslink_satisfaction.py --hyres hyres/runs/prod_rep1
"""
import sys, warnings
from pathlib import Path
from argparse import ArgumentParser
import numpy as np, pandas as pd, yaml
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
warnings.filterwarnings('ignore')

HERE = Path(__file__).resolve().parent
OUT = HERE / 'analysis'
OUT.mkdir(exist_ok=True)
sys.path.insert(0, str(HERE))
from analyze_binding import SETS, sysname, topology
import analyze_convergence as ac

# (parp9 resid, dtx3l resid, IMProv class)
XLINKS = [(557, 401, 'rigid-rigid'), (632, 401, 'mixed'), (557, 363, 'mixed')]
THRESH = {'25A (IMProv stable)': 2.5, '30A (BS3 conventional)': 3.0,
          '40A (IMProv flexible)': 4.0}
TIERED = {'rigid-rigid': 2.5, 'mixed': 4.0}
EQ_FRAMES = 500
CONTACT_NM = 1.0


def series_calvados(set_name, stride):
    """Per-crosslink CA-CA distance (nm) and per-lysine exposure, pooled over reps."""
    chains = SETS[set_name]
    i9, ix = chains.index('parp9'), chains.index('dtx3l')
    off, o = {}, 0
    for c in chains:
        off[c] = o; o += ac.LENGTHS[c]
    d = {k: [] for k in XLINKS}
    expo = {}
    for rep in sorted((HERE / 'binding' / set_name).glob('rep-*'),
                      key=lambda p: int(p.name.split('-')[1])):
        dcd = rep / f'{sysname(set_name)}.dcd'; top = topology(rep)
        if top is None or not dcd.exists():
            continue
        cfg = yaml.safe_load(open(rep / 'config.yaml'))
        b = np.array(cfg['box'], float)
        boxv = np.array([b[0]*10, b[1]*10, b[2]*10, 90., 90., 90.], np.float32)
        u = mda.Universe(str(top), str(dcd))
        if u.atoms.n_atoms != o:
            continue
        idx = {(c, r): off[c] + r - 1 for c in ('parp9', 'dtx3l') for r in
               set([x[0] for x in XLINKS] if c == 'parp9' else [x[1] for x in XLINKS])}
        for ts in u.trajectory[EQ_FRAMES::stride]:
            P = u.atoms.positions.astype(np.float32)
            for k in XLINKS:
                a, bb = idx[('parp9', k[0])], idx[('dtx3l', k[1])]
                d[k].append(float(distance_array(P[a:a+1], P[bb:bb+1], box=boxv)[0, 0]) / 10.0)
            for c, rr in (('parp9', [x[0] for x in XLINKS]), ('dtx3l', [x[1] for x in XLINKS])):
                lo, hi = off[c], off[c] + ac.LENGTHS[c]
                sub = P[lo:hi]
                for r in set(rr):
                    j = r - 1
                    nb = int((np.linalg.norm(sub - sub[j], axis=1) < CONTACT_NM * 10).sum())
                    expo.setdefault((c, r), []).append(nb)
    return {k: np.asarray(v) for k, v in d.items()}, expo


def series_hyres(run_dir, stride):
    """Same, from the HyRes trajectory. CA beads carry the residue numbering."""
    run_dir = Path(run_dir)
    u = mda.Universe(str(run_dir / 'complex_start.pdb'), str(run_dir / 'system.dcd'))
    ca = u.select_atoms('name CA')
    n9 = 854
    # chains were concatenated parp9 then dtx3l, each renumbered from 1
    if len(ca) != n9 + 740:
        raise SystemExit(f'expected {n9+740} CA, found {len(ca)}')
    L = 400.0
    boxv = np.array([L, L, L, 90., 90., 90.], np.float32)
    d = {k: [] for k in XLINKS}
    expo = {}
    for ts in u.trajectory[::stride]:
        P = ca.positions.astype(np.float32)
        for k in XLINKS:
            a, b = k[0] - 1, n9 + k[1] - 1
            d[k].append(float(distance_array(P[a:a+1], P[b:b+1], box=boxv)[0, 0]) / 10.0)
        for c, lo, hi, rr in (('parp9', 0, n9, [x[0] for x in XLINKS]),
                              ('dtx3l', n9, len(ca), [x[1] for x in XLINKS])):
            sub = P[lo:hi]
            for r in set(rr):
                nb = int((np.linalg.norm(sub - sub[r-1], axis=1) < CONTACT_NM * 10).sum())
                expo.setdefault((c, r), []).append(nb)
    return {k: np.asarray(v) for k, v in d.items()}, expo


def report(label, d, expo):
    n = len(next(iter(d.values())))
    med = {}
    for (c, r), v in expo.items():
        med.setdefault(c, []).extend(v)
    chain_med = {c: float(np.median(v)) for c, v in med.items()}
    rows = []
    for k in XLINKS:
        p9, dx, cls = k
        v = d[k]
        row = dict(model=label, crosslink=f'PARP9 K{p9}-DTX3L K{dx}', tier=cls,
                   min_nm=round(float(v.min()), 2), median_nm=round(float(np.median(v)), 1))
        for name, t in THRESH.items():
            row[name] = round(100 * float((v <= t).mean()), 2)
        row['tiered'] = round(100 * float((v <= TIERED[cls]).mean()), 2)
        e9 = np.asarray(expo[('parp9', p9)]); ex = np.asarray(expo[('dtx3l', dx)])
        row['K_exposed_%'] = (f"{100*(e9 <= chain_med['parp9']).mean():.0f}/"
                              f"{100*(ex <= chain_med['dtx3l']).mean():.0f}")
        rows.append(row)
    df = pd.DataFrame(rows)
    print(f'\n=== {label} ({n} frames analysed) ===')
    cols = ['crosslink', 'tier', 'min_nm', 'median_nm'] + list(THRESH) + ['tiered', 'K_exposed_%']
    print(df[cols].to_string(index=False))
    print(f"  SATISFACTION at 30 A (mean over 3):  {df['30A (BS3 conventional)'].mean():6.2f}% of frames")
    print(f"  SATISFACTION tiered (IMProv):        {df['tiered'].mean():6.2f}% of frames")
    rr = df[df.tier == 'rigid-rigid']
    print(f"  the one RIGID-RIGID restraint at 25 A: {rr['25A (IMProv stable)'].iloc[0]:6.2f}% of frames")
    return df


def main():
    ap = ArgumentParser()
    ap.add_argument('--set', default='p9_dtx3l')
    ap.add_argument('--hyres', default=None)
    ap.add_argument('--stride', type=int, default=10)
    a = ap.parse_args()
    out = []
    d, e = series_calvados(a.set, a.stride)
    out.append(report(f'CALVADOS {a.set}', d, e))
    if a.hyres:
        d, e = series_hyres(a.hyres, a.stride)
        out.append(report('HyRes (1 replicate, 500 ns)', d, e))
    pd.concat(out, ignore_index=True).to_csv(OUT / 'crosslink_satisfaction.csv', index=False)
    print(f'\nsaved {OUT/"crosslink_satisfaction.csv"}')


if __name__ == '__main__':
    main()
