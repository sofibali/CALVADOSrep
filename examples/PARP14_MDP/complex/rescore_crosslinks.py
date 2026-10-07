#!/usr/bin/env python3
"""Rescore the simulations against ALL published inter-protein crosslinks.

WHAT CHANGED FROM THE FIRST SCORING
-----------------------------------
The first pass used 3 crosslinks. That was an extraction error: the three
Summary sheets of the supplement use different column layouts, so two of them
silently returned nothing. There are 51 inter-protein crosslinks:

    FL PARP9 + D3        26     DTX3L truncated to 230-510
    FL PARP9 + D3RD      22     D3 + RING + DTC
    FL PARP9 + FL DTX3L   3     both full length

and the per-replicate 'x' marks now give reproducibility: 30 of 51 are seen in
>=2 of the 4 BS3 replicates.

The truncated constructs are the informative ones. Cutting DTX3L back to D3
removes the flexible N-terminus that otherwise spreads the crosslinking, so the
D3-PARP9 interface -- the nanomolar one -- is sampled far more densely. Those
crosslinks are still valid targets for a full-length simulation: the same
interface is expected to form, just diluted. What the truncation changes is
detection, not the structure being reported.

TWO METRICS, because they answer different questions
----------------------------------------------------
per-frame     mean over crosslinks of (% of frames within threshold).
              How much of the ensemble looks like the experiment.
ensemble      % of crosslinks satisfied in at least one frame (and, more
              robustly, in at least 1% of frames). The integrative-modelling
              question: can this ensemble explain the restraint at all?

A model can score near zero on the first and high on the second -- that is a
model that visits the right geometry rarely. The opposite is a model stuck in
one wrong place.

Usage:  python rescore_crosslinks.py --hyres hyres/runs/prod_rep1
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
XL = HERE / 'data' / 'published_xlinks_v2.csv'
sys.path.insert(0, str(HERE))
from analyze_binding import SETS, sysname, topology
import analyze_convergence as ac

THRESH = {'25A': 2.5, '30A': 3.0, '40A': 4.0}
EQ = 500


def load_xl():
    d = pd.read_csv(XL)
    d = d[d.block == 'P9-Dtx3L'].copy()
    d['repro_2of4'] = d.repro_2of4.fillna(False).astype(bool)
    return d.reset_index(drop=True)


def dists_calvados(set_name, pairs, stride, root='binding'):
    base = set_name if set_name in SETS else 'p9_dtx3l'
    chains = SETS[base]; off, o = {}, 0
    for c in chains: off[c] = o; o += ac.LENGTHS[c]
    acc = [[] for _ in pairs]
    for rep in sorted((HERE/root/set_name).glob('rep-*'),
                      key=lambda p: int(p.name.split('-')[1])):
        dcd = rep/f'{sysname(base)}.dcd'; top = topology(rep)
        if top is None or not dcd.exists(): continue
        cfg = yaml.safe_load(open(rep/'config.yaml')); b = np.array(cfg['box'], float)
        boxv = np.array([b[0]*10,b[1]*10,b[2]*10,90.,90.,90.], np.float32)
        u = mda.Universe(str(top), str(dcd))
        if u.atoms.n_atoms != o: continue
        ia = np.array([off['parp9']+a-1 for a,_ in pairs])
        ib = np.array([off['dtx3l']+b2-1 for _,b2 in pairs])
        for ts in u.trajectory[EQ::stride]:
            P = u.atoms.positions.astype(np.float32)
            d = np.array([distance_array(P[i:i+1], P[j:j+1], box=boxv)[0,0]
                          for i, j in zip(ia, ib)])/10.0
            for k, v in enumerate(d): acc[k].append(v)
    return [np.asarray(a) for a in acc]


def dists_hyres(run_dir, pairs, stride):
    run = Path(run_dir)
    top = run/'complex_start.pdb'
    if not top.exists():
        top = run/'start.pdb'          # the matched runs name it start.pdb
    u = mda.Universe(str(top), str(run/'system.dcd'))
    ca = u.select_atoms('name CA'); n9 = 854
    L = 400.0; boxv = np.array([L,L,L,90.,90.,90.], np.float32)
    ia = np.array([a-1 for a,_ in pairs]); ib = np.array([n9+b-1 for _,b in pairs])
    acc = [[] for _ in pairs]
    for ts in u.trajectory[::stride]:
        P = ca.positions.astype(np.float32)
        for k,(i,j) in enumerate(zip(ia, ib)):
            acc[k].append(float(distance_array(P[i:i+1], P[j:j+1], box=boxv)[0,0])/10.0)
    return [np.asarray(a) for a in acc]


def score(label, xl, dd):
    rows = []
    for (_, r), v in zip(xl.iterrows(), dd):
        row = dict(model=label, construct=r.construct, n_samples=r.n_samples,
                   repro=r.repro_2of4, xl=f'K{r.res_a}-K{r.res_b}',
                   min_nm=round(float(v.min()), 2))
        for n, t in THRESH.items():
            row[f'pct_{n}'] = round(100*float((v <= t).mean()), 3)
        rows.append(row)
    df = pd.DataFrame(rows)
    def block(sub, name):
        if not len(sub): return None
        return dict(subset=name, n=len(sub),
                    per_frame_30A=round(sub.pct_30A.mean(), 3),
                    ens_ever_30A=round(100*(sub.pct_30A > 0).mean(), 1),
                    ens_1pct_30A=round(100*(sub.pct_30A >= 1.0).mean(), 1),
                    ens_ever_40A=round(100*(sub.pct_40A > 0).mean(), 1))
    summ = [block(df, 'ALL inter-protein'),
            block(df[df.repro], 'reproducible (>=2/4)'),
            block(df[df.n_samples == 4], 'all 4 replicates')]
    for c in df.construct.unique():
        summ.append(block(df[df.construct == c], c))
    s = pd.DataFrame([x for x in summ if x])
    print(f'\n########## {label} ##########')
    print(s.to_string(index=False))
    return df


def main():
    ap = ArgumentParser()
    ap.add_argument('--set', default='p9_dtx3l')
    ap.add_argument('--root', default='binding')
    ap.add_argument('--label', default=None)
    ap.add_argument('--hyres', default=None, nargs='*')
    ap.add_argument('--stride', type=int, default=20)
    a = ap.parse_args()
    xl = load_xl()
    pairs = list(zip(xl.res_a, xl.res_b))
    print(f'{len(pairs)} inter-protein crosslinks ({xl.repro_2of4.sum()} reproducible)')
    allr = []
    allr.append(score(a.label or f'CALVADOS {a.set}', xl, dists_calvados(a.set, pairs, a.stride, a.root)))
    if a.hyres:
        pooled = None
        for h in a.hyres:
            dh = dists_hyres(h, pairs, a.stride)
            pooled = dh if pooled is None else [np.concatenate([x, y]) for x, y in zip(pooled, dh)]
        allr.append(score(f'HyRes ({len(a.hyres)} reps)', xl, pooled))
    out = pd.concat(allr, ignore_index=True)
    out.to_csv(OUT/'crosslink_rescore.csv', index=False)
    print(f'\nsaved {OUT/"crosslink_rescore.csv"}')
    best = out[out.model.str.startswith('CALVADOS')].nlargest(8, 'pct_30A')
    print('\nCALVADOS best-satisfied crosslinks:')
    print(best[['xl','construct','n_samples','min_nm','pct_30A','pct_40A']].to_string(index=False))


if __name__ == '__main__':
    main()
