#!/usr/bin/env python3
"""Predicted INTER-CHAIN crosslinks for XL-MS on pulled-down complexes.

MATCHED TO THE EXPERIMENT
-------------------------
The readout available is crosslinking + MS on protein pulled out of THP-1
macrophages polarised to M1 or M2. That is not a purified-protein affinity
measurement, so the useful prediction is not a Kd -- it is *which inter-chain
residue pairs are close enough to crosslink, and how often*. Those are peptide
pairs that either show up in the MS or do not.

CUTOFF, following the project convention in sim_analysis/analyze_lys_contacts.py
  LYS_LYS_CUTOFF     = 3.0 nm CA-CA  (DSS / BS3 / DSSO)
  LYS_ACIDIC_CUTOFF  = 3.0 nm CA-CA  (EDC / DMTMM, Lys-Asp/Glu)
CALVADOS is CA-only, so CA-CA is the native observable -- no side-chain
reconstruction is needed or possible. Note this is a Euclidean cutoff; the
single-chain pipeline additionally filters by solvent-accessible surface
distance (SASD), which removes pairs whose straight line passes through protein.
That correction is NOT applied here, so treat these as an upper bound on
feasible pairs.

TWO PERSISTENCE NUMBERS, AND WHY BOTH MATTER
  persist_all    fraction of ALL frames the pair is within the cutoff.
                 This is what crosslink YIELD scales with. Since these complexes
                 are bound only 6-18% of the time, absolute yields should be low
                 -- a prediction in its own right, and a warning that absence of
                 a crosslink is weak evidence here.
  persist_bound  fraction of BOUND frames. This is the interface register: which
                 pairs are close *given* that the chains are together. Use this
                 to rank candidate peptide pairs; use persist_all to set
                 expectations about detecting them.

Usage:  python predict_interchain_crosslinks.py [sets ...] [--target-frames 1000]
"""
import sys, warnings, itertools
from pathlib import Path
from argparse import ArgumentParser
import numpy as np, pandas as pd, yaml
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array, self_distance_array
warnings.filterwarnings('ignore')

HERE = Path(__file__).resolve().parent
ROOT = HERE / 'binding'
OUT = HERE / 'analysis'
OUT.mkdir(exist_ok=True)
sys.path.insert(0, str(HERE))
from analyze_binding import SETS, sysname, topology
import analyze_convergence as ac

LYS_LYS_NM = 3.0
LYS_ACIDIC_NM = 3.0
BOUND_NM = 1.0
EQ_FRAMES = 500
DEFAULT_SETS = ['p9_dtx3l', 'p14_dtx3l', 'p14_p9', 'ternary', 'p14_homo']


def run_set(s, target_frames):
    chains = SETS[s]
    sysn = sysname(s)
    acc, nall, nbound = {}, 0, 0
    for rep in sorted((ROOT / s).glob('rep-*'), key=lambda p: int(p.name.split('-')[1])):
        dcd = rep / f'{sysn}.dcd'; top = topology(rep)
        if top is None or not dcd.exists():
            continue
        cfg = yaml.safe_load(open(rep / 'config.yaml'))
        box = np.array(cfg['box'], float)
        boxv = np.array([box[0]*10, box[1]*10, box[2]*10, 90., 90., 90.], np.float32)
        u = mda.Universe(str(top), str(dcd))
        offs, o = [], 0
        for c in chains:
            offs.append((o, o + ac.LENGTHS[c])); o += ac.LENGTHS[c]
        if u.atoms.n_atoms != o:
            continue
        names = u.atoms.resnames; resids = u.atoms.resids
        # reactive-group indices per chain
        sel = []
        for (a, b) in offs:
            k = np.where(names[a:b] == 'LYS')[0] + a
            ac_ = np.where(np.isin(names[a:b], ('ASP', 'GLU')))[0] + a
            sel.append((k, ac_))
        cpairs = list(itertools.combinations(range(len(chains)), 2))
        usable = max(0, len(u.trajectory) - EQ_FRAMES)
        step = max(1, usable // target_frames) if target_frames else 1
        for ts in u.trajectory[EQ_FRAMES::step]:
            P = u.atoms.positions.astype(np.float32)
            gmin = np.inf
            for (i, j) in cpairs:
                D = distance_array(P[offs[i][0]:offs[i][1]],
                                   P[offs[j][0]:offs[j][1]], box=boxv) / 10.0
                gmin = min(gmin, float(D.min()))
            isb = gmin < BOUND_NM
            nall += 1; nbound += int(isb)
            for (i, j) in cpairs:
                for kind, (ai, aj), cut in (
                        ('DSS K-K', (sel[i][0], sel[j][0]), LYS_LYS_NM),
                        ('EDC K-acidic', (sel[i][0], sel[j][1]), LYS_ACIDIC_NM),
                        ('EDC acidic-K', (sel[i][1], sel[j][0]), LYS_ACIDIC_NM)):
                    if not len(ai) or not len(aj):
                        continue
                    D = distance_array(P[ai], P[aj], box=boxv) / 10.0
                    xi, yj = np.nonzero(D < cut)
                    for x, y in zip(xi, yj):
                        key = (kind, chains[i], int(resids[ai[x]]),
                               chains[j], int(resids[aj[y]]))
                        v = acc.setdefault(key, [0, 0])
                        v[0] += 1
                        v[1] += int(isb)
    if not nall:
        return None
    rows = []
    for (kind, ci, ri, cj, rj), (na, nb) in acc.items():
        rows.append(dict(set=s, linker=kind, chain_a=ci, resid_a=ri,
                         chain_b=cj, resid_b=rj,
                         persist_all=round(na / nall, 5),
                         persist_bound=round(nb / max(nbound, 1), 4),
                         n_frames=nall, n_bound=nbound))
    return pd.DataFrame(rows)


def main():
    ap = ArgumentParser()
    ap.add_argument('sets', nargs='*', default=DEFAULT_SETS)
    ap.add_argument('--target-frames', type=int, default=1000)
    ap.add_argument('--top', type=int, default=15)
    a = ap.parse_args()
    all_df = []
    for s in a.sets:
        if not (ROOT / s).is_dir():
            continue
        print(f'== {s}', flush=True)
        d = run_set(s, a.target_frames)
        if d is None or not len(d):
            print('   no crosslinkable pairs'); continue
        all_df.append(d)
        top = d.sort_values('persist_bound', ascending=False).head(a.top)
        print(f'   {len(d)} feasible pairs; bound in {d.n_bound.iloc[0]}/{d.n_frames.iloc[0]} frames')
        print(top[['linker', 'chain_a', 'resid_a', 'chain_b', 'resid_b',
                   'persist_bound', 'persist_all']].to_string(index=False))
    if all_df:
        df = pd.concat(all_df, ignore_index=True)
        df.to_csv(OUT / 'interchain_crosslinks.csv', index=False)
        print(f'\nsaved {OUT/"interchain_crosslinks.csv"}  ({len(df)} pairs)')


if __name__ == '__main__':
    main()
