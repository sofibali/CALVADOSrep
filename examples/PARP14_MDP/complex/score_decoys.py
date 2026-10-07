#!/usr/bin/env python3
"""Dimension-controlled crosslink scoring: true crosslinks against lysine decoys.

WHY THIS EXISTS. Every score in this campaign is a distance threshold, and the
two engines disagree about chain dimensions as much as about the interface:
from the same start CALVADOS expands DTX3L's Rg to 6.39 nm while HyRes holds it
at 3.13 nm. A chain half the size puts every lysine closer to every other one,
so it satisfies distance restraints better whether or not its interface is
right. "HyRes 88% vs CALVADOS 0.08%" therefore cannot be read as interface
recognition on its own.

THE CONTROL. Score every inter-chain lysine-lysine pair the model contains, not
just the published ones, and ask where the published ones fall WITHIN THAT
MODEL'S OWN distribution. BS3 reaches any solvent-exposed lysine pair, so the
unobserved pairs are the natural decoy set. Because the comparison is a rank
inside each model, overall compaction cancels:

    AUC ~ 0.5   the published crosslinks are no closer than random lysine
                pairs in this model -- it is compact, or wrong, but it is not
                recognising the interface
    AUC -> 1.0  the published crosslinks are specifically the close ones

This is the standard control in integrative XL-MS modelling, and it is the
only number here that compares the two force fields on equal footing.

The per-pair statistic is the same one the headline table uses: the fraction of
frames within the threshold. Pairs are labelled from data/published_xlinks_v2.csv.

Usage:
    python score_decoys.py --root binding --set p9_dtx3l
    python score_decoys.py --hyres hyres/runs/matched/docked_rep1 --label "..."
    python score_decoys.py --all          # every model in the campaign
"""
import sys, warnings
from pathlib import Path
from argparse import ArgumentParser
import numpy as np, pandas as pd, yaml
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
warnings.filterwarnings('ignore')

HERE = Path(__file__).resolve().parent
OUT = HERE / 'analysis'
sys.path.insert(0, str(HERE))
from analyze_binding import SETS, sysname, topology
import analyze_convergence as ac

XL_CSV = HERE / 'data' / 'published_xlinks_v2.csv'
THRESH_A = 30.0                      # the project's LYS_LYS_CUTOFF
EQ_FRAMES = 500                      # 25 ns equilibration, raw dcd frames
HYRES_BOX_A = 400.0


def lysines(chain):
    """resid of every LYS in a chain, from its input structure."""
    pdb = {'parp9': 'binding/p9_dtx3l/input/parp9.pdb',
           'dtx3l': 'binding/p9_dtx3l/input/dtx3l.pdb',
           'parp14': 'binding/p14_p9/input/parp14.pdb'}[chain]
    u = mda.Universe(str(HERE / pdb))
    return np.sort(u.select_atoms('name CA and resname LYS').resids)


def published(a, b):
    """{(resid_a, resid_b)} and the reproducible subset, for chain pair a-b."""
    if not XL_CSV.is_file():
        raise SystemExit(f'missing {XL_CSV}')
    d = pd.read_csv(XL_CSV)
    block = {('parp9', 'dtx3l'): 'P9-Dtx3L'}.get((a, b))
    if block is None:                       # PARP14 pairs: not yet available
        return set(), set()
    d = d[d.block == block]
    allp = {(int(r.res_a), int(r.res_b)) for _, r in d.iterrows()}
    rep = {(int(r.res_a), int(r.res_b)) for _, r in d.iterrows() if r.repro_2of4}
    return allp, rep


def _tally(P, boxv, ia, ib, hit, n):
    D = distance_array(P[ia], P[ib], box=boxv)
    hit += (D < THRESH_A)
    return n + 1


def collect(label, frames_iter, chains, target_frames):
    """frames_iter yields (positions_A, box_vector, chain_offsets)."""
    a, b = chains
    ka, kb = lysines(a), lysines(b)
    hit = np.zeros((len(ka), len(kb)), dtype=np.int32)
    n = 0
    for P, boxv, off in frames_iter:
        ia = off[a] + ka - 1
        ib = off[b] + kb - 1
        n = _tally(P, boxv, ia, ib, hit, n)
    if not n:
        return None
    pct = 100.0 * hit / n
    allp, rep = published(a, b)
    rows = []
    for i, ra in enumerate(ka):
        for j, rb in enumerate(kb):
            rows.append(dict(model=label, chain_a=a, chain_b=b,
                             res_a=int(ra), res_b=int(rb),
                             pct_within=round(float(pct[i, j]), 3),
                             is_xl=(int(ra), int(rb)) in allp,
                             is_repro=(int(ra), int(rb)) in rep))
    return pd.DataFrame(rows), n


def auc(df, col='is_xl'):
    """Rank-based AUC of true crosslinks vs decoys. No sklearn dependency."""
    y = df[col].to_numpy(bool)
    if y.sum() == 0 or (~y).sum() == 0:
        return np.nan
    r = pd.Series(df.pct_within).rank().to_numpy()
    n1, n0 = y.sum(), (~y).sum()
    return float((r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def frames_calvados(s, root, target_frames):
    base = s if s in SETS else 'p9_dtx3l'
    chains = SETS[base]
    off, o = {}, 0
    for c in chains:
        off[c] = o; o += ac.LENGTHS[c]
    step = None
    for rep in sorted((HERE / root / s).glob('rep-*'),
                      key=lambda p: int(p.name.split('-')[1])):
        dcd = rep / f'{sysname(base)}.dcd'
        top = topology(rep)
        if top is None or not dcd.exists():
            continue
        cfg = yaml.safe_load(open(rep / 'config.yaml'))
        bx = np.array(cfg['box'], float)
        boxv = np.array([bx[0]*10, bx[1]*10, bx[2]*10, 90., 90., 90.], np.float32)
        u = mda.Universe(str(top), str(dcd))
        if u.atoms.n_atoms != o:
            continue
        usable = max(0, len(u.trajectory) - EQ_FRAMES)
        if usable < 10:
            continue
        if step is None:
            step = max(1, usable // target_frames) if target_frames else 1
        for ts in u.trajectory[EQ_FRAMES::step]:
            yield u.atoms.positions.astype(np.float32), boxv, off


def frames_hyres(run_dirs, target_frames):
    off = {'parp9': 0, 'dtx3l': ac.LENGTHS['parp9']}
    boxv = np.array([HYRES_BOX_A]*3 + [90., 90., 90.], np.float32)
    for rd in run_dirs:
        run = HERE / rd
        top = run / 'complex_start.pdb'
        if not top.exists():
            top = run / 'start.pdb'
        if not top.exists() or not (run / 'system.dcd').exists():
            continue
        u = mda.Universe(str(top), str(run / 'system.dcd'))
        ca = u.select_atoms('name CA')
        if ca.n_atoms != ac.LENGTHS['parp9'] + ac.LENGTHS['dtx3l']:
            continue
        step = max(1, len(u.trajectory) // target_frames) if target_frames else 1
        for ts in u.trajectory[::step]:
            yield ca.positions.astype(np.float32), boxv, off


JOBS = [('CALVADOS unrestrained (separated)', dict(root='binding', s='p9_dtx3l')),
        ('CALVADOS unrestrained (docked)',    dict(root='binding', s='p9_dtx3l_docked')),
        ('CALVADOS ternary (docked)',         dict(root='binding', s='ternary_docked')),
        ('CALVADOS +3 XL harmonic k=20',      dict(root='binding_go', s='p9_dtx3l_xl_h20')),
        ('CALVADOS +3 XL Go r=1.5 k=100',     dict(root='binding_gotight', s='p9_dtx3l_got_r15_k100')),
        ('CALVADOS +21 XL IMProv-tiered',     dict(root='binding_xl', s='p9_dtx3l_xl30_k20')),
        ('HyRes unrestrained',                dict(hyres=['hyres/runs/prod_rep1'])),
        ('HyRes matched (separated)',         dict(hyres=[f'hyres/runs/matched/separated_rep{i}' for i in (1,2,3)])),
        ('HyRes matched (docked)',            dict(hyres=[f'hyres/runs/matched/docked_rep{i}' for i in (1,2,3)]))]


def run_one(label, spec, target_frames):
    if 'hyres' in spec:
        it = frames_hyres(spec['hyres'], target_frames)
        chains = ('parp9', 'dtx3l')
    else:
        it = frames_calvados(spec['s'], spec['root'], target_frames)
        base = spec['s'] if spec['s'] in SETS else 'p9_dtx3l'
        cs = SETS[base]
        chains = ('parp9', 'dtx3l') if {'parp9','dtx3l'} <= set(cs) else None
        if chains is None:
            print(f'   {label}: no PARP9+DTX3L pair, skipped'); return None
    got = collect(label, it, chains, target_frames)
    if got is None:
        print(f'   {label}: no usable frames'); return None
    df, n = got
    a_all, a_rep = auc(df, 'is_xl'), auc(df, 'is_repro')
    xl, dec = df[df.is_xl], df[~df.is_xl]
    print(f'   {label:34s} n={n:5d} frames  AUC(all)={a_all:.3f}  '
          f'AUC(repro)={a_rep:.3f}  XL median {xl.pct_within.median():6.2f}%  '
          f'decoy median {dec.pct_within.median():6.2f}%')
    return df, dict(model=label, n_frames=n, n_xl=int(df.is_xl.sum()),
                    n_decoy=int((~df.is_xl).sum()),
                    auc_all=round(a_all, 4), auc_repro=round(a_rep, 4),
                    xl_median_pct=round(float(xl.pct_within.median()), 3),
                    decoy_median_pct=round(float(dec.pct_within.median()), 3),
                    xl_mean_pct=round(float(xl.pct_within.mean()), 3),
                    decoy_mean_pct=round(float(dec.pct_within.mean()), 3))


def main():
    ap = ArgumentParser()
    ap.add_argument('--root', default='binding')
    ap.add_argument('--set', dest='s', default=None)
    ap.add_argument('--hyres', nargs='*', default=None)
    ap.add_argument('--label', default=None)
    ap.add_argument('--all', action='store_true')
    ap.add_argument('--target-frames', type=int, default=600)
    a = ap.parse_args()

    jobs = JOBS if a.all else [(a.label or (a.s or 'hyres'),
                                dict(hyres=a.hyres) if a.hyres
                                else dict(root=a.root, s=a.s))]
    pairs, summ = [], []
    print(f'threshold {THRESH_A:.0f} A, decoys = all other inter-chain Lys-Lys pairs')
    for label, spec in jobs:
        r = run_one(label, spec, a.target_frames)
        if r:
            pairs.append(r[0]); summ.append(r[1])
    if not summ:
        return
    OUT.mkdir(exist_ok=True)
    pd.concat(pairs, ignore_index=True).to_csv(OUT / 'decoy_pairs.csv', index=False)
    s = pd.DataFrame(summ)
    s.to_csv(OUT / 'decoy_auc.csv', index=False)
    print('\n' + s.to_string(index=False))
    print('\nwrote', OUT / 'decoy_auc.csv', 'and', OUT / 'decoy_pairs.csv')


if __name__ == '__main__':
    main()
