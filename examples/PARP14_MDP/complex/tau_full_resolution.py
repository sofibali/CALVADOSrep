#!/usr/bin/env python3
"""Contact lifetimes at FULL trajectory resolution (0.05 ns), all sets.

Why this exists: tau measured from the stride-10 stage-1 data is pinned to the
sampling interval, not the physics. Decimating one trajectory shows it directly:

    0.05 ns  113 episodes  tau 0.27 ns
    0.10 ns   79 episodes  tau 0.39 ns
    0.20 ns   51 episodes  tau 0.61 ns
    0.50 ns   37 episodes  tau 0.89 ns   <- what episode_stats.py reported

Every number moves monotonically with the stride, so none of them is converged
except the finest. Bound FRACTION is robust (0.061 -> 0.066 over that range);
tau and the event count are not. This recomputes the bound/unbound series at
stride 1 so the lifetimes are real, and reports the survival-fit lifetime
alongside the naive mean, because a mean of quantised run lengths is biased low
by censoring at the trajectory ends.

Parallel over replicates. ~62 replicates; use -j to match the machine.
"""
import os, sys, warnings, itertools, json
from pathlib import Path
from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor
import numpy as np, yaml
warnings.filterwarnings('ignore')

HERE = Path(__file__).resolve().parent
ROOT = HERE / 'binding'
OUT = HERE / 'analysis'
NS_PER_FRAME = 0.05
BOUND_NM = 1.0
EQ_NS = 25.0
SKIP = int(EQ_NS / NS_PER_FRAME)

sys.path.insert(0, str(HERE))
from analyze_binding import SETS, sysname, topology
import analyze_convergence as ac


def one_rep(args):
    s, rep_name = args
    import MDAnalysis as mda
    from MDAnalysis.analysis.distances import distance_array
    rep = ROOT / s / rep_name
    dcd = rep / f'{sysname(s)}.dcd'
    top = topology(rep)
    if top is None or not dcd.exists():
        return None
    chains = SETS[s]
    cfg = yaml.safe_load(open(rep / 'config.yaml'))
    box = np.array(cfg['box'], float)
    boxv = np.array([box[0]*10, box[1]*10, box[2]*10, 90., 90., 90.], np.float32)
    u = mda.Universe(str(top), str(dcd))
    offs, o = [], 0
    for c in chains:
        offs.append((o, o + ac.LENGTHS[c])); o += ac.LENGTHS[c]
    if u.atoms.n_atoms != o:
        return None
    pairs = list(itertools.combinations(range(len(chains)), 2))
    series = {p: [] for p in pairs}
    for ts in u.trajectory[SKIP:]:          # stride 1
        P = u.atoms.positions.astype(np.float32)
        for p in pairs:
            i, j = p
            D = distance_array(P[offs[i][0]:offs[i][1]], P[offs[j][0]:offs[j][1]], box=boxv)
            series[p].append(D.min() / 10.0)
    return s, rep_name, {f'{chains[p[0]]}-{chains[p[1]]}': np.asarray(v, np.float32)
                         for p, v in series.items()}


def episodes(b):
    d = np.diff(np.concatenate(([0], np.asarray(b, int), [0])))
    st, en = np.where(d == 1)[0], np.where(d == -1)[0]
    return st, en


def survival_tau(durs, censored):
    """Exponential MLE with right-censoring: tau = sum(all durations)/n_uncensored.

    Episodes running past the last frame are censored, not short. Ignoring that
    biases a naive mean low; this is the standard correction.
    """
    durs = np.asarray(durs, float)
    n_unc = int((~np.asarray(censored, bool)).sum())
    return float(durs.sum() / n_unc) if n_unc else np.nan


def main():
    ap = ArgumentParser()
    ap.add_argument('-j', '--workers', type=int, default=16)
    ap.add_argument('sets', nargs='*')
    a = ap.parse_args()
    todo = a.sets or [s for s in SETS if (ROOT / s).is_dir()]
    jobs = []
    for s in todo:
        for rep in sorted((ROOT / s).glob('rep-*'), key=lambda p: int(p.name.split('-')[1])):
            if (rep / f'{sysname(s)}.dcd').exists():
                jobs.append((s, rep.name))
    print(f'{len(jobs)} replicates, {a.workers} workers', flush=True)

    acc = {}
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for k, r in enumerate(ex.map(one_rep, jobs), 1):
            if r is None:
                continue
            s, rep_name, d = r
            for pair, m in d.items():
                acc.setdefault((s, pair), []).append(m)
            print(f'  [{k}/{len(jobs)}] {s}/{rep_name}', flush=True)

    rows = []
    for (s, pair), mats in sorted(acc.items()):
        durs, cens, fbs = [], [], []
        for m in mats:
            b = m < BOUND_NM
            fbs.append(float(b.mean()))
            st, en = episodes(b)
            durs += list((en - st) * NS_PER_FRAME)
            cens += [bool(e == len(b)) for e in en]
        durs = np.asarray(durs); fb = np.asarray(fbs)
        rows.append(dict(
            set=s, pair=pair, n_rep=len(mats),
            ns_per_rep=round(len(mats[0]) * NS_PER_FRAME, 1),
            bound_frac=round(float(fb.mean()), 4),
            bound_frac_ci95=round(float(1.96*fb.std(ddof=1)/np.sqrt(len(fb))), 4) if len(fb) > 1 else np.nan,
            n_episodes=int(durs.size),
            episodes_per_us=round(float(durs.size / (len(mats) * len(mats[0]) * NS_PER_FRAME / 1000)), 1),
            tau_mean_ns=round(float(durs.mean()), 3),
            tau_median_ns=round(float(np.median(durs)), 3),
            tau_survival_ns=round(survival_tau(durs, cens), 3),
            tau_p95_ns=round(float(np.percentile(durs, 95)), 2),
            tau_max_ns=round(float(durs.max()), 2)))
    import pandas as pd
    df = pd.DataFrame(rows).sort_values('bound_frac', ascending=False)
    df.to_csv(OUT / 'tau_full_resolution.csv', index=False)
    np.save(OUT / 'data' / 'mind_full.npy',
            np.array([{k: [np.asarray(x) for x in v] for k, v in acc.items()}], dtype=object),
            allow_pickle=True)
    print(df.to_string(index=False), flush=True)
    print(f"\nsaved {OUT/'tau_full_resolution.csv'}", flush=True)


if __name__ == '__main__':
    main()
