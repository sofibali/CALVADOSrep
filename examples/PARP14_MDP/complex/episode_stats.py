#!/usr/bin/env python3
"""Per-set binding-episode statistics, and which sets deserve deeper sampling.

The campaign's original follow-up rule was "extend replicates that haven't
dissociated". Measured on p9_dtx3l that rule never fires: across 1500 ns the
chains associated and dissociated ~35 times per replicate and the longest single
contact was 6 ns. No replicate is ever "still bound".

The observable that does discriminate is the episode lifetime tau. A set with
tau ~ 1 ns is a diffusive encounter complex -- more sampling just buys more
glancing collisions. A set with tau >> 1 ns is a real complex, and there more
sampling is worth it, because the bound-window contact map is built only from
bound frames and those are the scarce resource.

Reads analysis/data/raw.npy (written by analyze_binding.py), so it costs nothing.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd

OUT = Path(__file__).resolve().parent / 'analysis'
NS = 0.5           # ns per analysed sample
BOUND_NM = 1.0
TAU_REAL = 10.0    # ns; above this the contact is more than a diffusive touch
FRAC_REAL = 0.20   # bound fraction worth building a contact map from


def episodes(m):
    """Contiguous bound-episode durations (ns) for one replicate's min-distance series."""
    b = np.asarray(m) < BOUND_NM
    d = np.diff(np.concatenate(([0], b.astype(int), [0])))
    return (np.where(d == -1)[0] - np.where(d == 1)[0]) * NS


def main():
    f = OUT / 'data' / 'raw.npy'
    if not f.exists():
        print('no raw.npy yet'); return
    store = np.load(f, allow_pickle=True)[0]

    rows = []
    for key, reps in store.items():
        s = key.split('|')[0]
        chains = reps[0][1]['chains']
        for p in reps[0][1]['mind'].keys():
            pair = f'{chains[p[0]]}-{chains[p[1]]}'
            md = [r[1]['mind'][p] for r in reps]
            fb = np.array([(np.asarray(m) < BOUND_NM).mean() for m in md])
            ep = np.concatenate([episodes(m) for m in md]) if md else np.array([])
            ns_tot = sum(len(m) for m in md) * NS
            rows.append(dict(
                set=s, pair=pair, n_rep=len(md), ns_sampled=round(ns_tot),
                bound_frac=round(float(fb.mean()), 4),
                bound_frac_ci95=round(float(1.96 * fb.std(ddof=1) / np.sqrt(len(fb))), 4) if len(fb) > 1 else np.nan,
                n_episodes=int(ep.size),
                tau_mean_ns=round(float(ep.mean()), 2) if ep.size else 0.0,
                tau_median_ns=round(float(np.median(ep)), 2) if ep.size else 0.0,
                tau_max_ns=round(float(ep.max()), 1) if ep.size else 0.0,
                ep_over_10ns=int((ep > 10).sum()) if ep.size else 0))

    if not rows:
        print('no pairs'); return
    df = pd.DataFrame(rows).sort_values(['bound_frac', 'tau_mean_ns'], ascending=False)

    def verdict(r):
        if r.tau_mean_ns >= TAU_REAL and r.bound_frac >= FRAC_REAL:
            return 'REAL COMPLEX -> extend'
        if r.tau_mean_ns >= TAU_REAL or r.bound_frac >= FRAC_REAL:
            return 'borderline -> extend'
        return 'transient encounter -> done'

    df['verdict'] = df.apply(verdict, axis=1)
    df.to_csv(OUT / 'episode_stats.csv', index=False)
    print(df.to_string(index=False))

    ext = sorted(set(df.loc[df.verdict.str.contains('extend'), 'set']))
    (OUT / 'extend_sets.txt').write_text('\n'.join(ext) + ('\n' if ext else ''))
    print(f'\nsets flagged for deeper sampling: {len(ext)}'
          + (f'  -> {" ".join(ext)}' if ext else '  (none: every interface is a transient encounter)'))
    print(f'written: {OUT/"episode_stats.csv"}, {OUT/"extend_sets.txt"}')


if __name__ == '__main__':
    main()
