#!/usr/bin/env python
"""
Is the inter-chain contact network LIQUID or FROZEN?

This is the control that decides whether an inter-chain domain contact map can
be read as an interaction preference at all.

A liquid condensate continuously re-partners: the specific bead pairs in
contact at time t are mostly gone a few ns later, even though the total number
of contacts is steady. A kinetically trapped or interpenetrated arrangement
keeps the SAME pairs indefinitely.

The slab protocol compacts every system into a ~6-7 nm slab before production.
For a construct that does not condense, the chains then separate -- but
extended multi-domain chains that were threaded through each other at the end
of equilibration can stay interpenetrated, because the restrained domains are
rigid and cannot pass through one another. Those leftover contacts are present
in the very first production frame, never exchange, and dominate the raw
contact map. Reading them as "domain X likes domain Y" would be wrong.

METRIC
    retention(dt) = |pairs(t) & pairs(t+dt)| / |pairs(t)|

    ~0      fully exchanging  -> liquid, contact map is meaningful
    ~1      same pairs forever -> frozen, contact map is an initial condition

Measured at dt = 10 frames (~10 ns) and dt = end of run.

USAGE
    python analyze_contact_exchange.py --all --workers 12
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
PANEL = os.path.join(HERE, 'rnapanel')


def exchange(run_dir, cutoff=1.0, dt=10):
    import mdtraj as md
    from scipy.spatial import cKDTree

    run_dir = run_dir.rstrip('/')
    meta = yaml.safe_load(open(os.path.join(run_dir, 'slab_meta.yaml')))
    nres, nchain, sysname = meta['nres'], meta['nchain'], meta['sysname']
    dcd = os.path.join(run_dir, f'{sysname}.dcd')
    top = os.path.join(run_dir, 'checkpoint.pdb')
    if not (os.path.isfile(dcd) and os.path.isfile(top)):
        return None
    with md.open(dcd) as fh:
        n = len(fh)
    P = np.arange(nchain * nres)
    chain = P // nres

    def pairs(fi):
        fr = md.load_frame(dcd, int(fi), top=top)
        box = fr.unitcell_lengths[0].astype(float)
        xyz = np.mod(fr.xyz[0].astype(float), box)
        t = cKDTree(xyz[P], boxsize=box)
        pr = t.query_pairs(cutoff, output_type='ndarray')
        if not len(pr):
            return set()
        a, b = pr[:, 0], pr[:, 1]
        k = chain[a] != chain[b]
        return set(map(tuple, np.sort(np.c_[a[k], b[k]], axis=1)))

    mid = n // 2
    s1 = pairs(mid)
    if not s1:
        return None
    s2, s3 = pairs(min(mid + dt, n - 1)), pairs(n - 1)
    return dict(run=os.path.basename(run_dir), construct=meta['construct'],
                arm=meta['arm'], seed=meta['seed'],
                n_interchain_pairs=len(s1),
                retention_10ns=len(s1 & s2) / len(s1),
                retention_end=len(s1 & s3) / len(s1))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run_dir', nargs='?')
    ap.add_argument('--all', action='store_true')
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--out', default=os.path.join(HERE, 'analysis'))
    a = ap.parse_args()

    runs = ([d for d in sorted(glob.glob(os.path.join(PANEL, '*')))
             if os.path.isfile(os.path.join(d, 'slab_meta.yaml'))
             and not os.path.isfile(os.path.join(d, '.parked'))]
            if a.all else [os.path.abspath(a.run_dir)])

    from concurrent.futures import ProcessPoolExecutor, as_completed
    rows = []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(exchange, r): r for r in runs}
        for f in as_completed(futs):
            try:
                r = f.result()
            except Exception as e:
                print(f'  FAILED {os.path.basename(futs[f])}: {e}', flush=True)
                continue
            if r:
                rows.append(r)
                print(f"  {r['run']:<42} retain10ns={r['retention_10ns']:.3f}", flush=True)
    if rows:
        os.makedirs(a.out, exist_ok=True)
        p = os.path.join(a.out, 'contact_exchange.csv')
        pd.DataFrame(rows).sort_values(['construct', 'arm', 'seed']).to_csv(p, index=False)
        print(f'\nwrote {p}')


if __name__ == '__main__':
    main()
