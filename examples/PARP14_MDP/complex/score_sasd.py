#!/usr/bin/env python3
"""Re-score the published crosslinks with solvent-accessible surface distance.

Every number in this campaign so far uses the Euclidean Ca-Ca distance, which
lets the crosslinker pass straight through the protein. The field standard is
SASD: the shortest path that stays in solvent, which is what decides whether a
crosslink can actually form. sim_analysis/analyze_lys_contacts.compute_sasd
already implements it for this bead model and is reused here unchanged, so the
two analyses cannot drift apart.

Applied to the 36 published inter-protein crosslinks only. The 4,716 decoy
pairs are left on Euclidean distance: SASD costs ~0.26 s per pair-frame, so the
full decoy set would be weeks of CPU. The decoy AUC in score_decoys.py is the
dimension control; this is the correction to the headline numbers.

Both chains are placed in the same periodic image before any path search --
compute_sasd works in one image, and 90% of frames here need the shift.

Expect SASD to cost CALVADOS less than HyRes: a compact, extensively
self-contacted chain buries more of its lysines, and a buried lysine cannot be
crosslinked at any distance.

Usage:
    python score_sasd.py --all --target-frames 300 --workers 48
"""
import sys, warnings
from pathlib import Path
from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor
import numpy as np, pandas as pd
warnings.filterwarnings('ignore')

HERE = Path(__file__).resolve().parent
OUT = HERE / 'analysis'
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / 'sim_analysis'))
from analyze_lys_contacts import compute_sasd
from pbc_center import center_pair
import score_decoys as sd

THRESH_NM = 3.0                      # 30 A, the project's LYS_LYS_CUTOFF


def _frame_task(args):
    """SASD for every requested pair in one frame. Returns [(k, euclid, sasd)]."""
    pos_A, n_first, box_nm, idx_pairs = args
    # center_pair takes Angstrom; compute_sasd takes nm
    pos = center_pair(pos_A, n_first, box_nm)[0] / 10.0   # returns (pos, ref, best)
    out = []
    for k, (i, j) in enumerate(idx_pairs):
        e = float(np.linalg.norm(pos[i] - pos[j]))
        # only pairs that already pass on Euclidean can be rescued or lost;
        # SASD is never shorter than Euclidean, so a failing pair stays failing
        s = compute_sasd(pos, i, j) if e < THRESH_NM else np.inf
        out.append((k, e, s))
    return out


def run_model(label, spec, target_frames, workers):
    allp, rep = sd.published('parp9', 'dtx3l')
    pairs = sorted(allp)
    off_p9, off_dx = 0, sd.ac.LENGTHS['parp9']
    idx_pairs = [(off_p9 + a - 1, off_dx + b - 1) for a, b in pairs]

    if 'hyres' in spec:
        it = sd.frames_hyres(spec['hyres'], target_frames)
        box_nm = sd.HYRES_BOX_A / 10.0
    else:
        it = sd.frames_calvados(spec['s'], spec['root'], target_frames)
        box_nm = None

    tasks, n_first = [], sd.ac.LENGTHS['parp9']
    for P, boxv, off in it:
        b = box_nm if box_nm is not None else float(boxv[0]) / 10.0
        tasks.append((P.copy(), n_first, b, idx_pairs))
    if not tasks:
        print(f'   {label}: no frames'); return None

    hit_e = np.zeros(len(pairs), int); hit_s = np.zeros(len(pairs), int)
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for res in ex.map(_frame_task, tasks, chunksize=4):
            for k, e, s in res:
                hit_e[k] += e < THRESH_NM
                hit_s[k] += s < THRESH_NM
    n = len(tasks)
    rows = [dict(model=label, res_a=a, res_b=b, repro=(a, b) in rep,
                 n_frames=n,
                 pct_euclid=round(100.0 * hit_e[k] / n, 3),
                 pct_sasd=round(100.0 * hit_s[k] / n, 3))
            for k, (a, b) in enumerate(pairs)]
    df = pd.DataFrame(rows)
    r = df[df.repro]
    print(f'   {label:34s} {n:4d} frames | ALL euclid {df.pct_euclid.mean():6.2f}% '
          f'-> SASD {df.pct_sasd.mean():6.2f}%  | repro {r.pct_euclid.mean():6.2f}% '
          f'-> {r.pct_sasd.mean():6.2f}%')
    return df


def main():
    ap = ArgumentParser()
    ap.add_argument('--all', action='store_true')
    ap.add_argument('--target-frames', type=int, default=300)
    ap.add_argument('--workers', type=int, default=48)
    a = ap.parse_args()
    jobs = sd.JOBS if a.all else sd.JOBS[:1]
    print(f'SASD re-scoring, {THRESH_NM*10:.0f} A, {a.workers} workers')
    got = []
    for label, spec in jobs:
        if 'hyres' not in spec:
            base = spec['s'] if spec['s'] in sd.SETS else 'p9_dtx3l'
            if not {'parp9', 'dtx3l'} <= set(sd.SETS[base]):
                continue
        d = run_model(label, spec, a.target_frames, a.workers)
        if d is not None:
            got.append(d)
    if got:
        df = pd.concat(got, ignore_index=True)
        df.to_csv(OUT / 'sasd_rescore.csv', index=False)
        print('\nwrote', OUT / 'sasd_rescore.csv')


if __name__ == '__main__':
    main()
