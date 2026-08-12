#!/usr/bin/env python
"""
STAGE 2 - Select k representative states from a TICA projection (stage 1 output)
and write one CG PDB per state.

--rep-mode:
  spread  (default) farthest-point sampling across the populated 2D landscape +
          Voronoi reassignment -> states are EVENLY DISTRIBUTED and non-overlapping.
          Best when the system is diffusive (implied timescales don't plateau): the
          states are "distinct sampled conformations", not Markovian basins.
  pcca    MSM (at --msm-lag) + PCCA+ -> kinetic metastable macrostates. Use only when
          the implied timescales converge (a real spectral gap exists).

Example:
  python 02_select_states.py --tica-npz ../data/tica_md_full_pose.npz \
      --k 4 --rep-mode spread --out-dir ../states/md_full_pose
"""
import os, sys, json, argparse
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from cluster_states import (cluster_kmeans, build_msm_transition_matrix, pcca_plus,  # noqa: E402
                            find_spread_frames, find_centroid_frames, extract_pdb_for_frame)


def main():
    ap = argparse.ArgumentParser(description="Stage 2: TICA projection -> k representative states")
    ap.add_argument('--tica-npz', required=True, help='Stage-1 .npz')
    ap.add_argument('--k', type=int, default=4, help='Number of representative states')
    ap.add_argument('--rep-mode', choices=['spread', 'pcca'], default='spread')
    ap.add_argument('--n-microstates', type=int, default=200,
                    help='k-means microstates before macrostate selection')
    ap.add_argument('--msm-lag', type=int, default=200, help='MSM lag in frames (pcca mode)')
    ap.add_argument('--out-dir', required=True, help='Where to write state_*.pdb')
    a = ap.parse_args()

    d = np.load(a.tica_npz, allow_pickle=True)
    Y = d['Y']
    meta = [tuple(int(x) for x in m) for m in d['metadata']]
    rb = [tuple(int(x) for x in b) for b in d['rep_boundaries']]
    set_key = str(d['set_key'])
    print(f"loaded {a.tica_npz}: set='{set_key}' feature='{d['features']}' "
          f"frames={len(Y)} TICA-dims={Y.shape[1]}")

    print(f"[1] {a.n_microstates} microstates (k-means on TICA) ...")
    micro, _, _ = cluster_kmeans(Y, a.n_microstates)

    if a.rep_mode == 'pcca':
        print(f"[2] MSM (lag {a.msm_lag}) + PCCA+ -> {a.k} metastable macrostates ...")
        T, pi, C = build_msm_transition_matrix(micro, rb, a.msm_lag, a.n_microstates)
        m2M, _ = pcca_plus(T, a.k)
        labels = np.array([m2M[m] for m in micro])
        idx = find_centroid_frames(Y, labels, a.k)
    else:
        print(f"[2] SPREAD: farthest-point sampling -> {a.k} states (Voronoi tiling) ...")
        idx, labels = find_spread_frames(Y[:, :2], a.k)

    os.makedirs(a.out_dir, exist_ok=True)
    info = []
    print(f"{'state':6} {'rep frame':>26} {'pop%':>6}  status")
    for c, fi in enumerate(idx):
        if fi is None:
            continue
        seed, sample, frame = meta[fi]
        pop = int((labels == c).sum()); pct = pop / len(labels) * 100
        out = os.path.join(a.out_dir, f'state_{c + 1}.pdb')
        ok = extract_pdb_for_frame(set_key, seed, sample, frame, out)
        print(f"{c + 1:6} {f'seed-{seed}_sample-{sample} fr{frame}':>26} "
              f"{pct:6.1f}  {'OK' if ok else 'FAIL'}")
        info.append(dict(state=c + 1, seed=seed, sample=sample, frame=frame,
                         population=pop, fraction=pop / len(labels)))
    json.dump(dict(set=set_key, features=str(d['features']), k=a.k,
                   rep_mode=a.rep_mode, states=info),
              open(os.path.join(a.out_dir, 'state_info.json'), 'w'), indent=2)
    print(f"    wrote {len(info)} states -> {a.out_dir}")


if __name__ == '__main__':
    main()
