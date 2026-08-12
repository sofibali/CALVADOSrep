#!/usr/bin/env python
"""
STAGE 1 - Featurize a CALVADOS simulation set and project it with TICA.

Reads all replicate trajectories of --set, computes the chosen featurization, runs
TICA at --tica-lag, and saves the TICA projection (+ metadata) for stage 2.

Featurizations (--features):
  pose         rigid-body relative pose per domain pair (position + orientation) -- best
               for "the distinct arrangements you see by eye". (rotation-aware)
  pose_iface   pose + interface contact distances (arrangement + contact register)
  ca           pairwise CA-CA distances, every --ca-stride-th residue (rotation-blind)
  interface_ca domain-edge CA-CA contact distances
  com          inter-domain COM-COM distances (coarsest)
  orient/linker_ca/torsion/segments/inter_exposed  (see cluster_states.py)

Lag is in FRAMES (1 frame = 0.01 ns here); 200 = 2 ns.

Example:
  python 01_featurize_tica.py --set md_full --features pose --tica-lag 200 \
      --out ../data/tica_md_full_pose.npz
"""
import os, sys, argparse
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from cluster_states import collect_features, compute_tica  # noqa: E402

ALL_DOMAINS = ['RRM1', 'RRM2', 'RRM3', 'KH1-KH6', 'KH7a',
               'MD1L1', 'MD2', 'MD3', 'KHb-KH8', 'WWE', 'ART']


def main():
    ap = argparse.ArgumentParser(description="Stage 1: features -> TICA projection")
    ap.add_argument('--set', required=True, dest='set_key',
                    help='Simulation set (e.g. md_full, fl_optimized, or a fragment name)')
    ap.add_argument('--features', default='pose',
                    choices=['com', 'ca', 'pose', 'pose_iface', 'interface_ca',
                             'orient', 'linker_ca', 'torsion', 'segments', 'inter_exposed'])
    ap.add_argument('--ca-stride', type=int, default=25,
                    help='Residue stride for ca / interface / pose_iface edge size (default 25)')
    ap.add_argument('--tica-lag', type=int, default=200,
                    help='TICA lag in FRAMES (1 frame = 0.01 ns). Default 200 = 2 ns')
    ap.add_argument('--n-components', type=int, default=10, help='TICA components to keep')
    ap.add_argument('--domains', nargs='+', default=ALL_DOMAINS,
                    help='Domains to featurize (auto-filtered to those present in the set)')
    ap.add_argument('--workers', type=int, default=16)
    ap.add_argument('--out', required=True, help='Output .npz (TICA projection + metadata)')
    a = ap.parse_args()

    print(f"[1] Featurizing set='{a.set_key}' feature='{a.features}' (stride {a.ca_stride}) ...")
    d = collect_features(a.set_key, a.domains, False, a.workers,
                         features_mode=a.features, ca_stride=a.ca_stride)
    X = d['features']
    print(f"    feature matrix: {X.shape}  ({len(d['rep_boundaries'])} replicates)")

    print(f"[2] TICA at lag {a.tica_lag} frames ({a.tica_lag * 0.01:.2f} ns), "
          f"{a.n_components} components ...")
    Y, eigvals, eigvecs = compute_tica(X, d['rep_boundaries'],
                                       lag=a.tica_lag, n_components=a.n_components)

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    np.savez(a.out, Y=Y, X_raw=X,
             metadata=np.array(d['metadata']),
             rep_boundaries=np.array(d['rep_boundaries']),
             feature_names=np.array(d['feature_names']),
             eigvals=eigvals, set_key=a.set_key,
             features=a.features, tica_lag=a.tica_lag)
    print(f"    saved -> {a.out}   (TICA projection {Y.shape})")


if __name__ == '__main__':
    main()
