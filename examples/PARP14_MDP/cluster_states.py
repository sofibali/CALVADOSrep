#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cluster CALVADOS simulation frames into top conformational states.

Uses inter-domain COM-COM distances (among MD1L1, MD2, MD3, ART by default)
as features. Standardizes, optionally PCA-reduces, then K-means clusters
into K (default 5) states. Extracts representative centroid frames per
state for downstream PyMOL visualization.

Outputs:
    data/cluster_features_{set}.npz             (raw features + metadata)
    data/cluster_assignments_{set}.npz          (cluster labels per frame)
    figures/cluster_silhouette_{set}.png        (k-sweep validation)
    figures/cluster_pca_{set}.png               (frames in 2D PC space)
    figures/cluster_distance_profile_{set}.png  (mean features per state)
    figures/cluster_population_{set}.png        (cluster sizes)
    figures/cluster_replicate_{set}.png         (which replicates → which state)
    representative_frames/{set}/state_{k}.pdb   (centroid frames as PDB)

Usage:
    python cluster_states.py --set fl_optimized
    python cluster_states.py --set fl_optimized --k 5 --pca 3
    python cluster_states.py --set fl_optimized --domains MD1L1 MD2 MD3 ART
    python cluster_states.py --set fl_optimized --add-rg     # include domain Rg
"""

import os
import json
import argparse
import numpy as np
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ============================================================
# Configuration
# ============================================================

CWD = Path(__file__).resolve().parent
DATA_PATH = CWD / 'data'
FIG_PATH = CWD / 'figures'
REP_PATH = CWD / 'representative_frames'
for p in (DATA_PATH, FIG_PATH, REP_PATH):
    p.mkdir(exist_ok=True)

DOMAINS_FL = {
    'RRM1':    (6, 88),
    'RRM2':    (146, 224), 'RRM3': (225, 314),
    'KH1-KH6': (315, 737), 'KH7a': (738, 789),
    'MD1L1':   (790, 1004), 'MD2': (1005, 1193), 'MD3': (1207, 1388),
    'KHb-KH8': (1389, 1533), 'WWE': (1534, 1602), 'ART': (1603, 1801),
}

SEEDS = range(1, 6)
SAMPLES = range(0, 5)
SKIP_FRAMES = 50


# ============================================================
# Feature extraction
# ============================================================

def extract_features_one_replicate(args):
    """Compute inter-domain features for one replicate.

    Returns dict with:
      'features':   (n_frames, n_features) array
      'metadata':   per-frame (seed, sample, frame_idx) tuples
      'feature_names': list of feature labels
    """
    set_key, seed, sample, domains, add_rg = args
    import MDAnalysis as mda

    if set_key.startswith('frag_'):
        sim_dir = CWD / 'fragments' / set_key[5:] / f'seed-{seed}_sample-{sample}'
    else:
        sim_dir = CWD / set_key / f'seed-{seed}_sample-{sample}'

    pdb = sim_dir / 'top.pdb'
    dcd = sim_dir / 'parp14.dcd'
    if not pdb.is_file() or not dcd.is_file():
        return None

    u = mda.Universe(str(pdb), str(dcd))
    ags = {}
    for dn in domains:
        ds, de = DOMAINS_FL[dn]
        ag = u.select_atoms(f'resid {ds}:{de}')
        if len(ag) > 0:
            ags[dn] = ag
    if not ags:
        return None

    # Build pair list
    domain_list = [d for d in domains if d in ags]
    pairs = [(domain_list[i], domain_list[j])
             for i in range(len(domain_list))
             for j in range(i + 1, len(domain_list))]
    feature_names = [f'dist_{a}_{b}' for a, b in pairs]
    if add_rg:
        feature_names += [f'rg_{d}' for d in domain_list]

    rows = []
    metadata = []
    for frame_idx, ts in enumerate(u.trajectory[SKIP_FRAMES:]):
        row = []
        for a, b in pairs:
            com_a = ags[a].center_of_geometry() / 10.0
            com_b = ags[b].center_of_geometry() / 10.0
            row.append(float(np.linalg.norm(com_a - com_b)))
        if add_rg:
            for d in domain_list:
                positions = ags[d].positions / 10.0
                com = positions.mean(axis=0)
                rg = float(np.sqrt(np.mean(np.sum((positions - com) ** 2, axis=1))))
                row.append(rg)
        rows.append(row)
        metadata.append((seed, sample, frame_idx + SKIP_FRAMES))

    return {
        'set_key': set_key,
        'features': np.array(rows),
        'metadata': metadata,
        'feature_names': feature_names,
    }


def collect_features(set_key, domains, add_rg, workers):
    """Run extraction across all replicates and concatenate."""
    jobs = []
    for seed in SEEDS:
        for sample in SAMPLES:
            jobs.append((set_key, seed, sample, tuple(domains), add_rg))

    print(f"  Processing {len(jobs)} replicates...")
    all_feats = []
    all_meta = []
    feature_names = None

    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(extract_features_one_replicate, j): j
                   for j in jobs}
        for fut in as_completed(futures):
            res = fut.result()
            if res is None:
                continue
            all_feats.append(res['features'])
            all_meta.extend(res['metadata'])
            feature_names = res['feature_names']
            print(f"    seed-{res['metadata'][0][0]}_sample-"
                  f"{res['metadata'][0][1]}: {len(res['features'])} frames")

    if not all_feats:
        return None
    return {
        'features': np.concatenate(all_feats, axis=0),
        'metadata': all_meta,
        'feature_names': feature_names,
    }


# ============================================================
# Clustering
# ============================================================

def silhouette_sweep(X, k_range, sample_size=5000):
    """Compute silhouette score for each k in k_range."""
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score

    # Subsample for speed
    rng = np.random.default_rng(42)
    if len(X) > sample_size:
        idx = rng.choice(len(X), sample_size, replace=False)
        X_sub = X[idx]
    else:
        X_sub = X

    scores = []
    inertias = []
    for k in k_range:
        km = KMeans(n_clusters=k, n_init=10, random_state=42)
        labels = km.fit_predict(X_sub)
        score = silhouette_score(X_sub, labels) if k > 1 else 0.0
        scores.append(score)
        inertias.append(km.inertia_)
        print(f"    k={k}: silhouette={score:.4f}, inertia={km.inertia_:.1f}")
    return np.array(scores), np.array(inertias)


def cluster_kmeans(X, k, seed=42):
    """K-means clustering. Returns (labels, centroids)."""
    from sklearn.cluster import KMeans
    km = KMeans(n_clusters=k, n_init=20, random_state=seed)
    labels = km.fit_predict(X)
    return labels, km.cluster_centers_, km


def find_centroid_frames(X, labels, k):
    """For each cluster, find the frame closest to the cluster centroid.

    Returns dict {cluster_id: frame_index_in_X}.
    """
    centroids = []
    for c in range(k):
        mask = labels == c
        if mask.sum() == 0:
            centroids.append(None)
            continue
        cluster_pts = X[mask]
        cluster_mean = cluster_pts.mean(axis=0)
        # Find the frame closest to the mean
        dists = np.linalg.norm(cluster_pts - cluster_mean, axis=1)
        local_idx = np.argmin(dists)
        # Map back to global index
        global_idx = np.where(mask)[0][local_idx]
        centroids.append(int(global_idx))
    return centroids


# ============================================================
# Extract representative frame as PDB
# ============================================================

def extract_pdb_for_frame(set_key, seed, sample, frame_idx, out_pdb):
    """Extract a single frame as PDB from a replicate trajectory."""
    import MDAnalysis as mda

    if set_key.startswith('frag_'):
        sim_dir = CWD / 'fragments' / set_key[5:] / f'seed-{seed}_sample-{sample}'
    else:
        sim_dir = CWD / set_key / f'seed-{seed}_sample-{sample}'

    pdb = sim_dir / 'top.pdb'
    dcd = sim_dir / 'parp14.dcd'
    u = mda.Universe(str(pdb), str(dcd))
    if frame_idx >= len(u.trajectory):
        return False
    u.trajectory[frame_idx]
    with mda.Writer(str(out_pdb), n_atoms=u.atoms.n_atoms) as W:
        W.write(u.atoms)
    return True


# ============================================================
# Plots
# ============================================================

def plot_silhouette(k_range, sil, inertias, outpath):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
    ax1.plot(k_range, sil, 'o-', color='#9b59b6', linewidth=2,
             markersize=8)
    ax1.set_xlabel('Number of clusters k', fontsize=11)
    ax1.set_ylabel('Silhouette score', fontsize=11)
    ax1.set_title('Cluster validity (higher = better)', fontsize=12)
    ax1.set_xticks(list(k_range))
    ax1.grid(alpha=0.3)

    ax2.plot(k_range, inertias, 's-', color='#e67e22', linewidth=2,
             markersize=8)
    ax2.set_xlabel('Number of clusters k', fontsize=11)
    ax2.set_ylabel('Inertia (within-cluster SSE)', fontsize=11)
    ax2.set_title('Elbow plot (look for kink)', fontsize=12)
    ax2.set_xticks(list(k_range))
    ax2.grid(alpha=0.3)

    plt.tight_layout()
    fig.savefig(f'{outpath}.png', dpi=200, bbox_inches='tight')
    plt.close()


def plot_pca_clusters(X_pca, labels, k, centroid_indices, outpath):
    """2D PCA scatter colored by cluster."""
    fig, ax = plt.subplots(figsize=(10, 8))
    colors = plt.cm.tab10(np.linspace(0, 1, k))

    # Plot each cluster
    for c in range(k):
        mask = labels == c
        ax.scatter(X_pca[mask, 0], X_pca[mask, 1], s=2, alpha=0.3,
                   c=[colors[c]], label=f'State {c+1} (n={mask.sum()})',
                   edgecolors='none')

    # Mark centroid frames
    for c, idx in enumerate(centroid_indices):
        if idx is not None:
            ax.plot(X_pca[idx, 0], X_pca[idx, 1], marker='*',
                    markersize=20, color=colors[c],
                    markeredgecolor='black', markeredgewidth=1.5, zorder=10)

    ax.set_xlabel('PC1', fontsize=12)
    ax.set_ylabel('PC2', fontsize=12)
    ax.set_title('Clustered conformational states in 2D PCA\n'
                 '(stars = centroid frames per state)', fontsize=13)
    ax.legend(fontsize=10, loc='best')

    plt.tight_layout()
    fig.savefig(f'{outpath}.png', dpi=200, bbox_inches='tight')
    plt.close()


def plot_state_distance_profile(X_raw, labels, feature_names, k, outpath):
    """Mean feature values per state (rows = states, cols = features).

    Annotated heatmap.
    """
    means = np.zeros((k, X_raw.shape[1]))
    stds = np.zeros((k, X_raw.shape[1]))
    for c in range(k):
        mask = labels == c
        if mask.sum() > 0:
            means[c] = X_raw[mask].mean(axis=0)
            stds[c] = X_raw[mask].std(axis=0)

    fig, ax = plt.subplots(figsize=(max(8, len(feature_names) * 0.8),
                                     max(4, k * 0.6)))
    im = ax.imshow(means, cmap='RdYlBu_r', aspect='auto')
    ax.set_xticks(range(len(feature_names)))
    ax.set_xticklabels(feature_names, rotation=45, ha='right', fontsize=9)
    ax.set_yticks(range(k))
    ax.set_yticklabels([f'State {c+1}' for c in range(k)], fontsize=10)
    plt.colorbar(im, ax=ax, label='Mean distance (nm)', shrink=0.8)

    # Annotate
    for i in range(k):
        for j in range(len(feature_names)):
            val = means[i, j]
            color = 'white' if val < means.min() + 0.3 * (means.max() - means.min()) else 'black'
            ax.text(j, i, f'{val:.2f}\n±{stds[i,j]:.2f}',
                    ha='center', va='center', fontsize=7, color=color)

    ax.set_title('Mean feature values per cluster state', fontsize=13)
    plt.tight_layout()
    fig.savefig(f'{outpath}.png', dpi=200, bbox_inches='tight')
    plt.close()


def plot_population(labels, k, outpath):
    """Bar chart of cluster populations + percentages."""
    counts = np.bincount(labels, minlength=k)
    fracs = counts / counts.sum()

    fig, ax = plt.subplots(figsize=(8, 5))
    colors = plt.cm.tab10(np.linspace(0, 1, k))
    bars = ax.bar(range(1, k + 1), counts, color=colors, edgecolor='black',
                  linewidth=0.5)
    for bar, frac in zip(bars, fracs):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                f'{frac*100:.1f}%', ha='center', va='bottom', fontsize=10)
    ax.set_xlabel('State', fontsize=12)
    ax.set_ylabel('Number of frames', fontsize=12)
    ax.set_title(f'Cluster populations (total {counts.sum()} frames)',
                 fontsize=13)
    ax.set_xticks(range(1, k + 1))
    plt.tight_layout()
    fig.savefig(f'{outpath}.png', dpi=200, bbox_inches='tight')
    plt.close()


def plot_replicate_state_distribution(metadata, labels, k, outpath):
    """Heatmap showing how each replicate's frames distribute across states."""
    rep_state = {}  # (seed, sample) -> [count per state]
    for (seed, sample, _), c in zip(metadata, labels):
        key = (seed, sample)
        if key not in rep_state:
            rep_state[key] = np.zeros(k, dtype=int)
        rep_state[key][c] += 1

    rep_keys = sorted(rep_state.keys())
    mat = np.array([rep_state[k] for k in rep_keys])
    mat_norm = mat / mat.sum(axis=1, keepdims=True)

    fig, ax = plt.subplots(figsize=(8, max(6, len(rep_keys) * 0.25)))
    im = ax.imshow(mat_norm, cmap='YlGnBu', aspect='auto', vmin=0, vmax=1)
    ax.set_xticks(range(k))
    ax.set_xticklabels([f'State {c+1}' for c in range(k)], fontsize=10)
    ax.set_yticks(range(len(rep_keys)))
    ax.set_yticklabels([f'seed-{s}_sample-{p}' for s, p in rep_keys],
                       fontsize=7)
    plt.colorbar(im, ax=ax, label='Fraction of frames', shrink=0.8)

    # Annotate
    for i in range(len(rep_keys)):
        for j in range(k):
            v = mat_norm[i, j]
            if v > 0.05:
                ax.text(j, i, f'{v*100:.0f}', ha='center', va='center',
                        fontsize=6,
                        color='white' if v > 0.5 else 'black')

    ax.set_title('Replicate-by-state population (% of frames)', fontsize=12)
    plt.tight_layout()
    fig.savefig(f'{outpath}.png', dpi=200, bbox_inches='tight')
    plt.close()


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--set', required=True,
                        help='Set name (e.g. fl_optimized, fl, or fragment '
                             'name like md1l1_md2_md3)')
    parser.add_argument('--k', type=int, default=5,
                        help='Number of clusters (default 5)')
    parser.add_argument('--domains', nargs='+',
                        default=['MD1L1', 'MD2', 'MD3', 'ART'],
                        help='Domains to use for inter-domain distances')
    parser.add_argument('--add-rg', action='store_true',
                        help='Include per-domain Rg as additional features')
    parser.add_argument('--pca', type=int, default=None,
                        help='Reduce to N PCA components before clustering '
                             '(default: keep 95%% variance)')
    parser.add_argument('--no-silhouette', action='store_true',
                        help='Skip silhouette k-sweep (faster)')
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--skip-pdb', action='store_true',
                        help='Skip extracting representative PDBs')
    args = parser.parse_args()

    set_key = args.set
    # Auto-prefix fragments
    if set_key not in ('fl', 'fl_optimized', 'md', 'core', 'mka',
                       'norrm', 'noart', 'md3art'):
        if (CWD / 'fragments' / set_key).is_dir():
            set_key = f'frag_{set_key}'

    print(f"=" * 60)
    print(f"Clustering {set_key}")
    print(f"=" * 60)
    print(f"  Domains:     {args.domains}")
    print(f"  Add Rg:      {args.add_rg}")
    print(f"  K:           {args.k}")

    # ── Step 1: extract features ──
    print("\n[1] Extracting features per frame...")
    cache_npz = DATA_PATH / f'cluster_features_{set_key.replace("frag_","")}.npz'

    if cache_npz.exists():
        print(f"  Loading cached features from {cache_npz}")
        cached = np.load(cache_npz, allow_pickle=True)
        X_raw = cached['features']
        metadata = list(cached['metadata'])
        feature_names = list(cached['feature_names'])
    else:
        data = collect_features(set_key, args.domains, args.add_rg,
                                args.workers)
        if data is None:
            print("ERROR: no replicates found")
            return 1
        X_raw = data['features']
        metadata = data['metadata']
        feature_names = data['feature_names']
        np.savez(cache_npz,
                 features=X_raw,
                 metadata=np.array(metadata, dtype=object),
                 feature_names=np.array(feature_names))
        print(f"  Saved cache: {cache_npz}")

    print(f"  Feature matrix: {X_raw.shape}")
    print(f"  Features: {feature_names}")

    # ── Step 2: standardize ──
    from sklearn.preprocessing import StandardScaler
    scaler = StandardScaler()
    X_std = scaler.fit_transform(X_raw)

    # ── Step 3: PCA ──
    from sklearn.decomposition import PCA
    if args.pca:
        n_pca = min(args.pca, X_std.shape[1])
    else:
        # Keep 95% variance
        pca_full = PCA().fit(X_std)
        n_pca = int(np.argmax(np.cumsum(pca_full.explained_variance_ratio_)
                              >= 0.95) + 1)
        n_pca = max(2, n_pca)
    pca = PCA(n_components=n_pca, random_state=42)
    X_pca = pca.fit_transform(X_std)
    print(f"\n[2] PCA: {n_pca} components, "
          f"explained variance = "
          f"{pca.explained_variance_ratio_.sum()*100:.1f}%")
    print(f"  Components: "
          f"{[f'{r*100:.1f}%' for r in pca.explained_variance_ratio_]}")

    # ── Step 4: silhouette sweep ──
    suffix = set_key.replace('frag_', '')
    if not args.no_silhouette:
        print("\n[3] Silhouette k-sweep (k=2..10)...")
        k_range = range(2, 11)
        sil, inertias = silhouette_sweep(X_pca, k_range)
        best_k = list(k_range)[int(np.argmax(sil))]
        print(f"  Best k by silhouette: {best_k}")
        plot_silhouette(k_range, sil, inertias,
                         FIG_PATH / f'cluster_silhouette_{suffix}')

    # ── Step 5: K-means with chosen k ──
    print(f"\n[4] K-means with k={args.k}...")
    labels, centroids, km = cluster_kmeans(X_pca, args.k)
    print(f"  Cluster sizes: {np.bincount(labels)}")

    # Find centroid frames (closest to mean in PCA space)
    centroid_indices = find_centroid_frames(X_pca, labels, args.k)

    # Save assignments
    out_assign = DATA_PATH / f'cluster_assignments_{suffix}.npz'
    np.savez(out_assign,
             labels=labels,
             centroid_frames=np.array([
                 metadata[idx] if idx is not None else (-1, -1, -1)
                 for idx in centroid_indices]),
             feature_names=np.array(feature_names),
             X_raw=X_raw)
    print(f"  Saved: {out_assign}")

    # ── Step 6: plots ──
    print("\n[5] Generating plots...")
    plot_pca_clusters(X_pca, labels, args.k, centroid_indices,
                       FIG_PATH / f'cluster_pca_{suffix}')
    plot_state_distance_profile(X_raw, labels, feature_names, args.k,
                                 FIG_PATH / f'cluster_distance_profile_{suffix}')
    plot_population(labels, args.k,
                     FIG_PATH / f'cluster_population_{suffix}')
    plot_replicate_state_distribution(metadata, labels, args.k,
                                       FIG_PATH / f'cluster_replicate_{suffix}')

    # ── Step 7: extract representative PDBs ──
    if not args.skip_pdb:
        print("\n[6] Extracting representative frames...")
        rep_dir = REP_PATH / suffix
        rep_dir.mkdir(exist_ok=True)
        rep_info = []
        for c, idx in enumerate(centroid_indices):
            if idx is None:
                continue
            seed, sample, frame_idx = metadata[idx]
            pop = (labels == c).sum()
            pct = pop / len(labels) * 100
            out_pdb = rep_dir / f'state_{c+1}.pdb'
            ok = extract_pdb_for_frame(set_key, seed, sample, frame_idx,
                                        out_pdb)
            status = 'OK' if ok else 'FAIL'
            print(f"  State {c+1}: seed-{seed}_sample-{sample} "
                  f"frame {frame_idx} → {out_pdb.name} ({status}) "
                  f"[{pop} frames, {pct:.1f}%]")
            rep_info.append({
                'state': c + 1,
                'population': int(pop),
                'fraction': float(pct / 100),
                'representative_seed': int(seed),
                'representative_sample': int(sample),
                'representative_frame': int(frame_idx),
                'representative_pdb': str(out_pdb.relative_to(CWD)),
                'mean_features': {fn: float(X_raw[labels == c][:, j].mean())
                                  for j, fn in enumerate(feature_names)},
            })

        # Save state info JSON
        out_info = rep_dir / 'state_info.json'
        with open(out_info, 'w') as f:
            json.dump({
                'set': set_key,
                'k': args.k,
                'feature_names': feature_names,
                'states': rep_info,
            }, f, indent=2)
        print(f"  Saved: {out_info}")

    # ── Summary ──
    print(f"\n{'=' * 60}")
    print(f"State summary ({args.k} clusters, {len(labels)} frames total)")
    print(f"{'=' * 60}")
    print(f"  {'State':<7} {'N':>6} {'%':>6}  Mean features (nm)")
    for c in range(args.k):
        mask = labels == c
        n = mask.sum()
        pct = n / len(labels) * 100
        means = X_raw[mask].mean(axis=0)
        means_str = ' '.join(f'{fn[:8]}:{v:.2f}'
                              for fn, v in zip(feature_names, means))
        print(f"  {c+1:<7} {n:>6} {pct:>5.1f}%  {means_str}")


if __name__ == '__main__':
    main()
