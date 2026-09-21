#!/usr/bin/env python3
"""Representative bound states of each complex, by PCA or TICA + k-means.

WHAT IS CLUSTERED, AND WHY
--------------------------
Features are INTER-CHAIN domain-domain minimum distances -- which domain of one
chain is near which domain of the other. Two deliberate choices:

* Not intra-domain contacts. Every folded domain here is held by harmonic
  restraints, so its internal contact map is frozen by construction; the
  project's own featurization notes say as much ("mostly frozen by restraints;
  usually scores low"). The slow, informative coordinate in a binding run is the
  interface REGISTER between chains, not anything inside a domain.

* Bound frames only. These complexes are unbound 75-94% of the time. Clustering
  every frame would mostly partition "how far apart", and the representative
  structures would be pictures of two chains in different corners of the box.
  Restricting to frames with a sub-cutoff inter-chain contact makes the states
  mean "which interface", which is the question.

TICA vs PCA
-----------
PCA finds directions of maximum variance; TICA finds the SLOWEST ones, which is
what you want for metastable states. TICA needs time-lagged pairs, and the bound
frames are not contiguous -- so lagged pairs are formed only WITHIN a continuous
bound episode, never across a dissociation. Episodes shorter than the lag
contribute nothing, which for a ~1 ns-lifetime system means TICA leans on the
rare long episodes. With --reduce pca that problem does not arise. Both are
provided; PCA is the default because it is the honest one here.

No scikit-learn in the CALVADOS env, and installing into a shared env during a
live campaign is not worth it -- PCA, TICA, k-means++ and silhouette are
implemented directly below.

Usage:
    python cluster_complex_states.py                 # every set, PCA
    python cluster_complex_states.py ternary_docked --reduce tica --k 4
"""
import sys, warnings, itertools
from pathlib import Path
from argparse import ArgumentParser
import numpy as np, yaml
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array

# Validated project palette (see the style guide's PROTEIN/FACE/CONSTRUCT sets).
STATE_COLORS = ['#4A3AA7', '#C8481A', '#117C60', '#EDA100', '#A05195', '#2A78D6']
warnings.filterwarnings('ignore')

HERE = Path(__file__).resolve().parent
ROOT = HERE / 'binding'
OUT = HERE / 'analysis'
FIG = OUT / 'figures'
STATES = OUT / 'states'
for d in (FIG, STATES, OUT / 'data'):
    d.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(HERE))
from analyze_binding import SETS, sysname, topology
from figure_face_separation import load_domains, DOMAIN_NAMES
import analyze_convergence as ac

NS_PER_FRAME = 0.05
BOUND_NM = 1.0
EQ_FRAMES = 500          # 25 ns
DEFAULT_STRIDE = 2       # 0.1 ns


# ---------------------------------------------------------------- linear algebra
def pca(X, n=2):
    Xc = X - X.mean(0)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    var = S**2 / max(len(X) - 1, 1)
    return Xc @ Vt[:n].T, Vt[:n], var[:n] / var.sum()


def tica(X, episodes, lag, n=2):
    """Time-lagged independent components; lagged pairs only within an episode."""
    Xc = X - X.mean(0)
    pairs = [(a, b) for (s, e) in episodes for a, b in
             zip(range(s, e - lag), range(s + lag, e))]
    if len(pairs) < 10 * X.shape[1]:
        return None
    i = np.array([p[0] for p in pairs]); j = np.array([p[1] for p in pairs])
    C0 = (Xc.T @ Xc) / len(Xc)
    Ct = (Xc[i].T @ Xc[j] + Xc[j].T @ Xc[i]) / (2 * len(i))
    C0 += 1e-6 * np.eye(len(C0))
    w, V = np.linalg.eigh(C0)
    keep = w > 1e-10
    W = V[:, keep] / np.sqrt(w[keep])
    Cw = W.T @ Ct @ W
    ew, ev = np.linalg.eigh((Cw + Cw.T) / 2)
    order = np.argsort(ew)[::-1][:n]
    comps = (W @ ev[:, order]).T
    return Xc @ comps.T, comps, ew[order]


def kmeans(X, k, seed=0, iters=200):
    rng = np.random.default_rng(seed)
    c = [X[rng.integers(len(X))]]
    for _ in range(k - 1):                       # k-means++
        d = np.min([((X - ci)**2).sum(1) for ci in c], axis=0)
        c.append(X[rng.choice(len(X), p=d / d.sum())] if d.sum() > 0 else X[rng.integers(len(X))])
    C = np.array(c)
    for _ in range(iters):
        lab = np.argmin(((X[:, None, :] - C[None]) ** 2).sum(2), axis=1)
        newC = np.array([X[lab == i].mean(0) if (lab == i).any() else C[i] for i in range(k)])
        if np.allclose(newC, C):
            break
        C = newC
    return lab, C


def silhouette(X, lab, cap=2000, seed=0):
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(X), min(cap, len(X)), replace=False)
    Xs, ls = X[idx], lab[idx]
    D = np.sqrt(((Xs[:, None, :] - Xs[None]) ** 2).sum(2))
    out = []
    for i in range(len(Xs)):
        same = ls == ls[i]; same[i] = False
        if not same.any():
            continue
        a = D[i, same].mean()
        b = min(D[i, ls == o].mean() for o in set(ls) - {ls[i]}) if len(set(ls)) > 1 else np.inf
        out.append((b - a) / max(a, b))
    return float(np.mean(out)) if out else np.nan


# ---------------------------------------------------------------- featurization
def featurize(s, stride):
    chains = SETS[s]
    blocks = load_domains(s, chains)
    names = [DOMAIN_NAMES.get(c, []) for c in chains]
    sysn = sysname(s)
    feats, mind_all, epi, srcs = [], [], [], []
    labels = None
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
        resids = u.atoms.resids
        # per-chain domain residue index slices
        dom_idx = []
        for ci, c in enumerate(chains):
            a, b = offs[ci]
            r = resids[a:b]
            dom_idx.append([(np.where((r >= lo) & (r <= hi))[0] + a)
                            for lo, hi in blocks[ci]])
        cpairs = list(itertools.combinations(range(len(chains)), 2))
        if labels is None:
            labels = [f'{chains[i]}:{(names[i][x] if x < len(names[i]) else x)}'
                      f'-{chains[j]}:{(names[j][y] if y < len(names[j]) else y)}'
                      for (i, j) in cpairs
                      for x in range(len(dom_idx[i])) for y in range(len(dom_idx[j]))]
        start = len(feats)
        for k, ts in enumerate(u.trajectory[EQ_FRAMES::stride]):
            P = u.atoms.positions.astype(np.float32)
            # One distance_array per CHAIN pair, then slice out each domain
            # block's minimum. Doing it per domain pair instead means ~286
            # distance_array calls per frame for the ternary set; this is one
            # call per chain pair and numpy slicing for the rest.
            row, gmin = [], np.inf
            for (i, j) in cpairs:
                D = distance_array(P[offs[i][0]:offs[i][1]],
                                   P[offs[j][0]:offs[j][1]], box=boxv) / 10.0
                bi, bj = offs[i][0], offs[j][0]
                for xi in dom_idx[i]:
                    Dx = D[xi - bi]
                    for yj in dom_idx[j]:
                        d = float(Dx[:, yj - bj].min())
                        row.append(d); gmin = min(gmin, d)
            feats.append(row); mind_all.append(gmin)
            srcs.append((rep.name, EQ_FRAMES + k * stride))
        epi.append((start, len(feats)))
    return (np.asarray(feats, np.float32), np.asarray(mind_all, np.float32),
            labels, srcs, epi)


# ---------------------------------------------------------------- main per set
def bound_episodes(mask, offsets):
    """Contiguous runs of bound frames, not crossing a replicate boundary."""
    out = []
    for (s, e) in offsets:
        m = mask[s:e]
        d = np.diff(np.concatenate(([0], m.astype(int), [0])))
        for a, b in zip(np.where(d == 1)[0], np.where(d == -1)[0]):
            out.append((s + a, s + b))
    return out


def write_state_pdb(s, src, out_pdb):
    """Write the representative frame as a PDB."""
    rep_name, frame = src
    rep = ROOT / s / rep_name
    u = mda.Universe(str(topology(rep)), str(rep / f'{sysname(s)}.dcd'))
    u.trajectory[frame]
    u.atoms.write(str(out_pdb))


def run_set(s, reduce_mode, kmax, stride, lag_ns, fixed_k=None):
    X, mind, labels, srcs, offs = featurize(s, stride)
    if not len(X):
        print(f'  {s}: no frames'); return None
    mask = mind < BOUND_NM
    nb = int(mask.sum())
    print(f'  {s}: {len(X)} frames, {nb} bound ({100*nb/len(X):.1f}%)', flush=True)
    if nb < 50:
        print(f'  {s}: too few bound frames to cluster'); return None

    eps_all = bound_episodes(mask, offs)
    idx = np.concatenate([np.arange(a, b) for a, b in eps_all])
    # Distance -> smooth contact score. Raw distances make every unengaged pair
    # sit at whatever cap is chosen, so the clustering separates "how far" rather
    # than "which contact", and the state profile saturates. A sigmoid centred on
    # the 1.0 nm contact cutoff gives ~1 when engaged and ~0 when not, so a state
    # profile reads directly as per-domain-pair contact probability.
    Xb = 1.0 / (1.0 + np.exp((X[idx] - BOUND_NM) / 0.2))
    keep = Xb.std(0) > 1e-4          # pairs that never vary carry no information
    Xb = Xb[:, keep]
    kept_labels = [l for l, kp in zip(labels, keep) if kp]
    mu, sd = Xb.mean(0), Xb.std(0) + 1e-6
    Z = (Xb - mu) / sd

    # remap episodes onto the compacted bound array for TICA's lagged pairs
    pos, eps_b, c = {v: i for i, v in enumerate(idx)}, [], 0
    for a, b in eps_all:
        eps_b.append((c, c + (b - a))); c += (b - a)

    lag = max(1, int(round(lag_ns / (NS_PER_FRAME * stride))))
    used = reduce_mode
    if reduce_mode == 'tica':
        r = tica(Z, eps_b, lag, n=2)
        if r is None:
            print(f'  {s}: too few within-episode lagged pairs at lag {lag_ns} ns -> PCA')
            Y, comps, expl = pca(Z, 2); used = 'pca (tica fell back)'
        else:
            Y, comps, expl = r
    else:
        Y, comps, expl = pca(Z, 2)

    ks = range(2, kmax + 1)
    sils = {}
    for k in ks:
        lab, _ = kmeans(Y, k, seed=0)
        sils[k] = silhouette(Y, lab)
    best = fixed_k or max(sils, key=lambda k: sils[k])
    lab, C = kmeans(Y, best, seed=0)
    print(f'  {s}: {used}, k={best} (silhouette {sils[best]:.3f})', flush=True)

    # representative frame per state = nearest to the centroid
    reps = []
    for i in range(best):
        m = lab == i
        if not m.any():
            continue
        j = np.where(m)[0][np.argmin(((Y[m] - C[i])**2).sum(1))]
        src = srcs[idx[j]]
        pdb = STATES / f'{s}_state{i+1}.pdb'
        try:
            write_state_pdb(s, src, pdb)
        except Exception as e:
            print(f'    ! state {i+1} pdb: {e}')
        reps.append(dict(state=i+1, pop=float(m.mean()), n=int(m.sum()),
                         rep=src[0], frame=int(src[1]), pdb=pdb.name))
    return dict(set=s, Y=Y, lab=lab, best=best, sils=sils, used=used,
                Xb=Xb, labels=kept_labels, reps=reps, expl=expl,
                bound_pct=100*nb/len(X))


def plot_set(r):
    s = r['set']
    fig = plt.figure(figsize=(13.5, 4.2))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.5], wspace=0.35)

    ax = fig.add_subplot(gs[0, 0])
    for i in range(r['best']):
        m = r['lab'] == i
        ax.scatter(r['Y'][m, 0], r['Y'][m, 1], s=5, alpha=0.55,
                   color=STATE_COLORS[i % len(STATE_COLORS)],
                   label=f'S{i+1} ({100*m.mean():.0f}%)')
    ax.set_xlabel(f"{r['used'].split()[0].upper()} 1"); ax.set_ylabel(f"{r['used'].split()[0].upper()} 2")
    ax.set_title(f'{s} — bound-frame states', fontsize=9, loc='left')
    ax.legend(fontsize=6, frameon=False, markerscale=2)
    ax.spines[['top', 'right']].set_visible(False)

    ax = fig.add_subplot(gs[0, 1])
    ks = sorted(r['sils']); ax.plot(ks, [r['sils'][k] for k in ks], 'o-', color='#4A3AA7')
    ax.axvline(r['best'], color='#C8481A', ls='--', lw=1)
    ax.set_xlabel('k'); ax.set_ylabel('silhouette')
    ax.set_title('cluster-count sweep', fontsize=9, loc='left')
    ax.spines[['top', 'right']].set_visible(False)

    # per-state interface profile: the domain pairs that distinguish the states
    ax = fig.add_subplot(gs[0, 2])
    M = np.array([r['Xb'][r['lab'] == i].mean(0) for i in range(r['best'])])
    spread = M.max(0) - M.min(0)
    top = np.argsort(spread)[::-1][:12]
    im = ax.imshow(M[:, top], aspect='auto', cmap='magma_r', vmin=0, vmax=1)
    ax.set_xticks(range(len(top)))
    ax.set_xticklabels([r['labels'][t] for t in top], rotation=90, fontsize=5)
    ax.set_yticks(range(r['best']))
    ax.set_yticklabels([f'S{i+1}' for i in range(r['best'])], fontsize=7)
    cb = plt.colorbar(im, ax=ax, fraction=0.046); cb.set_label('contact probability', fontsize=6)
    cb.ax.tick_params(labelsize=5)
    ax.set_title('domain pairs that distinguish the states', fontsize=9, loc='left')
    fig.savefig(FIG / f'10_{s}_states.png', dpi=150, bbox_inches='tight')
    plt.close(fig)


def main():
    ap = ArgumentParser()
    ap.add_argument('sets', nargs='*')
    ap.add_argument('--reduce', default='pca', choices=['pca', 'tica'])
    ap.add_argument('--kmax', type=int, default=6)
    ap.add_argument('--k', type=int, default=None)
    ap.add_argument('--stride', type=int, default=DEFAULT_STRIDE)
    ap.add_argument('--lag-ns', type=float, default=1.0)
    a = ap.parse_args()
    todo = a.sets or [s for s in SETS if (ROOT / s).is_dir()]
    import pandas as pd
    rows = []
    for s in todo:
        r = run_set(s, a.reduce, a.kmax, a.stride, a.lag_ns, a.k)
        if r is None:
            continue
        plot_set(r)
        for d in r['reps']:
            rows.append(dict(set=s, reduce=r['used'], silhouette=round(r['sils'][r['best']], 3),
                             bound_pct=round(r['bound_pct'], 1), **d))
    if rows:
        df = pd.DataFrame(rows)
        df.to_csv(OUT / 'representative_states.csv', index=False)
        print(df.to_string(index=False))
        print(f"\nstructures in {STATES}")


if __name__ == '__main__':
    main()
