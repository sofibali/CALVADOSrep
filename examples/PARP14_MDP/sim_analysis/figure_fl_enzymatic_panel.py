#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
How do PARP14's enzymatic domains behave in the full-length simulation?

Several analyses converge on MD1's active site being the least reachable
(Q3: lowest solid-angle accessibility in every construct; Q4: writer-leaning
dominance everywhere; Q7: ART keeps its own catalytic face free). This script
looks at the FULL-LENGTH ensemble specifically and asks what the enzymatic
domains are actually doing to each other.

THE ORIENTATION IDEA
--------------------
Each enzymatic domain gets an active-site VECTOR:

    v_A = unit(catalytic_centroid_A - restrained_COM_A)

i.e. which way the pocket faces from the middle of the folded domain. Then for
an ordered pair A -> B:

    theta(A->B) = angle between v_A and unit(COM_B - COM_A)

    theta ~   0 deg  A's pocket points straight AT B      -> occluded by B
    theta ~  90 deg  A's pocket points past B             -> unaffected
    theta ~ 180 deg  A's pocket points away from B        -> free

This is the geometric content of "alternate access": a pocket can be buried not
because the protein is compact but because one particular partner sits in front
of it. theta is direction-only and so is independent of how far apart the
domains are -- distance is reported separately, and the occlusion score combines
the two.

FIGURES (FL only)
-----------------
  1 site_exposure          per-site SAA and inter-domain contact load
  2 site_orientation       theta distributions, every ordered enzymatic pair
  3 occlusion_matrix       fraction of frames each pocket is pointed at a partner
  4 face_contact_matrix    which face of A meets which face of B
  5 defining_crosslinks    the Lys pairs that would report each interaction

Usage:
    python figure_fl_enzymatic_panel.py
    python figure_fl_enzymatic_panel.py --set fl_optimized --target-frames 2000
"""

import os
import sys
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

import os as _os_boot, sys as _sys_boot  # __ROOTBOOT__ (script lives in a subfolder; root=parent)
_sys_boot.path.insert(0, _os_boot.path.dirname(_os_boot.path.dirname(_os_boot.path.abspath(__file__))))
import sim_registry as reg

CWD = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_PATH = os.path.join(CWD, 'data')
sys.path.insert(0, CWD)
from _fig_layout import get_fig_dir as _get_fig_dir
from figure_face_contacts import restrained_domains, FACES, FACE_COLOR

# Enzymatic / ligand-binding domains, in sequence order.
ENZ = ['MD1', 'MD2', 'MD3', 'WWE', 'ART']
ROLE = {'MD1': 'eraser', 'MD2': 'reader', 'MD3': 'reader',
        'WWE': 'PAR reader', 'ART': 'writer'}
SITE_COLOR = reg.SITE_COLORS
# sasa_face_per_residue.csv names MD1's block MD1L1; everything else matches.
FACE_DOMAIN = {'MD1': 'MD1L1', 'MD2': 'MD2', 'MD3': 'MD3', 'ART': 'ART'}

CONTACT_CUTOFF = 1.0     # nm, CA-CA
POINTS_AT_DEG = 60.0     # theta below this = pocket pointed at the partner
NEAR_NM = 6.0            # and partner COM within this = actually occluding
SKIP_FRAMES = 50


def _resolve(set_key, seed, sample):
    sysname = reg.SETS[set_key]['sysname']
    d = os.path.join(CWD, set_key, f'seed-{seed}_sample-{sample}')
    dcd = os.path.join(d, f'{sysname}.dcd')
    pdb = next((os.path.join(d, n) for n in ('top.pdb', 'checkpoint.pdb', 'restart.pdb')
                if os.path.isfile(os.path.join(d, n))), os.path.join(d, 'top.pdb'))
    return pdb, dcd


def load_faces():
    """FL resid -> (face_domain, face)."""
    path = os.path.join(DATA_PATH, 'sasa_face_per_residue.csv')
    if not os.path.isfile(path):
        sys.exit("ERROR: data/sasa_face_per_residue.csv missing. "
                 "Run figure_sasa_faces.py first.")
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            if r['face']:
                out[int(r['resid'])] = (r['domain'], r['face'])
    if not any(v[1] == 'rim' for v in out.values()):
        sys.exit("ERROR: face CSV predates the three-face split -- rerun "
                 "figure_sasa_faces.py.")
    return out


def worker(args):
    set_key, seed, sample, target_frames = args
    import MDAnalysis as mda
    import warnings
    warnings.filterwarnings('ignore')
    from scipy.spatial import cKDTree

    pdb, dcd = _resolve(set_key, seed, sample)
    if not os.path.isfile(dcd):
        return None
    faces = load_faces()
    rd = restrained_domains()

    l2f = reg.build_construct_to_fl_map(reg.get_units(set_key), set_key=set_key)
    u = mda.Universe(pdb, dcd)
    fl = np.array([l2f(int(r)) or -1 for r in u.atoms.resids])

    # restrained-domain index per bead (-1 = flexible)
    names = list(rd)
    dom = np.full(len(fl), -1, int)
    for i, n in enumerate(names):
        lo, hi = rd[n]
        dom[(fl >= lo) & (fl <= hi)] = i

    # per-enzymatic-domain masks: whole restrained body, and catalytic residues
    body, cat = {}, {}
    for d in ENZ:
        key = 'md1l1' if d == 'MD1' else d.lower()
        key = {'md1l1': 'md1l1', 'md2': 'md2', 'md3': 'md3',
               'wwe': 'wwe', 'art': 'art'}.get(key, key)
        if key not in rd:
            continue
        lo, hi = rd[key]
        b = (fl >= lo) & (fl <= hi)
        c = np.isin(fl, reg.ACTIVE_SITES_FL[d]['catalytic'])
        if b.sum() and c.sum():
            body[d], cat[d] = b, c
    present = [d for d in ENZ if d in body]
    if len(present) < 2:
        return None

    # face membership per bead, for the enzymatic domains that have faces
    face_idx = {}
    for d, fd in FACE_DOMAIN.items():
        if d not in body:
            continue
        for f in FACES:
            face_idx[(d, f)] = np.array(
                [faces.get(int(x), (None, None))[1] == f
                 and faces.get(int(x), (None, None))[0] == fd for x in fl])

    step = 1
    if target_frames:
        usable = max(0, len(u.trajectory) - SKIP_FRAMES)
        step = max(1, usable // target_frames)

    theta = defaultdict(list)      # (A,B) -> angle A's site vector vs A->B
    dist = defaultdict(list)       # (A,B) -> COM-COM nm
    sitesite = defaultdict(list)   # (A,B) -> angle between the two site vectors
    facepair = defaultdict(float)  # (A,fa,B,fb) -> contact count
    n_frames = 0

    for ts in u.trajectory[SKIP_FRAMES::step]:
        pos = u.atoms.positions / 10.0
        com = {d: pos[body[d]].mean(0) for d in present}
        vec = {}
        for d in present:
            v = pos[cat[d]].mean(0) - com[d]
            n = np.linalg.norm(v)
            if n > 1e-6:
                vec[d] = v / n
        for a in present:
            if a not in vec:
                continue
            for b in present:
                if b == a:
                    continue
                w = com[b] - com[a]
                dd = np.linalg.norm(w)
                if dd < 1e-6:
                    continue
                theta[(a, b)].append(
                    np.degrees(np.arccos(np.clip(np.dot(vec[a], w / dd), -1, 1))))
                dist[(a, b)].append(dd)
                if b in vec and a < b:
                    sitesite[(a, b)].append(
                        np.degrees(np.arccos(np.clip(np.dot(vec[a], vec[b]), -1, 1))))

        # face-resolved inter-domain contacts
        pairs = cKDTree(pos).query_pairs(CONTACT_CUTOFF, output_type='ndarray')
        if len(pairs):
            p, q = pairs[:, 0], pairs[:, 1]
            m = (dom[p] >= 0) & (dom[q] >= 0) & (dom[p] != dom[q])
            p, q = p[m], q[m]
            for (i, j) in ((p, q), (q, p)):
                for d in FACE_DOMAIN:
                    if d not in body:
                        continue
                    for fa in FACES:
                        sel = face_idx[(d, fa)][i]
                        if not sel.any():
                            continue
                        js = j[sel]
                        for d2 in FACE_DOMAIN:
                            if d2 == d or d2 not in body:
                                continue
                            for fb in FACES:
                                facepair[(d, fa, d2, fb)] += int(
                                    face_idx[(d2, fb)][js].sum())
        n_frames += 1

    if not n_frames:
        return None
    return dict(
        theta={k: (float(np.mean(v)), float(np.median(v)),
                   float(np.mean(np.array(v) < POINTS_AT_DEG)))
               for k, v in theta.items()},
        theta_raw={k: np.array(v, dtype=np.float32) for k, v in theta.items()},
        dist={k: float(np.mean(v)) for k, v in dist.items()},
        occl={k: float(np.mean((np.array(theta[k]) < POINTS_AT_DEG)
                               & (np.array(dist[k]) < NEAR_NM)))
              for k in theta},
        sitesite={k: float(np.mean(v)) for k, v in sitesite.items()},
        facepair={k: v / n_frames for k, v in facepair.items()},
        present=present, n_frames=n_frames)


# ============================================================
# Figures
# ============================================================

def fig_site_exposure(set_key, res, figdir):
    """Per-site SAA (cached) alongside inter-domain contact load."""
    npz = os.path.join(DATA_PATH, 'accessibility_stats.npz')
    saa = {}
    if os.path.isfile(npz):
        with np.load(npz, allow_pickle=True) as d:
            for s in ENZ:
                k = f'{set_key}_{s}_saa'
                if k in d.files:
                    saa[s] = np.asarray(d[k])
    sites = [s for s in ENZ if s in saa]
    if not sites:
        print("  (no cached SAA for this set -- skipping exposure figure)")
        return
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.3))
    parts = ax1.violinplot([saa[s] for s in sites], showmeans=True, widths=0.75)
    for pc, s in zip(parts['bodies'], sites):
        pc.set_facecolor(SITE_COLOR[s]); pc.set_alpha(0.65)
    ax1.set_xticks(range(1, len(sites) + 1))
    ax1.set_xticklabels([f'{s}\n{ROLE[s]}' for s in sites], fontsize=9)
    ax1.set_ylabel('solid-angle accessibility')
    ax1.set_title('How open is each pocket?', fontsize=11)
    ax1.grid(axis='y', alpha=0.25)

    load = {s: sum(v for (a, fa, b, fb), v in res['facepair'].items() if a == s)
            for s in sites}
    ax2.bar(range(len(sites)), [load[s] for s in sites],
            color=[SITE_COLOR[s] for s in sites], edgecolor='black', linewidth=0.5)
    ax2.set_xticks(range(len(sites)))
    ax2.set_xticklabels(sites, fontsize=9)
    ax2.set_ylabel('inter-domain contacts per frame')
    ax2.set_title('How much does each domain touch the others?', fontsize=11)
    ax2.grid(axis='y', alpha=0.25)
    fig.suptitle(f'{set_key} — enzymatic site exposure', fontsize=12)
    _save(fig, figdir, '1_site_exposure')


def fig_orientation(res, figdir):
    """theta distributions: does A's pocket point at B?"""
    pairs = sorted(res['theta_raw'])
    doms = sorted({a for a, b in pairs}, key=ENZ.index)
    fig, axes = plt.subplots(1, len(doms), figsize=(3.0 * len(doms), 4.2),
                             sharey=True)
    if len(doms) == 1:
        axes = [axes]
    for ax, a in zip(axes, doms):
        partners = [b for (x, b) in pairs if x == a]
        data = [res['theta_raw'][(a, b)] for b in partners]
        parts = ax.violinplot(data, showmedians=True, widths=0.8)
        for pc, b in zip(parts['bodies'], partners):
            pc.set_facecolor(SITE_COLOR[b]); pc.set_alpha(0.6)
        ax.axhline(POINTS_AT_DEG, ls='--', lw=1, color='#c0392b')
        ax.axhline(90, ls=':', lw=1, color='#555')
        ax.set_xticks(range(1, len(partners) + 1))
        ax.set_xticklabels(partners, fontsize=8.5, rotation=45)
        ax.set_title(f'{a} pocket\n({ROLE[a]})', fontsize=10)
        ax.set_ylim(0, 180); ax.set_yticks([0, 60, 90, 120, 180])
        ax.grid(axis='y', alpha=0.2)
    axes[0].set_ylabel('angle between pocket vector and direction to partner (°)')
    fig.suptitle('Which way does each pocket face?\n'
                 f'below the red line ({POINTS_AT_DEG:.0f}°) the pocket points AT that partner; '
                 '180° = points away', fontsize=11)
    fig.tight_layout()
    _save(fig, figdir, '2_site_orientation', tight=False)


def fig_occlusion(res, figdir):
    """Fraction of frames each pocket is pointed at AND close to a partner."""
    doms = sorted({a for a, b in res['occl']}, key=ENZ.index)
    mat = np.full((len(doms), len(doms)), np.nan)
    for i, a in enumerate(doms):
        for j, b in enumerate(doms):
            if (a, b) in res['occl']:
                mat[i, j] = res['occl'][(a, b)]
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    cmap = plt.get_cmap('magma_r').copy(); cmap.set_bad('#e8e8e8')
    im = ax.imshow(np.ma.masked_invalid(mat), cmap=cmap, vmin=0)
    ax.set_xticks(range(len(doms))); ax.set_xticklabels(doms, fontsize=10)
    ax.set_yticks(range(len(doms))); ax.set_yticklabels(
        [f'{d}\n{ROLE[d]}' for d in doms], fontsize=9)
    ax.set_xlabel('occluding partner'); ax.set_ylabel('pocket')
    for i in range(len(doms)):
        for j in range(len(doms)):
            if np.isnan(mat[i, j]):
                ax.text(j, i, '—', ha='center', va='center', color='#999')
            else:
                ax.text(j, i, f'{mat[i, j]:.2f}', ha='center', va='center',
                        fontsize=9,
                        color='white' if mat[i, j] > np.nanmax(mat) * 0.55 else 'black')
    plt.colorbar(im, ax=ax, shrink=0.8,
                 label=f'fraction of frames (<{POINTS_AT_DEG:.0f}° and <{NEAR_NM:.0f} nm)')
    ax.set_title('Who blocks whose pocket?\n'
                 'row pocket points at column partner AND the partner is close',
                 fontsize=11)
    _save(fig, figdir, '3_occlusion_matrix')


def fig_face_matrix(res, figdir):
    """Which face of A meets which face of B."""
    fp = res['facepair']
    doms = sorted({a for a, fa, b, fb in fp}, key=lambda d: ENZ.index(d)
                  if d in ENZ else 99)
    pairs = [(a, b) for i, a in enumerate(doms) for b in doms[i + 1:]]
    pairs = [p for p in pairs
             if any(fp.get((p[0], fa, p[1], fb), 0) for fa in FACES for fb in FACES)]
    if not pairs:
        print("  (no face-resolved contacts -- skipping face matrix)")
        return
    n = len(pairs)
    ncol = min(4, n); nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.0 * ncol, 3.1 * nrow),
                             squeeze=False)
    vmax = max(fp.values()) or 1
    for k, (a, b) in enumerate(pairs):
        ax = axes[k // ncol][k % ncol]
        m = np.array([[fp.get((a, fa, b, fb), 0.0) for fb in FACES]
                      for fa in FACES])
        im = ax.imshow(m, cmap='viridis', vmin=0, vmax=vmax)
        ax.set_xticks(range(3)); ax.set_xticklabels(FACES, fontsize=8, rotation=45)
        ax.set_yticks(range(3)); ax.set_yticklabels(FACES, fontsize=8)
        ax.set_xlabel(b, fontsize=9); ax.set_ylabel(a, fontsize=9)
        for i in range(3):
            for j in range(3):
                if m[i, j] > 0:
                    ax.text(j, i, f'{m[i, j]:.0f}', ha='center', va='center',
                            fontsize=7.5,
                            color='white' if m[i, j] < vmax * 0.55 else 'black')
    for k in range(n, nrow * ncol):
        axes[k // ncol][k % ncol].axis('off')
    fig.suptitle('Face-to-face contacts between enzymatic domains\n'
                 '(contacts per frame; row = face of the left domain)', fontsize=11)
    fig.tight_layout()
    _save(fig, figdir, '4_face_contact_matrix', tight=False)


def fig_crosslinks(set_key, figdir):
    """The Lys pairs that would report each enzymatic-domain interaction."""
    path = os.path.join(DATA_PATH, f'{set_key}_lys_sasd.csv')
    if not os.path.isfile(path):
        print(f"  (no {set_key}_lys_sasd.csv -- skipping crosslink figure)")
        return
    faces = load_faces()
    rows = []
    with open(path) as fh:
        for r in csv.DictReader(fh):
            if r.get('inter_domain') != '1' or r['reachable'] != '1':
                continue
            di, dj = r['domain_i'], r['domain_j']
            if di not in ENZ or dj not in ENZ:
                continue
            fi = faces.get(int(r['resid_i']), (None, None))[1]
            fj = faces.get(int(r['resid_j']), (None, None))[1]
            rows.append((tuple(sorted((di, dj))), int(r['resid_i']),
                         int(r['resid_j']), float(r['euclid_persistence']),
                         float(r['sasd_min_nm']), fi, fj))
    if not rows:
        print("  (no enzymatic-domain crosslinks -- skipping)")
        return
    bypair = defaultdict(list)
    for p, i, j, per, sasd, fi, fj in rows:
        bypair[p].append((per, sasd, i, j, fi, fj))
    order = sorted(bypair, key=lambda p: -len(bypair[p]))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 4.6),
                                   gridspec_kw={'width_ratios': [1, 1.35]})
    lbl = [f'{a}-{b}' for a, b in order]
    ax1.barh(range(len(order)), [len(bypair[p]) for p in order],
             color='#4363d8', edgecolor='black', linewidth=0.5)
    ax1.set_yticks(range(len(order))); ax1.set_yticklabels(lbl, fontsize=9)
    ax1.invert_yaxis()
    ax1.set_xlabel('reachable crosslinks')
    ax1.set_title('Crosslinks per enzymatic domain pair', fontsize=11)
    ax1.grid(axis='x', alpha=0.25)

    top = []
    for p in order:
        for per, sasd, i, j, fi, fj in sorted(bypair[p], key=lambda t: (-t[0], t[1]))[:3]:
            top.append((f'{p[0]}-{p[1]}', f'K{i}–K{j}', per, sasd, fi, fj))
    top = top[:14]
    y = np.arange(len(top))
    ax2.barh(y, [t[2] for t in top], color='#2ca02c', edgecolor='black',
             linewidth=0.5, height=0.62)
    ax2.set_yticks(y)
    ax2.set_yticklabels([f'{t[0]}  {t[1]}' for t in top], fontsize=8)
    ax2.invert_yaxis()
    ax2.set_xlabel('persistence (fraction of frames within DSS reach)')
    ax2.set_xlim(0, 1.05)
    for yi, t in zip(y, top):
        ax2.text(1.02, yi, f'SASD {t[3]:.2f} nm · {t[4] or "?"}→{t[5] or "?"}',
                 va='center', fontsize=7, color='#555')
    ax2.set_title('Best reporters, with the faces they connect', fontsize=11)
    fig.suptitle(f'{set_key} — crosslinks that would define these interactions',
                 fontsize=12)
    fig.tight_layout()
    _save(fig, figdir, '5_defining_crosslinks', tight=False)

    out = os.path.join(DATA_PATH, f'{set_key}_enzymatic_crosslinks.csv')
    with open(out, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['domain_pair', 'resid_i', 'resid_j', 'face_i', 'face_j',
                    'persistence', 'sasd_min_nm'])
        for p, i, j, per, sasd, fi, fj in sorted(rows, key=lambda t: (t[0], -t[3])):
            w.writerow([f'{p[0]}-{p[1]}', i, j, fi or '', fj or '',
                        f'{per:.4f}', f'{sasd:.3f}'])
    print(f"  Saved: {os.path.relpath(out, CWD)}")


def _save(fig, figdir, name, tight=True):
    if tight:
        fig.tight_layout()
    for ext in ('png', 'svg'):
        fig.savefig(os.path.join(figdir, f'{name}.{ext}'), dpi=150,
                    bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {name}.png/.svg")


def merge(results):
    """Pool replicate results."""
    out = dict(theta_raw=defaultdict(list), occl=defaultdict(list),
               dist=defaultdict(list), sitesite=defaultdict(list),
               facepair=defaultdict(list))
    for r in results:
        for k, v in r['theta_raw'].items():
            out['theta_raw'][k].append(v)
        for f in ('occl', 'dist', 'sitesite'):
            for k, v in r[f].items():
                out[f][k].append(v)
        for k, v in r['facepair'].items():
            out['facepair'][k].append(v)
    return dict(
        theta_raw={k: np.concatenate(v) for k, v in out['theta_raw'].items()},
        occl={k: float(np.mean(v)) for k, v in out['occl'].items()},
        dist={k: float(np.mean(v)) for k, v in out['dist'].items()},
        sitesite={k: float(np.mean(v)) for k, v in out['sitesite'].items()},
        facepair={k: float(np.mean(v)) for k, v in out['facepair'].items()})


def main():
    ap = ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--set', default='fl')
    ap.add_argument('--target-frames', type=int, default=1000)
    ap.add_argument('--workers', type=int, default=14)
    args = ap.parse_args()

    jobs = [(args.set, s, m, args.target_frames)
            for s in reg.SEEDS for m in reg.SAMPLES
            if os.path.isfile(_resolve(args.set, s, m)[1])]
    if not jobs:
        sys.exit(f"ERROR: no trajectories for '{args.set}'")
    print(f"  {args.set}: {len(jobs)} replicates, ~{args.target_frames} frames each")

    results = []
    with ProcessPoolExecutor(max_workers=min(args.workers, len(jobs))) as pool:
        for r in pool.map(worker, jobs):
            if r:
                results.append(r)
    if not results:
        sys.exit("ERROR: no replicate produced results")
    res = merge(results)
    print(f"  pooled {len(results)} replicates")

    figdir = str(_get_fig_dir('02_main_analysis', subname='enzymatic_panel',
                              sims=[args.set]))
    fig_site_exposure(args.set, res, figdir)
    fig_orientation(res, figdir)
    fig_occlusion(res, figdir)
    fig_face_matrix(res, figdir)
    fig_crosslinks(args.set, figdir)

    print(f"\n{'=' * 82}\nPOCKET ORIENTATION — {args.set}\n{'=' * 82}")
    print(f"{'pocket':8s}{'partner':9s}{'median θ':>10}{'mean d (nm)':>13}"
          f"{'occluded':>11}")
    print('-' * 82)
    for (a, b), th in sorted(res['theta_raw'].items(),
                             key=lambda kv: -res['occl'][kv[0]]):
        print(f"{a:8s}{b:9s}{np.median(th):>10.1f}{res['dist'][(a, b)]:>13.2f}"
              f"{res['occl'][(a, b)]:>11.3f}")
    print(f"\n  'occluded' = fraction of frames the pocket points within "
          f"{POINTS_AT_DEG:.0f}° of a partner that is also within {NEAR_NM:.0f} nm.")
    print(f"  Figures -> {os.path.relpath(figdir, CWD)}")


if __name__ == '__main__':
    main()
