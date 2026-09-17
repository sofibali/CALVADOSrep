#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Which face of each domain makes the inter-domain contacts?

THE QUESTION
------------
figure_sasa_faces.py splits every active-site domain into three faces along the
axis from the domain centre to its catalytic site:

    active   the active-site pole
    rim      the equatorial band  (deliberately wide -- near the equator the
             sign of the projection carries no geometric information)
    back     the opposite pole

The rim is treated as a face in its own right precisely because it is the
lateral surface a domain presents to its neighbours while leaving both poles
free. This script tests that: across the MD ensemble, which face actually
carries the inter-domain contacts?

Two readings matter:

  - If contacts concentrate on the RIM, the architecture leaves both catalytic
    and back poles free, and domain packing does not occlude the active site.
  - If contacts land on the ACTIVE face, a neighbouring domain is sitting over
    the catalytic pocket -- a structural explanation for the low accessibility
    of MD1 found in Q3/Q4.

ENRICHMENT, NOT RAW COUNTS
--------------------------
The rim holds ~45-50% of residues, so it would carry the most contacts even
with no preference at all. Every comparison here is therefore an enrichment:

    enrichment(face) = (contacts on face / all contacts)
                       / (residues in face / all residues)

1.0 = exactly as many contacts as the face's size predicts. >1 = preferred.

Contacts are CA-CA within CONTACT_CUTOFF between residues of DIFFERENT domains,
matching the 1.0 nm convention used by figure_md_distances.py.

Usage:
    python figure_face_contacts.py                      # default construct set
    python figure_face_contacts.py --sets fl core norrm
    python figure_face_contacts.py --target-frames 500  # faster
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

CONTACT_CUTOFF = 1.0      # nm, CA-CA; matches figure_md_distances.py

# Folded domains are defined by the RESTRAINED residues, not the full unit
# ranges. Those residues are held rigid by harmonic restraints in the
# simulation; everything else is free to flex, and a flexible linker brushing
# past a domain is not a domain-domain interface. Using the restraint list
# keeps this analysis consistent with what the simulation actually holds fixed.
RESTRAINT_YAML = os.path.join(CWD, 'fl_optimized', 'input', 'domains.yaml')
SKIP_FRAMES = 50
DEFAULT_TARGET_FRAMES = 1000
DEFAULT_SETS = ['fl', 'core', 'norrm', 'kh1_art_full', 'mka_full', 'md_full']

FACES = ['active', 'rim', 'back']
FACE_COLOR = {'active': '#27ae60', 'rim': '#f0b323', 'back': '#7f8c8d'}
SITE_DOMAINS = ['MD1L1', 'MD2', 'MD3', 'ART']


def load_faces():
    """FL resid -> (domain, face, rsa). From figure_sasa_faces.py."""
    path = os.path.join(DATA_PATH, 'sasa_face_per_residue.csv')
    if not os.path.isfile(path):
        sys.exit(f"ERROR: {path} not found. Run figure_sasa_faces.py first.")
    with open(path) as fh:
        head = fh.readline().rstrip('\n').split(',')
    if 'face' not in head:
        sys.exit(f"ERROR: {path} has no 'face' column.")
    out = {}
    n_rim = 0
    with open(path) as fh:
        for r in csv.DictReader(fh):
            if not r['face'] or r['domain'] not in SITE_DOMAINS:
                continue
            rsa = float(r['rsa']) if r['rsa'] not in ('', 'None') else np.nan
            out[int(r['resid'])] = (r['domain'], r['face'], rsa)
            n_rim += r['face'] == 'rim'
    if n_rim == 0:
        sys.exit("ERROR: no residues classified as 'rim'. The face CSV predates "
                 "the three-face split -- rerun figure_sasa_faces.py.")
    return out


def load_restraint_domains():
    """canonical domain name -> (start, end) of its restrained block(s), FL numbering.

    domains.yaml is an unnamed list of restrained blocks, so each is matched to
    the canonical unit it overlaps. Sub-blocks of one unit (fl_optimized splits
    KH1-6 into six) are merged, so "inter-domain" keeps the same meaning it has
    everywhere else in the project.
    """
    try:
        import yaml
        blocks = list(yaml.safe_load(open(RESTRAINT_YAML)).values())[0]
    except Exception as e:
        print(f"  NOTE: restraint domains unavailable ({e}); using FL_DOMAINS")
        return dict(reg.FL_DOMAINS)
    out = {}
    for name, (lo, hi) in reg.DOMAIN_UNITS.items():
        ov = [b for b in blocks if b[1] >= lo and b[0] <= hi]
        if ov:
            out[name] = (min(b[0] for b in ov), max(b[1] for b in ov))
    return out


_RESTRAINED = None


def restrained_domains():
    global _RESTRAINED
    if _RESTRAINED is None:
        _RESTRAINED = load_restraint_domains()
    return _RESTRAINED


def _resolve(set_key, seed, sample):
    sysname = reg.SETS[set_key]['sysname']
    d = os.path.join(CWD, set_key, f'seed-{seed}_sample-{sample}')
    dcd = os.path.join(d, f'{sysname}.dcd')
    pdb = next((os.path.join(d, n) for n in ('top.pdb', 'checkpoint.pdb', 'restart.pdb')
                if os.path.isfile(os.path.join(d, n))), os.path.join(d, 'top.pdb'))
    return pdb, dcd


def worker(args):
    """Per-residue inter-domain contact counts for one replicate, FL numbering."""
    set_key, seed, sample, target_frames = args
    import MDAnalysis as mda
    import warnings
    warnings.filterwarnings('ignore')
    from scipy.spatial import cKDTree

    pdb, dcd = _resolve(set_key, seed, sample)
    if not os.path.isfile(dcd):
        return None

    local_to_fl = reg.build_construct_to_fl_map(reg.get_units(set_key),
                                                set_key=set_key)
    u = mda.Universe(pdb, dcd)
    local_resids = u.atoms.resids
    fl = np.array([local_to_fl(int(r)) or -1 for r in local_resids])

    # Domain index per bead, from the RESTRAINED blocks. -1 = flexible/linker,
    # excluded from the contact count entirely.
    rd = restrained_domains()
    dom = np.full(len(fl), -1, dtype=int)
    names = list(rd)
    for i, n in enumerate(names):
        lo, hi = rd[n]
        dom[(fl >= lo) & (fl <= hi)] = i

    step = 1
    if target_frames:
        usable = max(0, len(u.trajectory) - SKIP_FRAMES)
        step = max(1, usable // target_frames)

    counts = defaultdict(float)
    n_frames = 0
    for ts in u.trajectory[SKIP_FRAMES::step]:
        pos = u.atoms.positions / 10.0
        pairs = cKDTree(pos).query_pairs(CONTACT_CUTOFF, output_type='ndarray')
        if len(pairs):
            a, b = pairs[:, 0], pairs[:, 1]
            # inter-DOMAIN only: both in a domain, and different domains
            m = (dom[a] >= 0) & (dom[b] >= 0) & (dom[a] != dom[b])
            for i in a[m]:
                counts[int(fl[i])] += 1
            for j in b[m]:
                counts[int(fl[j])] += 1
        n_frames += 1
    if not n_frames:
        return None
    return {k: v / n_frames for k, v in counts.items()}


def analyse_set(set_key, faces, target_frames, workers):
    jobs = [(set_key, s, m, target_frames)
            for s in reg.SEEDS for m in reg.SAMPLES
            if os.path.isfile(_resolve(set_key, s, m)[1])]
    if not jobs:
        return None
    per_res = defaultdict(list)
    with ProcessPoolExecutor(max_workers=min(workers, len(jobs))) as pool:
        for res in pool.map(worker, jobs):
            if res:
                for k, v in res.items():
                    per_res[k].append(v)
    mean = {k: float(np.mean(v)) for k, v in per_res.items()}

    # aggregate per (domain, face)
    contacts = defaultdict(float)
    sizes = defaultdict(int)
    rsa_sum = defaultdict(list)
    for resid, (dname, face, rsa) in faces.items():
        # only count residues the construct actually contains
        if resid not in mean and resid not in per_res:
            pass
        sizes[(dname, face)] += 1
        contacts[(dname, face)] += mean.get(resid, 0.0)
        if np.isfinite(rsa):
            rsa_sum[(dname, face)].append(rsa)
    return dict(contacts=dict(contacts), sizes=dict(sizes),
                rsa={k: float(np.mean(v)) for k, v in rsa_sum.items()},
                n_rep=len(jobs))


def enrichment(res, dname):
    """Contacts-per-residue on each face, normalised to the domain mean."""
    tot_c = sum(res['contacts'].get((dname, f), 0.0) for f in FACES)
    tot_n = sum(res['sizes'].get((dname, f), 0) for f in FACES)
    if tot_c <= 0 or tot_n <= 0:
        return {f: np.nan for f in FACES}
    out = {}
    for f in FACES:
        c = res['contacts'].get((dname, f), 0.0)
        n = res['sizes'].get((dname, f), 0)
        out[f] = (c / tot_c) / (n / tot_n) if n else np.nan
    return out


# ============================================================
# Figures
# ============================================================

def fig_composition(faces, figdir):
    """How many residues sit on each face, per domain."""
    counts = defaultdict(int)
    for _, (d, f, _) in faces.items():
        counts[(d, f)] += 1
    fig, ax = plt.subplots(figsize=(7.5, 4.4))
    x = np.arange(len(SITE_DOMAINS))
    bottom = np.zeros(len(SITE_DOMAINS))
    for f in FACES:
        v = np.array([counts[(d, f)] for d in SITE_DOMAINS], float)
        ax.bar(x, v, bottom=bottom, color=FACE_COLOR[f], label=f,
               edgecolor='black', linewidth=0.4)
        for xi, (vi, bi) in enumerate(zip(v, bottom)):
            if vi > 6:
                ax.text(xi, bi + vi / 2, f'{int(vi)}', ha='center', va='center',
                        fontsize=8, color='white' if f != 'rim' else 'black')
        bottom += v
    ax.set_xticks(x); ax.set_xticklabels(SITE_DOMAINS)
    ax.set_ylabel('residues')
    ax.set_title('Face composition per domain\n'
                 'rim is wide by construction (~45-50%) — see enrichment for preference',
                 fontsize=10)
    ax.legend(frameon=False, fontsize=9)
    _save(fig, figdir, 'face_composition')


def fig_rsa_by_face(results, figdir):
    """Mean RSA per face — is the active pole more exposed than the back?"""
    any_res = next(iter(results.values()))
    fig, ax = plt.subplots(figsize=(7.5, 4.4))
    w = 0.26
    x = np.arange(len(SITE_DOMAINS))
    for i, f in enumerate(FACES):
        v = [any_res['rsa'].get((d, f), np.nan) for d in SITE_DOMAINS]
        ax.bar(x + (i - 1) * w, v, w, color=FACE_COLOR[f], label=f,
               edgecolor='black', linewidth=0.4)
    ax.axhline(0.20, ls='--', lw=1, color='#555')
    ax.text(len(SITE_DOMAINS) - 0.45, 0.205, 'surface threshold', fontsize=7.5,
            color='#555', ha='right')
    ax.set_xticks(x); ax.set_xticklabels(SITE_DOMAINS)
    ax.set_ylabel('mean RSA')
    ax.set_title('Solvent exposure by face (static AF2 structure)', fontsize=11)
    ax.legend(frameon=False, fontsize=9)
    _save(fig, figdir, 'face_rsa')


def fig_contact_enrichment(results, figdir):
    """The headline: which face carries the inter-domain contacts?"""
    sets = list(results)
    fig, axes = plt.subplots(1, len(SITE_DOMAINS),
                             figsize=(3.1 * len(SITE_DOMAINS), 4.6),
                             sharey=True)
    if len(SITE_DOMAINS) == 1:
        axes = [axes]
    for ax, dname in zip(axes, SITE_DOMAINS):
        x = np.arange(len(sets))
        w = 0.26
        for i, f in enumerate(FACES):
            v = [enrichment(results[s], dname)[f] for s in sets]
            ax.bar(x + (i - 1) * w, v, w, color=FACE_COLOR[f],
                   label=f if ax is axes[0] else None,
                   edgecolor='black', linewidth=0.4)
        ax.axhline(1.0, color='black', lw=1.1)
        ax.set_xticks(x)
        ax.set_xticklabels(sets, rotation=45, ha='right', fontsize=8)
        ax.set_title(dname, fontsize=11)
        ax.grid(axis='y', alpha=0.25)
    axes[0].set_ylabel('contact enrichment\n(1.0 = as expected from face size)')
    axes[0].legend(frameon=False, fontsize=9, loc='upper left')
    fig.suptitle('Which face makes the inter-domain contacts?\n'
                 'above 1.0 = that face is preferred;  '
                 'active-face enrichment means a neighbour sits over the pocket',
                 fontsize=11)
    fig.tight_layout()
    _save(fig, figdir, 'face_contact_enrichment', tight=False)


def fig_active_face_load(results, figdir):
    """Absolute contact load on the active face, per construct."""
    sets = list(results)
    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    x = np.arange(len(sets))
    w = 0.2
    for i, dname in enumerate(SITE_DOMAINS):
        v = []
        for s in sets:
            c = results[s]['contacts'].get((dname, 'active'), 0.0)
            n = results[s]['sizes'].get((dname, 'active'), 0)
            v.append(c / n if n else np.nan)
        ax.bar(x + (i - 1.5) * w, v, w, label=dname, edgecolor='black',
               linewidth=0.4)
    ax.set_xticks(x); ax.set_xticklabels(sets, rotation=45, ha='right', fontsize=9)
    ax.set_ylabel('inter-domain contacts per active-face residue')
    ax.set_title('Absolute contact load on the active-site face', fontsize=11)
    ax.legend(frameon=False, fontsize=9, ncol=4)
    ax.grid(axis='y', alpha=0.25)
    _save(fig, figdir, 'active_face_contact_load')


def _save(fig, figdir, name, tight=True):
    if tight:
        fig.tight_layout()
    for ext in ('png', 'svg'):
        fig.savefig(os.path.join(figdir, f'{name}.{ext}'), dpi=150,
                    bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {name}.png/.svg")


def main():
    ap = ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--sets', nargs='+', default=DEFAULT_SETS)
    ap.add_argument('--target-frames', type=int, default=DEFAULT_TARGET_FRAMES)
    ap.add_argument('--workers', type=int, default=12)
    args = ap.parse_args()

    faces = load_faces()
    nrim = sum(1 for v in faces.values() if v[1] == 'rim')
    print(f"  Face assignment: {len(faces)} residues across {len(SITE_DOMAINS)} "
          f"domains ({nrim} rim)")

    results = {}
    for s in args.sets:
        if s not in reg.SETS:
            print(f"  WARN: unknown set '{s}'")
            continue
        print(f"  {s} ...", flush=True)
        r = analyse_set(s, faces, args.target_frames, args.workers)
        if r:
            results[s] = r
            print(f"    {r['n_rep']} replicates")
        else:
            print(f"    no trajectories, skipped")
    if not results:
        sys.exit("ERROR: no sets had trajectories.")

    figdir = str(_get_fig_dir('01_static_FL', subname='face_contacts'))
    fig_composition(faces, figdir)
    fig_rsa_by_face(results, figdir)
    fig_contact_enrichment(results, figdir)
    fig_active_face_load(results, figdir)

    out = os.path.join(DATA_PATH, 'face_contacts.csv')
    with open(out, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['set', 'domain', 'face', 'n_residues', 'mean_rsa',
                    'contacts_total', 'contacts_per_residue', 'enrichment'])
        for s, r in results.items():
            for d in SITE_DOMAINS:
                e = enrichment(r, d)
                for f in FACES:
                    n = r['sizes'].get((d, f), 0)
                    c = r['contacts'].get((d, f), 0.0)
                    w.writerow([s, d, f, n,
                                f"{r['rsa'].get((d, f), float('nan')):.4f}",
                                f'{c:.3f}', f'{c / n:.4f}' if n else '',
                                f'{e[f]:.3f}' if np.isfinite(e[f]) else ''])
    print(f"  Saved: {os.path.relpath(out, CWD)}")

    # console summary
    print(f"\n{'=' * 78}\nCONTACT ENRICHMENT BY FACE  (1.0 = as expected from size)\n{'=' * 78}")
    print(f"{'set':16s}{'domain':8s}" + ''.join(f'{f:>10s}' for f in FACES)
          + '   preferred')
    print('-' * 78)
    for s, r in results.items():
        for d in SITE_DOMAINS:
            e = enrichment(r, d)
            if not np.isfinite(e['rim']):
                continue
            best = max(FACES, key=lambda f: (e[f] if np.isfinite(e[f]) else -1))
            print(f"{s:16s}{d:8s}" + ''.join(f'{e[f]:>10.2f}' for f in FACES)
                  + f'   {best}')
    print(f"\n  Figures -> {os.path.relpath(figdir, CWD)}")


if __name__ == '__main__':
    main()
