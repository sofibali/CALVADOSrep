#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Why is MD1's face preference context-dependent?

Q7 found that ART and MD3 consistently present their BACK faces to the rest of
the protein, and that the rim is depleted everywhere. MD1 is the exception: it
prefers its ACTIVE face in fl / core / norrm / kh1_art_full and its BACK face in
mka_full / md_full.

That split tracked KH7a presence exactly (4/4 vs 2/2), which looked like a
mechanism -- especially since Q4 showed adding KH7a drops MD1 accessibility. A
per-partner breakdown on `core` refuted it: KH7a contributes only ~12% of MD1's
active-face contacts, behind MD3, ART and WWE.

This script does that breakdown for EVERY construct, so the question gets a real
answer rather than a single-construct anecdote. Three outcomes are possible:

  1. One partner dominates in the active-preferring constructs and is absent or
     minor in the back-preferring ones -> that partner is the mechanism.
  2. No single partner does, but the active-preferring constructs simply have
     MORE partners able to reach MD1 -> collective packing, not a specific
     interaction.
  3. The split does not correspond to partner composition at all -> the face
     preference is driven by something else (chain-length entropy, linker
     geometry) and should be reported as unexplained.

Usage:
    python diagnose_md1_faces.py
    python diagnose_md1_faces.py --domain MD3 --target-frames 600
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
from figure_face_contacts import restrained_domains, FACES, FACE_COLOR, _resolve

CONTACT_CUTOFF = 1.0
SKIP_FRAMES = 50
DEFAULT_SETS = ['fl', 'core', 'norrm', 'kh1_art_full', 'mka_full', 'md_full']
# sasa_face_per_residue.csv calls MD1's block MD1L1
FACE_NAME = {'MD1': 'MD1L1', 'MD2': 'MD2', 'MD3': 'MD3', 'ART': 'ART'}


def load_faces(face_domain):
    path = os.path.join(DATA_PATH, 'sasa_face_per_residue.csv')
    if not os.path.isfile(path):
        sys.exit("ERROR: run figure_sasa_faces.py first")
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            if r['domain'] == face_domain and r['face']:
                out[int(r['resid'])] = r['face']
    if not any(v == 'rim' for v in out.values()):
        sys.exit("ERROR: face CSV predates the three-face split")
    return out


def worker(args):
    """Face-resolved contacts on one domain, broken down BY PARTNER."""
    set_key, seed, sample, domain, target_frames = args
    import MDAnalysis as mda
    import warnings
    warnings.filterwarnings('ignore')
    from scipy.spatial import cKDTree

    pdb, dcd = _resolve(set_key, seed, sample)
    if not os.path.isfile(dcd):
        return None
    faces = load_faces(FACE_NAME[domain])
    rd = restrained_domains()
    unit = {'MD1': 'md1l1'}.get(domain, domain.lower())
    if unit not in rd:
        return None

    l2f = reg.build_construct_to_fl_map(reg.get_units(set_key), set_key=set_key)
    u = mda.Universe(pdb, dcd)
    fl = np.array([l2f(int(r)) or -1 for r in u.atoms.resids])

    names = list(rd)
    dom = np.full(len(fl), -1, int)
    for i, n in enumerate(names):
        lo, hi = rd[n]
        dom[(fl >= lo) & (fl <= hi)] = i
    me = names.index(unit)

    face_of = np.array([faces.get(int(x), '') for x in fl], dtype=object)

    step = 1
    if target_frames:
        usable = max(0, len(u.trajectory) - SKIP_FRAMES)
        step = max(1, usable // target_frames)

    acc = defaultdict(float)   # (face, partner) -> contacts
    nf = 0
    for ts in u.trajectory[SKIP_FRAMES::step]:
        pos = u.atoms.positions / 10.0
        pr = cKDTree(pos).query_pairs(CONTACT_CUTOFF, output_type='ndarray')
        if len(pr):
            a, b = pr[:, 0], pr[:, 1]
            for i, j in ((a, b), (b, a)):
                m = (dom[i] == me) & (dom[j] >= 0) & (dom[j] != me)
                for ii, jj in zip(i[m], j[m]):
                    f = face_of[ii]
                    if f:
                        acc[(f, names[dom[jj]])] += 1
        nf += 1
    if not nf:
        return None
    return {k: v / nf for k, v in acc.items()}


def main():
    ap = ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--domain', default='MD1', choices=list(FACE_NAME))
    ap.add_argument('--sets', nargs='+', default=DEFAULT_SETS)
    ap.add_argument('--target-frames', type=int, default=400)
    ap.add_argument('--workers', type=int, default=14)
    args = ap.parse_args()

    faces = load_faces(FACE_NAME[args.domain])
    size = {f: sum(1 for v in faces.values() if v == f) for f in FACES}
    print(f"  {args.domain} face sizes: " +
          "  ".join(f"{f}={size[f]}" for f in FACES))

    results = {}
    for s in args.sets:
        if s not in reg.SETS:
            print(f"  WARN: unknown set {s}")
            continue
        jobs = [(s, sd, sm, args.domain, args.target_frames)
                for sd in reg.SEEDS for sm in reg.SAMPLES
                if os.path.isfile(_resolve(s, sd, sm)[1])]
        if not jobs:
            continue
        print(f"  {s} ({len(jobs)} reps) ...", flush=True)
        acc = defaultdict(list)
        with ProcessPoolExecutor(max_workers=min(args.workers, len(jobs))) as pool:
            for r in pool.map(worker, jobs):
                if r:
                    for k, v in r.items():
                        acc[k].append(v)
        if acc:
            results[s] = {k: float(np.mean(v)) for k, v in acc.items()}
    if not results:
        sys.exit("ERROR: no results")

    partners = sorted({p for r in results.values() for (_, p) in r})

    # enrichment of the active face, per construct (same definition as Q7)
    print(f"\n{'=' * 92}")
    print(f"{args.domain} FACE PREFERENCE AND WHO CONTACTS THE ACTIVE FACE")
    print('=' * 92)
    print(f"{'construct':16s}{'act.enr':>9}{'pref':>7}{'act/frame':>11}   "
          f"top partners on the ACTIVE face")
    print('-' * 92)
    summary = {}
    for s, r in results.items():
        tot_c = sum(r.values())
        tot_n = sum(size.values())
        enr = {}
        for f in FACES:
            c = sum(v for (ff, _), v in r.items() if ff == f)
            enr[f] = (c / tot_c) / (size[f] / tot_n) if size[f] and tot_c else np.nan
        pref = max(FACES, key=lambda f: enr[f])
        act = {p: v for (f, p), v in r.items() if f == 'active'}
        act_tot = sum(act.values()) or 1
        top = sorted(act.items(), key=lambda t: -t[1])[:3]
        summary[s] = dict(enr=enr, pref=pref, act=act, act_tot=act_tot)
        print(f"{s:16s}{enr['active']:>9.2f}{pref:>7}{act_tot:>11.1f}   "
              + ", ".join(f"{p} {100*v/act_tot:.0f}%" for p, v in top))

    # does any single partner separate the two groups?
    act_pref = [s for s in results if summary[s]['pref'] == 'active']
    back_pref = [s for s in results if summary[s]['pref'] == 'back']
    print(f"\n  active-preferring: {', '.join(act_pref) or '-'}")
    print(f"  back-preferring:   {', '.join(back_pref) or '-'}")

    if act_pref and back_pref:
        print(f"\n{'partner':12s}{'mean % of active-face contacts':>32}   verdict")
        print(f"{'':12s}{'active-pref':>16}{'back-pref':>16}")
        print('-' * 92)
        for p in partners:
            fa = [100 * summary[s]['act'].get(p, 0) / summary[s]['act_tot']
                  for s in act_pref]
            fb = [100 * summary[s]['act'].get(p, 0) / summary[s]['act_tot']
                  for s in back_pref]
            ma, mb = np.mean(fa), np.mean(fb)
            # a partner "explains" the split only if it is substantial in one
            # group and near-absent in the other
            if ma >= 20 and mb <= 5:
                v = 'EXPLAINS the split'
            elif abs(ma - mb) >= 20:
                v = 'differs strongly'
            else:
                v = ''
            print(f"{p:12s}{ma:>16.1f}{mb:>16.1f}   {v}")

        n_act = [len([p for p in partners if summary[s]['act'].get(p, 0) > 0])
                 for s in act_pref]
        n_back = [len([p for p in partners if summary[s]['act'].get(p, 0) > 0])
                  for s in back_pref]
        print(f"\n  partners able to reach the active face: "
              f"active-pref {np.mean(n_act):.1f}, back-pref {np.mean(n_back):.1f}")

    make_figure(args.domain, results, summary, partners)

    out = os.path.join(DATA_PATH, f'{args.domain.lower()}_face_partners.csv')
    with open(out, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['set', 'face', 'partner', 'contacts_per_frame'])
        for s, r in results.items():
            for (f, p), v in sorted(r.items()):
                w.writerow([s, f, p, f'{v:.4f}'])
    print(f"\n  Saved: {os.path.relpath(out, CWD)}")


def make_figure(domain, results, summary, partners):
    sets = list(results)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 4.6),
                                 gridspec_kw={'width_ratios': [1, 1.5]})
    x = np.arange(len(sets)); w = 0.26
    for i, f in enumerate(FACES):
        a1.bar(x + (i - 1) * w, [summary[s]['enr'][f] for s in sets], w,
               color=FACE_COLOR[f], label=f, edgecolor='black', linewidth=0.4)
    a1.axhline(1.0, color='black', lw=1.1)
    a1.set_xticks(x); a1.set_xticklabels(sets, rotation=45, ha='right', fontsize=8)
    a1.set_ylabel('contact enrichment')
    a1.set_title(f'{domain} face preference', fontsize=11)
    a1.legend(frameon=False, fontsize=8)

    bottom = np.zeros(len(sets))
    cmap = plt.get_cmap('tab20')
    for k, p in enumerate(partners):
        v = np.array([100 * summary[s]['act'].get(p, 0) / summary[s]['act_tot']
                      for s in sets])
        a2.bar(x, v, bottom=bottom, label=p, color=cmap(k % 20),
               edgecolor='black', linewidth=0.3)
        bottom += v
    a2.set_xticks(x); a2.set_xticklabels(sets, rotation=45, ha='right', fontsize=8)
    a2.set_ylabel(f'% of {domain} active-face contacts')
    a2.set_title(f'Who contacts {domain}\'s active face?', fontsize=11)
    a2.legend(frameon=False, fontsize=7, ncol=2, loc='upper right')
    fig.tight_layout()
    figdir = str(_get_fig_dir('01_static_FL', subname='face_contacts'))
    for ext in ('png', 'svg'):
        fig.savefig(os.path.join(figdir, f'{domain.lower()}_face_partners.{ext}'),
                    dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {domain.lower()}_face_partners.png/.svg")


if __name__ == '__main__':
    main()
